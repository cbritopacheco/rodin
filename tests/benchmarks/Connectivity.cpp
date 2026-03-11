/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>

#include <array>
#include <string>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Benchmarks
{
  class ConnectivityBenchmark : public benchmark::Fixture
  {
    public:
      struct MeshSpec
      {
        Polytope::Type type;
        std::string name;
        std::array<size_t, 3> resolution;
        size_t topologicalDimension;
      };

      static constexpr MeshSpec kMeshes[] =
      {
        { Polytope::Type::Segment,       "Edge",          { 4096, 1, 1 }, 1 },
        { Polytope::Type::Triangle,      "Triangle",      {   64, 64, 1 }, 2 },
        { Polytope::Type::Quadrilateral, "Quadrilateral", {   64, 64, 1 }, 2 },
        { Polytope::Type::Tetrahedron,   "Tetrahedron",   {    16,  16, 16 }, 3 },
        { Polytope::Type::Hexahedron,    "Hexahedron",    {    16,  16, 16 }, 3 },
        { Polytope::Type::Wedge,         "Wedge",         {    16,  16, 16 }, 3 }
      };

      static LocalMesh makeMesh(const MeshSpec& spec)
      {
        switch (spec.topologicalDimension)
        {
          case 1:
            return LocalMesh::UniformGrid(spec.type, { spec.resolution[0] });

          case 2:
            return LocalMesh::UniformGrid(spec.type, { spec.resolution[0], spec.resolution[1] });

          case 3:
            return LocalMesh::UniformGrid(spec.type, { spec.resolution[0], spec.resolution[1], spec.resolution[2] });

          default:
            assert(false);
            return LocalMesh();
        }
      }

      static void setCommonCounters(benchmark::State& st, const LocalMesh& mesh)
      {
        const auto& conn = mesh.getConnectivity();
        const size_t D = mesh.getDimension();

        st.counters["D"]        = static_cast<double>(D);
        st.counters["vertices"] = static_cast<double>(conn.getCount(0));
        st.counters["cells"]    = static_cast<double>(conn.getCount(D));

        if (D >= 1)
          st.counters["entities_1"] = static_cast<double>(conn.getCount(1));
        if (D >= 2)
          st.counters["entities_2"] = static_cast<double>(conn.getCount(2));
        if (D >= 3)
          st.counters["entities_3"] = static_cast<double>(conn.getCount(3));
      }

      static void computeAll(LocalMesh& mesh, Connectivity<Context::Local>::Mode mode =
            Connectivity<Context::Local>::Mode::Discover)
      {
        const size_t D = mesh.getDimension();
        auto& conn = mesh.getConnectivity();
        for (size_t d = 0; d <= D; ++d)
          for (size_t dp = 0; dp <= D; ++dp)
            conn.compute(d, dp, mode);
      }

      static void clearAllDerived(LocalMesh& mesh)
      {
        const size_t D = mesh.getDimension();
        auto& conn = mesh.getConnectivity();

        for (size_t d = 0; d <= D; ++d)
        {
          for (size_t dp = 0; dp <= D; ++dp)
          {
            if ((d == D && dp == 0) || (d == D && dp == D))
              continue;
            conn.clear(d, dp);
          }
        }
      }

      static void prepareForTranspose(LocalMesh& mesh, size_t d, size_t dp)
      {
        assert(d < dp);
        auto& conn = mesh.getConnectivity();

        // Need dp -> d first, then benchmark transpose(d, dp)
        conn.compute(dp, d);
        conn.clear(d, dp);
      }

      static void prepareForIntersection(LocalMesh& mesh, size_t d, size_t dp, size_t dpp)
      {
        auto& conn = mesh.getConnectivity();

        // intersection(d, dp, dpp) requires d -> dpp and dpp -> dp
        conn.compute(d, dpp);
        conn.compute(dpp, dp);
        conn.clear(d, dp);
      }

      static void benchmarkColdCompute(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          size_t dp,
          Connectivity<Context::Local>::Mode mode = Connectivity<Context::Local>::Mode::Discover)
      {
        for (auto _ : st)
        {
          st.PauseTiming();
          auto mesh = makeMesh(spec);
          st.ResumeTiming();

          mesh.getConnectivity().compute(d, dp, mode);
        }

        auto mesh = makeMesh(spec);
        mesh.getConnectivity().compute(d, dp, mode);
        setCommonCounters(st, mesh);
        st.counters["from"] = static_cast<double>(d);
        st.counters["to"]   = static_cast<double>(dp);
      }

      static void benchmarkWarmCompute(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          size_t dp,
          Connectivity<Context::Local>::Mode mode = Connectivity<Context::Local>::Mode::Discover)
      {
        auto mesh = makeMesh(spec);
        mesh.getConnectivity().compute(d, dp, mode);

        for (auto _ : st)
        {
          auto& result = mesh.getConnectivity().compute(d, dp, mode);
          benchmark::DoNotOptimize(result);
          benchmark::ClobberMemory();
        }

        setCommonCounters(st, mesh);
        st.counters["from"] = static_cast<double>(d);
        st.counters["to"]   = static_cast<double>(dp);
      }

      static void benchmarkColdAllPairs(
          benchmark::State& st,
          const MeshSpec& spec,
          Connectivity<Context::Local>::Mode mode = Connectivity<Context::Local>::Mode::Discover)
      {
        for (auto _ : st)
        {
          st.PauseTiming();
          auto mesh = makeMesh(spec);
          st.ResumeTiming();

          computeAll(mesh, mode);
          benchmark::DoNotOptimize(mesh.getConnectivity());
        }

        auto mesh = makeMesh(spec);
        computeAll(mesh, mode);
        setCommonCounters(st, mesh);
      }

      static void benchmarkWarmAllPairs(
          benchmark::State& st,
          const MeshSpec& spec,
          Connectivity<Context::Local>::Mode mode = Connectivity<Context::Local>::Mode::Discover)
      {
        auto mesh = makeMesh(spec);
        computeAll(mesh, mode);

        for (auto _ : st)
        {
          computeAll(mesh, mode);
          benchmark::DoNotOptimize(mesh.getConnectivity());
          benchmark::ClobberMemory();
        }

        setCommonCounters(st, mesh);
      }

      static void benchmarkColdBuild(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          Connectivity<Context::Local>::Mode mode = Connectivity<Context::Local>::Mode::Discover)
      {
        for (auto _ : st)
        {
          st.PauseTiming();
          auto mesh = makeMesh(spec);
          st.ResumeTiming();

          mesh.getConnectivity().build(d, mode);
        }

        auto mesh = makeMesh(spec);
        mesh.getConnectivity().build(d, mode);
        setCommonCounters(st, mesh);
        st.counters["build_dim"] = static_cast<double>(d);
      }

      static void benchmarkWarmBuild(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          Connectivity<Context::Local>::Mode mode = Connectivity<Context::Local>::Mode::Discover)
      {
        auto mesh = makeMesh(spec);
        mesh.getConnectivity().build(d, mode);

        for (auto _ : st)
        {
          st.PauseTiming();
          mesh.getConnectivity().clear(mesh.getDimension(), d);
          st.ResumeTiming();

          mesh.getConnectivity().build(d, mode);
        }

        setCommonCounters(st, mesh);
        st.counters["build_dim"] = static_cast<double>(d);
      }

      static void benchmarkColdTranspose(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          size_t dp)
      {
        assert(d < dp);

        for (auto _ : st)
        {
          st.PauseTiming();
          auto mesh = makeMesh(spec);
          prepareForTranspose(mesh, d, dp);
          st.ResumeTiming();

          mesh.getConnectivity().transpose(d, dp);
        }

        auto mesh = makeMesh(spec);
        prepareForTranspose(mesh, d, dp);
        mesh.getConnectivity().transpose(d, dp);
        setCommonCounters(st, mesh);
        st.counters["from"] = static_cast<double>(d);
        st.counters["to"]   = static_cast<double>(dp);
      }

      static void benchmarkWarmTranspose(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          size_t dp)
      {
        assert(d < dp);

        auto mesh = makeMesh(spec);
        prepareForTranspose(mesh, d, dp);

        for (auto _ : st)
        {
          st.PauseTiming();
          mesh.getConnectivity().clear(d, dp);
          st.ResumeTiming();

          mesh.getConnectivity().transpose(d, dp);
        }

        setCommonCounters(st, mesh);
        st.counters["from"] = static_cast<double>(d);
        st.counters["to"]   = static_cast<double>(dp);
      }

      static void benchmarkColdIntersection(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          size_t dp,
          size_t dpp)
      {
        for (auto _ : st)
        {
          st.PauseTiming();
          auto mesh = makeMesh(spec);
          prepareForIntersection(mesh, d, dp, dpp);
          st.ResumeTiming();

          mesh.getConnectivity().intersection(d, dp, dpp);
        }

        auto mesh = makeMesh(spec);
        prepareForIntersection(mesh, d, dp, dpp);
        mesh.getConnectivity().intersection(d, dp, dpp);
        setCommonCounters(st, mesh);
        st.counters["from"] = static_cast<double>(d);
        st.counters["to"]   = static_cast<double>(dp);
        st.counters["via"]  = static_cast<double>(dpp);
      }

      static void benchmarkWarmIntersection(
          benchmark::State& st,
          const MeshSpec& spec,
          size_t d,
          size_t dp,
          size_t dpp)
      {
        auto mesh = makeMesh(spec);
        prepareForIntersection(mesh, d, dp, dpp);

        for (auto _ : st)
        {
          st.PauseTiming();
          mesh.getConnectivity().clear(d, dp);
          st.ResumeTiming();

          mesh.getConnectivity().intersection(d, dp, dpp);
        }

        setCommonCounters(st, mesh);
        st.counters["from"] = static_cast<double>(d);
        st.counters["to"]   = static_cast<double>(dp);
        st.counters["via"]  = static_cast<double>(dpp);
      }
  };

  constexpr ConnectivityBenchmark::MeshSpec ConnectivityBenchmark::kMeshes[];

  // --------------------------------------------------------------------------
  // 1D: EDGE / SEGMENT
  // --------------------------------------------------------------------------

  BENCHMARK_F(ConnectivityBenchmark, Edge_Cold_AllPairs)(benchmark::State& st)
  {
    benchmarkColdAllPairs(st, kMeshes[0]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Warm_AllPairs)(benchmark::State& st)
  {
    benchmarkWarmAllPairs(st, kMeshes[0]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Cold_Compute_0_0)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[0], 0, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Cold_Compute_0_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[0], 0, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Cold_Compute_1_0)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[0], 1, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Cold_Compute_1_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[0], 1, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Warm_Compute_0_0)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[0], 0, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Warm_Compute_0_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[0], 0, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Warm_Compute_1_0)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[0], 1, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Edge_Warm_Compute_1_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[0], 1, 1);
  }

  // --------------------------------------------------------------------------
  // 2D: TRIANGLE
  // --------------------------------------------------------------------------

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_AllPairs)(benchmark::State& st)
  {
    benchmarkColdAllPairs(st, kMeshes[1]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_AllPairs)(benchmark::State& st)
  {
    benchmarkWarmAllPairs(st, kMeshes[1]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Compute_2_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[1], 2, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Compute_1_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[1], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Compute_2_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[1], 2, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Compute_1_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[1], 1, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Compute_2_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[1], 2, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Compute_1_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[1], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Compute_2_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[1], 2, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Compute_1_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[1], 1, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Build_1)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[1], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Build_1)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[1], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Transpose_1_2)(benchmark::State& st)
  {
    benchmarkColdTranspose(st, kMeshes[1], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Transpose_1_2)(benchmark::State& st)
  {
    benchmarkWarmTranspose(st, kMeshes[1], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Intersection_2_2_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[1], 2, 2, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Intersection_2_2_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[1], 2, 2, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Cold_Intersection_1_1_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[1], 1, 1, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Triangle_Warm_Intersection_1_1_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[1], 1, 1, 0);
  }

  // --------------------------------------------------------------------------
  // 2D: QUADRILATERAL
  // --------------------------------------------------------------------------

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_AllPairs)(benchmark::State& st)
  {
    benchmarkColdAllPairs(st, kMeshes[2]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_AllPairs)(benchmark::State& st)
  {
    benchmarkWarmAllPairs(st, kMeshes[2]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Compute_2_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[2], 2, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Compute_1_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[2], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Compute_2_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[2], 2, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Compute_1_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[2], 1, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Compute_2_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[2], 2, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Compute_1_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[2], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Compute_2_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[2], 2, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Compute_1_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[2], 1, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Build_1)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[2], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Build_1)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[2], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Transpose_1_2)(benchmark::State& st)
  {
    benchmarkColdTranspose(st, kMeshes[2], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Transpose_1_2)(benchmark::State& st)
  {
    benchmarkWarmTranspose(st, kMeshes[2], 1, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Intersection_2_2_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[2], 2, 2, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Intersection_2_2_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[2], 2, 2, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Cold_Intersection_1_1_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[2], 1, 1, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Quadrilateral_Warm_Intersection_1_1_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[2], 1, 1, 0);
  }

  // --------------------------------------------------------------------------
  // 3D: TETRAHEDRON
  // --------------------------------------------------------------------------

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_AllPairs)(benchmark::State& st)
  {
    benchmarkColdAllPairs(st, kMeshes[3]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_AllPairs)(benchmark::State& st)
  {
    benchmarkWarmAllPairs(st, kMeshes[3]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Compute_3_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[3], 3, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Compute_3_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[3], 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Compute_2_3)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[3], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Compute_3_3)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[3], 3, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Compute_3_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[3], 3, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Compute_3_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[3], 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Compute_2_3)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[3], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Compute_3_3)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[3], 3, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Build_1)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[3], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Build_1)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[3], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Build_2)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[3], 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Build_2)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[3], 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Transpose_2_3)(benchmark::State& st)
  {
    benchmarkColdTranspose(st, kMeshes[3], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Transpose_2_3)(benchmark::State& st)
  {
    benchmarkWarmTranspose(st, kMeshes[3], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Intersection_3_3_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[3], 3, 3, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Intersection_3_3_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[3], 3, 3, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Cold_Intersection_3_3_via_2)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[3], 3, 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Tetrahedron_Warm_Intersection_3_3_via_2)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[3], 3, 3, 2);
  }

  // --------------------------------------------------------------------------
  // 3D: HEXAHEDRON
  // --------------------------------------------------------------------------

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_AllPairs)(benchmark::State& st)
  {
    benchmarkColdAllPairs(st, kMeshes[4]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_AllPairs)(benchmark::State& st)
  {
    benchmarkWarmAllPairs(st, kMeshes[4]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Compute_3_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[4], 3, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Compute_3_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[4], 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Compute_2_3)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[4], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Compute_3_3)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[4], 3, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Compute_3_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[4], 3, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Compute_3_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[4], 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Compute_2_3)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[4], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Compute_3_3)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[4], 3, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Build_1)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[4], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Build_1)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[4], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Build_2)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[4], 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Build_2)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[4], 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Transpose_2_3)(benchmark::State& st)
  {
    benchmarkColdTranspose(st, kMeshes[4], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Transpose_2_3)(benchmark::State& st)
  {
    benchmarkWarmTranspose(st, kMeshes[4], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Intersection_3_3_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[4], 3, 3, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Intersection_3_3_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[4], 3, 3, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Cold_Intersection_3_3_via_2)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[4], 3, 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Hexahedron_Warm_Intersection_3_3_via_2)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[4], 3, 3, 2);
  }

  // --------------------------------------------------------------------------
  // 3D: WEDGE
  // --------------------------------------------------------------------------

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_AllPairs)(benchmark::State& st)
  {
    benchmarkColdAllPairs(st, kMeshes[5]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_AllPairs)(benchmark::State& st)
  {
    benchmarkWarmAllPairs(st, kMeshes[5]);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Compute_3_1)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[5], 3, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Compute_3_2)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[5], 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Compute_2_3)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[5], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Compute_3_3)(benchmark::State& st)
  {
    benchmarkColdCompute(st, kMeshes[5], 3, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Compute_3_1)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[5], 3, 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Compute_3_2)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[5], 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Compute_2_3)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[5], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Compute_3_3)(benchmark::State& st)
  {
    benchmarkWarmCompute(st, kMeshes[5], 3, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Build_1)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[5], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Build_1)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[5], 1);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Build_2)(benchmark::State& st)
  {
    benchmarkColdBuild(st, kMeshes[5], 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Build_2)(benchmark::State& st)
  {
    benchmarkWarmBuild(st, kMeshes[5], 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Transpose_2_3)(benchmark::State& st)
  {
    benchmarkColdTranspose(st, kMeshes[5], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Transpose_2_3)(benchmark::State& st)
  {
    benchmarkWarmTranspose(st, kMeshes[5], 2, 3);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Intersection_3_3_via_0)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[5], 3, 3, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Intersection_3_3_via_0)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[5], 3, 3, 0);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Cold_Intersection_3_3_via_2)(benchmark::State& st)
  {
    benchmarkColdIntersection(st, kMeshes[5], 3, 3, 2);
  }

  BENCHMARK_F(ConnectivityBenchmark, Wedge_Warm_Intersection_3_3_via_2)(benchmark::State& st)
  {
    benchmarkWarmIntersection(st, kMeshes[5], 3, 3, 2);
  }
}

BENCHMARK_MAIN();
