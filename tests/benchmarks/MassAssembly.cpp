/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Full Eigen mass/reaction operator and projection-load assembly. */
#include "Rodin/Assembly.h"
#include "MassWorkload.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool massCheckFailed = false;

  /**
   * @brief Times one public operator or load assembly, with independent oracles.
   * Architecture: setup, original-loop assembly, constant-field energy and
   * projection identities are untimed. Repeated assembly must replace output.
   * A doubled reaction coefficient is checked after timing to detect stale
   * expression values. No projection solve or constraint elimination is timed.
   */
  template <size_t K, class Scalar, bool Vector, MassFamily Family>
  class MassAssembly
  {
    public:
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        bool reaction, bool rhs)
      {
        using namespace Variational;
        using Workload = MassWorkload<K, Scalar, Vector, Family>;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        auto space = Workload::makeSpace(mesh, dim);
        TrialFunction u(space);
        TestFunction v(space);
        BilinearForm matrix(u, v);
        LinearForm load(v);
        Real scale = 1;
        const Scalar expected =
          Workload::configure(u, v, matrix, load, dim, reaction, scale);
        matrix.assemble();
        load.assemble();
        BilinearForm reference(u, v);
        LinearForm referenceLoad(v);
        Workload::template configure<true>(
          u, v, reference, referenceLoad, dim, reaction, scale);
        reference.assemble();
        const auto baseline = reference.getOperator();
        const auto baselineLoad = load.getVector();
        const auto& c = u.getSolution().getData();
        const auto valid = [&]() {
          const auto& a = matrix.getOperator();
          const auto& b = load.getVector();
          const Scalar energy = c.dot(a * c);
          return std::isfinite(std::abs(energy)) &&
            std::abs(energy / (scale * expected) - Scalar(1)) < 1e-10 &&
            (a - scale * baseline).norm() <
            1e-11 * std::max(Real(1), baseline.norm() * scale) &&
            (a * c - b).norm() < 1e-11 * std::max(Real(1), b.norm()) &&
            (b - scale * baselineLoad).norm() < 1e-11 * std::max(Real(1), b.norm());
        };
        matrix.assemble();
        load.assemble();
        if (!valid())
        {
          massCheckFailed = true;
          state.SkipWithError("Mass/reaction operator, load, or analytic oracle failed");
          return;
        }
        for (auto _ : state)
        {
          if (rhs)
          {
            load.assemble();
            benchmark::DoNotOptimize(load.getVector().data());
          }
          else
          {
            matrix.assemble();
            benchmark::DoNotOptimize(matrix.getOperator().nonZeros());
          }
          benchmark::ClobberMemory();
        }
        if (!valid())
        {
          massCheckFailed = true;
          state.SkipWithError("Timed assembly changed operator or load");
        }
        if (reaction)
        {
          scale = 2;
          matrix.assemble();
          load.assemble();
          if (!valid())
          {
            massCheckFailed = true;
            state.SkipWithError("Coefficient update did not replace operator and load");
          }
        }
        state.counters["cells"] = mesh.getPolytopeCount(dim);
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = matrix.getOperator().nonZeros();
        state.counters["degree"] = K;
        state.counters["components"] = Vector ? dim : 1;
        state.counters["quadrature_order"] = 8;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(8, geometry).getSize();
      }

      static void registerCases()
      {
        MassWorkload<K, Scalar, Vector, Family>::registerCases(run, false);
      }
  };

  template <class Scalar, bool Vector>
  void registerMass()
  {
    MassAssembly<1, Scalar, Vector, MassFamily::H1>::registerCases();
    MassAssembly<2, Scalar, Vector, MassFamily::H1>::registerCases();
    MassAssembly<3, Scalar, Vector, MassFamily::H1>::registerCases();
    MassAssembly<1, Scalar, Vector, MassFamily::P1>::registerCases();
    MassAssembly<0, Scalar, Vector, MassFamily::P0>::registerCases();
    MassAssembly<0, Scalar, Vector, MassFamily::P0g>::registerCases();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerMass<Rodin::Real, false>();
  registerMass<Rodin::Real, true>();
  registerMass<Rodin::Complex, false>();
  registerMass<Rodin::Complex, true>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential");
#endif
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Warmed operator or projection-load assembly; no solve");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !cases || massCheckFailed ? 1 : 0;
}
