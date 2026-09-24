/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <benchmark/benchmark.h>

#include "Rodin/Assembly.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    LocalMesh unitSquare(size_t n)
    {
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    LocalMesh unitCube(size_t n)
    {
      LocalMesh mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    template <class FES>
    void integralMass(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = Integral(u, v);
      form.assemble();

      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void namedMass(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      MassForm form(u, v);

      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void integralDiffusion(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      form = Integral(Grad(u), Grad(v));
      form.assemble();

      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void namedDiffusion(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);
      DiffusionForm form(u, v);

      for (auto _ : state)
      {
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void integralMassFirst(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);

      for (auto _ : state)
      {
        BilinearForm form(u, v);
        form = Integral(u, v);
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void namedMassFirst(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);

      for (auto _ : state)
      {
        MassForm form(u, v);
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void integralDiffusionFirst(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);

      for (auto _ : state)
      {
        BilinearForm form(u, v);
        form = Integral(Grad(u), Grad(v));
        form.assemble();
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }

    template <class FES>
    void namedDiffusionFirst(benchmark::State& state, FES& fes)
    {
      TrialFunction u(fes);
      TestFunction v(fes);

      for (auto _ : state)
      {
        DiffusionForm form(u, v);
        benchmark::DoNotOptimize(form.getOperator().nonZeros());
      }
    }
  }

  /// @brief Measures 2D reassembly of the P1 mass form built through Integral.
  static void P1IntegralMass(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    integralMass(state, fes);
  }
  /// @brief Registers the P1IntegralMass benchmark.
  BENCHMARK(P1IntegralMass);

  /// @brief Measures 2D reassembly of the P1 mass named form.
  static void P1NamedMass(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    namedMass(state, fes);
  }
  /// @brief Registers the P1NamedMass benchmark.
  BENCHMARK(P1NamedMass);

  /// @brief Measures 2D reassembly of the P1 diffusion form built through Integral.
  static void P1IntegralDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    integralDiffusion(state, fes);
  }
  /// @brief Registers the P1IntegralDiffusion benchmark.
  BENCHMARK(P1IntegralDiffusion);

  /// @brief Measures 2D reassembly of the P1 diffusion named form.
  static void P1NamedDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    namedDiffusion(state, fes);
  }
  /// @brief Registers the P1NamedDiffusion benchmark.
  BENCHMARK(P1NamedDiffusion);

  /// @brief Measures 2D reassembly of the H1 order 2 mass form built through Integral.
  static void H1P2IntegralMass(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralMass(state, fes);
  }
  /// @brief Registers the H1P2IntegralMass benchmark.
  BENCHMARK(H1P2IntegralMass);

  /// @brief Measures 2D reassembly of the H1 order 2 mass named form.
  static void H1P2NamedMass(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedMass(state, fes);
  }
  /// @brief Registers the H1P2NamedMass benchmark.
  BENCHMARK(H1P2NamedMass);

  /// @brief Measures 2D reassembly of the H1 order 2 diffusion form built through Integral.
  static void H1P2IntegralDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralDiffusion(state, fes);
  }
  /// @brief Registers the H1P2IntegralDiffusion benchmark.
  BENCHMARK(H1P2IntegralDiffusion);

  /// @brief Measures 2D reassembly of the H1 order 2 diffusion named form.
  static void H1P2NamedDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedDiffusion(state, fes);
  }
  /// @brief Registers the H1P2NamedDiffusion benchmark.
  BENCHMARK(H1P2NamedDiffusion);

  /// @brief Measures 3D reassembly of the P1 mass form built through Integral.
  static void P1IntegralMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    integralMass(state, fes);
  }
  /// @brief Registers the P1IntegralMass3D benchmark.
  BENCHMARK(P1IntegralMass3D);

  /// @brief Measures 3D reassembly of the P1 mass named form.
  static void P1NamedMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    namedMass(state, fes);
  }
  /// @brief Registers the P1NamedMass3D benchmark.
  BENCHMARK(P1NamedMass3D);

  /// @brief Measures 3D reassembly of the P1 diffusion form built through Integral.
  static void P1IntegralDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    integralDiffusion(state, fes);
  }
  /// @brief Registers the P1IntegralDiffusion3D benchmark.
  BENCHMARK(P1IntegralDiffusion3D);

  /// @brief Measures 3D reassembly of the P1 diffusion named form.
  static void P1NamedDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    namedDiffusion(state, fes);
  }
  /// @brief Registers the P1NamedDiffusion3D benchmark.
  BENCHMARK(P1NamedDiffusion3D);

  /// @brief Measures 3D reassembly of the H1 order 2 mass form built through Integral.
  static void H1P2IntegralMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralMass(state, fes);
  }
  /// @brief Registers the H1P2IntegralMass3D benchmark.
  BENCHMARK(H1P2IntegralMass3D);

  /// @brief Measures 3D reassembly of the H1 order 2 mass named form.
  static void H1P2NamedMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedMass(state, fes);
  }
  /// @brief Registers the H1P2NamedMass3D benchmark.
  BENCHMARK(H1P2NamedMass3D);

  /// @brief Measures 3D reassembly of the H1 order 2 diffusion form built through Integral.
  static void H1P2IntegralDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralDiffusion(state, fes);
  }
  /// @brief Registers the H1P2IntegralDiffusion3D benchmark.
  BENCHMARK(H1P2IntegralDiffusion3D);

  /// @brief Measures 3D reassembly of the H1 order 2 diffusion named form.
  static void H1P2NamedDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedDiffusion(state, fes);
  }
  /// @brief Registers the H1P2NamedDiffusion3D benchmark.
  BENCHMARK(H1P2NamedDiffusion3D);
  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 mass form built through Integral.
  static void P1IntegralMassFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    integralMassFirst(state, fes);
  }
  /// @brief Registers the P1IntegralMassFirst benchmark.
  BENCHMARK(P1IntegralMassFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 mass named form.
  static void P1NamedMassFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    namedMassFirst(state, fes);
  }
  /// @brief Registers the P1NamedMassFirst benchmark.
  BENCHMARK(P1NamedMassFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 diffusion form built through Integral.
  static void P1IntegralDiffusionFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    integralDiffusionFirst(state, fes);
  }
  /// @brief Registers the P1IntegralDiffusionFirst benchmark.
  BENCHMARK(P1IntegralDiffusionFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 diffusion named form.
  static void P1NamedDiffusionFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    namedDiffusionFirst(state, fes);
  }
  /// @brief Registers the P1NamedDiffusionFirst benchmark.
  BENCHMARK(P1NamedDiffusionFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 mass form built through Integral.
  static void H1P2IntegralMassFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralMassFirst(state, fes);
  }
  /// @brief Registers the H1P2IntegralMassFirst benchmark.
  BENCHMARK(H1P2IntegralMassFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 mass named form.
  static void H1P2NamedMassFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedMassFirst(state, fes);
  }
  /// @brief Registers the H1P2NamedMassFirst benchmark.
  BENCHMARK(H1P2NamedMassFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 diffusion form built through Integral.
  static void H1P2IntegralDiffusionFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralDiffusionFirst(state, fes);
  }
  /// @brief Registers the H1P2IntegralDiffusionFirst benchmark.
  BENCHMARK(H1P2IntegralDiffusionFirst);

  /// @brief Measures 2D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 diffusion named form.
  static void H1P2NamedDiffusionFirst(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedDiffusionFirst(state, fes);
  }
  /// @brief Registers the H1P2NamedDiffusionFirst benchmark.
  BENCHMARK(H1P2NamedDiffusionFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 mass form built through Integral.
  static void P1IntegralMass3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    integralMassFirst(state, fes);
  }
  /// @brief Registers the P1IntegralMass3DFirst benchmark.
  BENCHMARK(P1IntegralMass3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 mass named form.
  static void P1NamedMass3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    namedMassFirst(state, fes);
  }
  /// @brief Registers the P1NamedMass3DFirst benchmark.
  BENCHMARK(P1NamedMass3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 diffusion form built through Integral.
  static void P1IntegralDiffusion3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    integralDiffusionFirst(state, fes);
  }
  /// @brief Registers the P1IntegralDiffusion3DFirst benchmark.
  BENCHMARK(P1IntegralDiffusion3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the P1 diffusion named form.
  static void P1NamedDiffusion3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    namedDiffusionFirst(state, fes);
  }
  /// @brief Registers the P1NamedDiffusion3DFirst benchmark.
  BENCHMARK(P1NamedDiffusion3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 mass form built through Integral.
  static void H1P2IntegralMass3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralMassFirst(state, fes);
  }
  /// @brief Registers the H1P2IntegralMass3DFirst benchmark.
  BENCHMARK(H1P2IntegralMass3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 mass named form.
  static void H1P2NamedMass3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedMassFirst(state, fes);
  }
  /// @brief Registers the H1P2NamedMass3DFirst benchmark.
  BENCHMARK(H1P2NamedMass3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 diffusion form built through Integral.
  static void H1P2IntegralDiffusion3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralDiffusionFirst(state, fes);
  }
  /// @brief Registers the H1P2IntegralDiffusion3DFirst benchmark.
  BENCHMARK(H1P2IntegralDiffusion3DFirst);

  /// @brief Measures 3D construction and first assembly, on a mesh whose quadrature cache is already warm, of the H1 order 2 diffusion named form.
  static void H1P2NamedDiffusion3DFirst(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedDiffusionFirst(state, fes);
  }
  /// @brief Registers the H1P2NamedDiffusion3DFirst benchmark.
  BENCHMARK(H1P2NamedDiffusion3DFirst);
  /// @brief Measures 2D construction and first assembly of the P1 mass form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1IntegralMassCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(128);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(u, v);
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1IntegralMassCold benchmark.
  BENCHMARK(P1IntegralMassCold);

  /// @brief Measures 2D construction and first assembly of the P1 mass named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1NamedMassCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(128);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          MassForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1NamedMassCold benchmark.
  BENCHMARK(P1NamedMassCold);

  /// @brief Measures 2D construction and first assembly of the P1 diffusion form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1IntegralDiffusionCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(128);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(Grad(u), Grad(v));
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1IntegralDiffusionCold benchmark.
  BENCHMARK(P1IntegralDiffusionCold);

  /// @brief Measures 2D construction and first assembly of the P1 diffusion named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1NamedDiffusionCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(128);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          DiffusionForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1NamedDiffusionCold benchmark.
  BENCHMARK(P1NamedDiffusionCold);

  /// @brief Measures 2D construction and first assembly of the H1 order 2 mass form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2IntegralMassCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(32);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(u, v);
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2IntegralMassCold benchmark.
  BENCHMARK(H1P2IntegralMassCold);

  /// @brief Measures 2D construction and first assembly of the H1 order 2 mass named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2NamedMassCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(32);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          MassForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2NamedMassCold benchmark.
  BENCHMARK(H1P2NamedMassCold);

  /// @brief Measures 2D construction and first assembly of the H1 order 2 diffusion form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2IntegralDiffusionCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(32);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(Grad(u), Grad(v));
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2IntegralDiffusionCold benchmark.
  BENCHMARK(H1P2IntegralDiffusionCold);

  /// @brief Measures 2D construction and first assembly of the H1 order 2 diffusion named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2NamedDiffusionCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitSquare(32);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          DiffusionForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2NamedDiffusionCold benchmark.
  BENCHMARK(H1P2NamedDiffusionCold);

  /// @brief Measures 3D construction and first assembly of the P1 mass form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1IntegralMass3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(24);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(u, v);
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1IntegralMass3DCold benchmark.
  BENCHMARK(P1IntegralMass3DCold);

  /// @brief Measures 3D construction and first assembly of the P1 mass named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1NamedMass3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(24);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          MassForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1NamedMass3DCold benchmark.
  BENCHMARK(P1NamedMass3DCold);

  /// @brief Measures 3D construction and first assembly of the P1 diffusion form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1IntegralDiffusion3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(24);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(Grad(u), Grad(v));
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1IntegralDiffusion3DCold benchmark.
  BENCHMARK(P1IntegralDiffusion3DCold);

  /// @brief Measures 3D construction and first assembly of the P1 diffusion named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void P1NamedDiffusion3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(24);
        P1 fes(mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          DiffusionForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the P1NamedDiffusion3DCold benchmark.
  BENCHMARK(P1NamedDiffusion3DCold);

  /// @brief Measures 3D construction and first assembly of the H1 order 2 mass form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2IntegralMass3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(10);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(u, v);
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2IntegralMass3DCold benchmark.
  BENCHMARK(H1P2IntegralMass3DCold);

  /// @brief Measures 3D construction and first assembly of the H1 order 2 mass named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2NamedMass3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(10);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          MassForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2NamedMass3DCold benchmark.
  BENCHMARK(H1P2NamedMass3DCold);

  /// @brief Measures 3D construction and first assembly of the H1 order 2 diffusion form built through Integral
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2IntegralDiffusion3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(10);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          BilinearForm form(u, v);
          form = Integral(Grad(u), Grad(v));
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2IntegralDiffusion3DCold benchmark.
  BENCHMARK(H1P2IntegralDiffusion3DCold);

  /// @brief Measures 3D construction and first assembly of the H1 order 2 diffusion named form
  /// on a fresh mesh, so that every polytope's quadrature is mapped during the
  /// measurement; building the mesh and the space is not timed.
  static void H1P2NamedDiffusion3DCold(benchmark::State& state)
  {
    for (auto _ : state)
    {
      state.PauseTiming();
      {
        auto mesh = unitCube(10);
        H1 fes(std::integral_constant<size_t, 2>{}, mesh);
        TrialFunction u(fes);
        TestFunction v(fes);
        state.ResumeTiming();
        {
          DiffusionForm form(u, v);
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
        }
        state.PauseTiming();
      }
      state.ResumeTiming();
    }
  }
  /// @brief Registers the H1P2NamedDiffusion3DCold benchmark.
  BENCHMARK(H1P2NamedDiffusion3DCold);
}
