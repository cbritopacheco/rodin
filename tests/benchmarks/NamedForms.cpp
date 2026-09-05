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
      LocalMesh mesh =
        LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    LocalMesh unitCube(size_t n)
    {
      LocalMesh mesh =
        LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
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
  }

  static void P1IntegralMass(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    integralMass(state, fes);
  }
  BENCHMARK(P1IntegralMass);

  static void P1NamedMass(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    namedMass(state, fes);
  }
  BENCHMARK(P1NamedMass);

  static void P1IntegralDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    integralDiffusion(state, fes);
  }
  BENCHMARK(P1IntegralDiffusion);

  static void P1NamedDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(128);
    P1 fes(mesh);
    namedDiffusion(state, fes);
  }
  BENCHMARK(P1NamedDiffusion);

  static void H1P2IntegralMass(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralMass(state, fes);
  }
  BENCHMARK(H1P2IntegralMass);

  static void H1P2NamedMass(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedMass(state, fes);
  }
  BENCHMARK(H1P2NamedMass);

  static void H1P2IntegralDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralDiffusion(state, fes);
  }
  BENCHMARK(H1P2IntegralDiffusion);

  static void H1P2NamedDiffusion(benchmark::State& state)
  {
    auto mesh = unitSquare(32);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedDiffusion(state, fes);
  }
  BENCHMARK(H1P2NamedDiffusion);

  static void P1IntegralMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    integralMass(state, fes);
  }
  BENCHMARK(P1IntegralMass3D);

  static void P1NamedMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    namedMass(state, fes);
  }
  BENCHMARK(P1NamedMass3D);

  static void P1IntegralDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    integralDiffusion(state, fes);
  }
  BENCHMARK(P1IntegralDiffusion3D);

  static void P1NamedDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(24);
    P1 fes(mesh);
    namedDiffusion(state, fes);
  }
  BENCHMARK(P1NamedDiffusion3D);

  static void H1P2IntegralMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralMass(state, fes);
  }
  BENCHMARK(H1P2IntegralMass3D);

  static void H1P2NamedMass3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedMass(state, fes);
  }
  BENCHMARK(H1P2NamedMass3D);

  static void H1P2IntegralDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    integralDiffusion(state, fes);
  }
  BENCHMARK(H1P2IntegralDiffusion3D);

  static void H1P2NamedDiffusion3D(benchmark::State& state)
  {
    auto mesh = unitCube(10);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    namedDiffusion(state, fes);
  }
  BENCHMARK(H1P2NamedDiffusion3D);
}
