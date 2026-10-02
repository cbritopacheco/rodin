/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Paired expression-reuse/original-loop binding benchmarks. */
#include <benchmark/benchmark.h>
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P0g.h"
#include "../convergence/Convergence.h"
#include "../QuadratureReference.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  bool expressionCheckFailed = false;

  template <size_t K, class Scalar, bool Vector, size_t Family>
  void expression(benchmark::State& state, Polytope::Type geometry, bool original)
  {
    auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
    const size_t d = mesh.getDimension();
    using Range = std::conditional_t<Vector, Math::SpatialVector<Scalar>, Scalar>;
    auto space = [&]() {
      if constexpr (Family == 0)
      {
        if constexpr (Vector)
          return H1<K, Range, LocalMesh>(std::integral_constant<size_t, K>{}, mesh, d);
        else
          return H1<K, Range, LocalMesh>(std::integral_constant<size_t, K>{}, mesh);
      }
      else
      {
        using F = std::conditional_t<Family == 1, P1<Range, LocalMesh>,
          std::conditional_t<Family == 2, P0<Range, LocalMesh>, P0g<Range, LocalMesh>>>;
        if constexpr (Vector)
          return F(mesh, d);
        else
          return F(mesh);
      }
    }();
    TrialFunction u(space);
    TestFunction v(space);
    auto integral = [&]() {
      if constexpr (Vector && Family < 2)
        return Integral(
          0.5 * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
      else
        return Integral(1.25 * (u + u), v - 0.25 * v);
    }();
    integral.setOrder(6);
    QuadratureReference baseline(integral);
    const auto valid = [&]() {
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        integral.setPolytope(*cell);
        baseline.assemble(*cell, 6);
        const auto& expected = baseline.getOperator();
        for (Eigen::Index te = 0; te < expected.rows(); ++te)
          for (Eigen::Index tr = 0; tr < expected.cols(); ++tr)
            if (integral.integrate(tr, te) != expected(te, tr))
              return false;
      }
      return true;
    };
    if (!valid())
    {
      expressionCheckFailed = true;
      state.SkipWithError("Expression kernel differs from original local entries");
      return;
    }
    for (auto _ : state)
    {
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        if (original)
        {
          baseline.assemble(*cell, 6);
          benchmark::DoNotOptimize(baseline.getOperator().data());
        }
        else
        {
          integral.setPolytope(*cell);
          benchmark::DoNotOptimize(integral.integrate(0, 0));
        }
      }
      benchmark::ClobberMemory();
    }
    if (!valid())
    {
      expressionCheckFailed = true;
      state.SkipWithError("Timed expression kernel differs from original local entries");
    }
    state.counters["cells"] = mesh.getPolytopeCount(d);
    state.counters["dofs"] = space.getSize();
    state.counters["degree"] = K;
    state.counters["components"] = Vector ? d : 1;
    state.counters["quadrature_order"] = 6;
    state.counters["quadrature_points_per_cell"] =
      QF::PolytopeQuadratureFormula::get(6, geometry).getSize();
  }

  template <size_t K, class Scalar, bool Vector, size_t Family>
  void registerExpressions()
  {
    for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      for (bool original : {false, true})
      {
        const std::string family = Family == 0 ? "H1/P" + std::to_string(K)
          : Family == 1                        ? "P1"
          : Family == 2                        ? "P0"
                                               : "P0g";
        const std::string name = std::string(original ? "Original/" : "Reuse/") + family +
          "/" + (std::is_same_v<Scalar, Complex> ? "Complex" : "Real") +
          (Vector ? "Vector/" : "Scalar/") +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        benchmark::RegisterBenchmark(name.c_str(),
          [geometry, original](auto& state) {
            expression<K, Scalar, Vector, Family>(state, geometry, original);
          })
          ->Arg(2)
          ->Arg(3)
          ->Arg(5)
          ->UseRealTime();
      }
  }

  template <class Scalar, bool Vector>
  void registerRanges()
  {
    registerExpressions<1, Scalar, Vector, 0>();
    registerExpressions<2, Scalar, Vector, 0>();
    registerExpressions<3, Scalar, Vector, 0>();
    registerExpressions<1, Scalar, Vector, 1>();
    registerExpressions<0, Scalar, Vector, 2>();
    registerExpressions<0, Scalar, Vector, 3>();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerRanges<Real, false>();
  registerRanges<Real, true>();
  registerRanges<Complex, false>();
  registerRanges<Complex, true>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("scope", "Sequential local binding; no scatter or solve");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const auto matched = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !matched || expressionCheckFailed ? 1 : 0;
}
