/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <string>
#include <type_traits>

#include "Rodin/Geometry.h"
#include "Rodin/Solid.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  /// @brief Measures complete assembly with a full-matrix oracle outside timing.
  /// @details Kind: scalar mass, diffusion, Helmholtz, vector mass, elasticity.
  /// First includes construction, assembly and destruction; mesh quadrature is warm.
  template <size_t Order, size_t Kind, bool FieldCoefficient>
  void namedEvaluation(benchmark::State& state)
  {
    const size_t dimension = state.range(0);
    auto mesh = dimension == 2 ?
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {17, 17}) :
      LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {5, 5, 5});
    for (size_t d = dimension; d > 0; --d)
      mesh.getConnectivity().compute(d, d - 1);
    using Range = std::conditional_t<(Kind >= 3), Math::SpatialVector<Real>, Real>;
    using Space = std::conditional_t<Order == 1, P1<Range>, H1<Order, Range>>;
    auto makeSpace = [&] {
      if constexpr (Order == 1 && Kind >= 3)
        return Space(mesh, dimension);
      else if constexpr (Order == 1)
        return Space(mesh);
      else if constexpr (Kind >= 3)
        return Space(std::integral_constant<size_t, Order>{}, mesh, dimension);
      else
        return Space(std::integral_constant<size_t, Order>{}, mesh);
    };
    auto fes = makeSpace();
    using CoefficientSpace = std::conditional_t<Order == 1, P1<Real>, H1<Order, Real>>;
    auto coefficientSpace = [&] {
      if constexpr (Order == 1)
        return CoefficientSpace(mesh);
      else
        return CoefficientSpace(std::integral_constant<size_t, Order>{}, mesh);
    }();
    GridFunction field(coefficientSpace);
    field.project(RealFunction([](const Point& p) {
      Real value = 1;
      for (size_t d = 0; d < p.getCoordinates().size(); ++d)
        value += (d + 1) * p.getCoordinates()(d) * p.getCoordinates()(d);
      return value;
    }));
    const auto constant = RealFunction(Real(2));
    const auto& coefficient = [&]() -> decltype(auto) {
      if constexpr (FieldCoefficient)
        return (field);
      else
        return (constant);
    }();
    TrialFunction u(fes);
    TestFunction v(fes);
    auto makeNamed = [&] {
      if constexpr (Kind == 0 || Kind == 3)
        return MassForm(coefficient, u, v);
      else if constexpr (Kind == 1)
        return DiffusionForm(coefficient, u, v);
      else if constexpr (Kind == 2)
        return HelmholtzForm(coefficient, coefficient, u, v);
      else
        return LinearElasticityForm(coefficient, coefficient, u, v);
    };
    auto makeIntegral = [&] {
      BilinearForm form(u, v);
      if constexpr (Kind == 0 || Kind == 3)
        form = Integral(coefficient * u, v);
      else if constexpr (Kind == 1)
        form = Integral(coefficient * Grad(u), Grad(v));
      else if constexpr (Kind == 2)
        form = Integral(coefficient * Grad(u), Grad(v)) + Integral(coefficient * u, v);
      else
        form = LinearElasticityIntegral(u, v)(coefficient, coefficient);
      return form;
    };
    auto named = makeNamed(); // Named constructors assemble automatically.
    auto integral = makeIntegral();
    integral.assemble();
    const Math::SparseMatrix<Real> expected = integral.getOperator();
    const auto error = [&](const auto& matrix) {
      return (matrix - expected).norm() / std::max(Real(1), expected.norm());
    };
    const Real initialError = error(named.getOperator());
    if (!std::isfinite(initialError) || initialError > 1e-11 || expected.norm() == 0)
    {
      state.SkipWithError("Named operator differs from Integral before timing");
      return;
    }
    Math::SparseMatrix<Real> last;
    const bool first = state.range(1) != 0;
    const bool useIntegral = state.range(2) != 0;
    for (auto _ : state)
    {
      if (first)
      {
        if (useIntegral)
        {
          auto form = makeIntegral();
          form.assemble();
          last.swap(form.getOperator());
        }
        else
        {
          auto form = makeNamed();
          last.swap(form.getOperator());
        }
        benchmark::DoNotOptimize(last.valuePtr());
      }
      else if (useIntegral)
      {
        integral.assemble();
        benchmark::DoNotOptimize(integral.getOperator().valuePtr());
      }
      else
      {
        named.assemble();
        benchmark::DoNotOptimize(named.getOperator().valuePtr());
      }
      benchmark::ClobberMemory();
    }
    const Real finalError = first ? error(last) :
      useIntegral ? error(integral.getOperator()) : error(named.getOperator());
    if (!std::isfinite(finalError) || finalError > 1e-11)
      state.SkipWithError("Final assembled operator differs from frozen Integral oracle");
    state.counters["relative_error"] = std::max(initialError, finalError);
    state.counters["cells"] = mesh.getCellCount();
    state.counters["dofs"] = expected.rows();
    state.counters["nonzeros"] = expected.nonZeros();
    state.counters["named_nonzeros"] = named.getOperator().nonZeros();
    state.counters["scatter_slots"] = mesh.getCellCount() *
      fes.getFiniteElement(dimension, 0).getCount() * fes.getFiniteElement(dimension, 0).getCount();
  }

  template <size_t Order, size_t Kind>
  void registerEvaluation()
  {
    const char* kinds[] = {"Mass", "Diffusion", "Helmholtz", "VectorMass", "Elasticity"};
    for (size_t coefficient = 0; coefficient < 2; ++coefficient)
      for (size_t first = 0; first < 2; ++first)
        for (size_t integral = 0; integral < 2; ++integral)
        {
          const std::string name = "NamedEvaluation/P" + std::to_string(Order) + "/" +
            kinds[Kind] + (coefficient ? "/FieldCoefficient" : "/Constant") +
            (integral ? "/Integral" : "/Named") + (first ? "/First" : "/Warm");
          auto* registration = coefficient ?
            benchmark::RegisterBenchmark(name.c_str(), &namedEvaluation<Order, Kind, true>) :
            benchmark::RegisterBenchmark(name.c_str(), &namedEvaluation<Order, Kind, false>);
          registration->Args({2, static_cast<int64_t>(first), static_cast<int64_t>(integral)})
            ->Args({3, static_cast<int64_t>(first), static_cast<int64_t>(integral)})
            ->ArgNames({"dimension", "first", "integral"})->Unit(benchmark::kMillisecond);
        }
  }

  const bool registeredEvaluation = [] {
    registerEvaluation<1, 0>(); registerEvaluation<2, 0>();
    registerEvaluation<1, 1>(); registerEvaluation<2, 1>();
    registerEvaluation<1, 2>(); registerEvaluation<2, 2>();
    registerEvaluation<1, 3>(); registerEvaluation<2, 3>();
    registerEvaluation<1, 4>(); registerEvaluation<2, 4>();
    return true;
  }();
}
