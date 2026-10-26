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
  bool namedEvaluationFailed = false;

  /// @brief Existing relative Frobenius budget for the assembled operator oracle.
  constexpr Real OperatorTolerance = 1e-11;

  enum class GeometryScope
  {
    Retained,
    PauseOnly,
    Flushed
  };

  /**
   * @brief Measures complete assembly against a frozen full-matrix oracle.
   *
   * Contexts are scalar mass, diffusion, Helmholtz, vector mass and linear
   * elasticity, with constant or represented field coefficients. First includes
   * form construction, assembly and destruction; Warm retains the form.
   * Ordinary scopes retain mapped geometry. ColdGeometry flushes the mesh's
   * transformation and mapped-quadrature caches before each assembly, outside
   * timing. It includes lazy geometry reconstruction, not merely index-slot
   * allocation, and does not alter topology, spaces or coefficients. These
   * controls must not be interpreted as pure cache allocation latencies.
   * PauseOnly uses identical timer boundaries without changing mesh state,
   * providing a negative control for the measurement protocol.
   * Every final sparse operator is compared with the same precomputed Integral
   * oracle. No solve or error integration is timed. Profilers are opt-in and
   * are not CTest or CI timing gates.
   *
   * @tparam Order Field degree, using native P1 or H1 degree two.
   * @tparam Kind Physical form context.
   * @tparam FieldCoefficient Whether the coefficient is a grid function.
   * @param state Dimension, construction and form-family selection.
   * @param geometryScope Retained geometry, timer control or flushed geometry.
   */
  template <size_t Order, size_t Kind, bool FieldCoefficient>
  void namedEvaluation(benchmark::State& state, GeometryScope geometryScope)
  {
    const size_t dimension = state.range(0);
    auto mesh = dimension == 2
      ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {17, 17})
      : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {5, 5, 5});
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
    if (!std::isfinite(initialError) || initialError > OperatorTolerance ||
      expected.norm() == 0)
    {
      namedEvaluationFailed = true;
      state.SkipWithError("Named operator differs from Integral before timing");
      return;
    }
    Math::SparseMatrix<Real> last;
    const bool first = state.range(1) != 0;
    const bool useIntegral = state.range(2) != 0;
    for (auto _ : state)
    {
      if (geometryScope != GeometryScope::Retained)
      {
        state.PauseTiming();
        if (geometryScope == GeometryScope::Flushed)
          mesh.flush();
        state.ResumeTiming();
      }
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
    const Real finalError = first ? error(last)
      : useIntegral               ? error(integral.getOperator())
                                  : error(named.getOperator());
    if (!std::isfinite(finalError) || finalError > OperatorTolerance)
    {
      namedEvaluationFailed = true;
      state.SkipWithError("Final assembled operator differs from frozen Integral oracle");
    }
    state.counters["relative_error"] = std::max(initialError, finalError);
    state.counters["cells"] = mesh.getCellCount();
    state.counters["dofs"] = expected.rows();
    state.counters["nonzeros"] = expected.nonZeros();
    state.counters["named_nonzeros"] = named.getOperator().nonZeros();
    const auto cell = *mesh.getCell();
    const size_t coefficientDegree = coefficient.getOrder(cell).value_or(0);
    const size_t massOrder = 2 * Order + coefficientDegree;
    const size_t diffusionOrder = 2 * Order - 2 + coefficientDegree;
    const size_t massPoints = Kind == 0 || Kind == 2 || Kind == 3
      ? QF::QuadratureFormula::get(massOrder, cell.getGeometry()).getSize()
      : 0;
    // LinearElasticityForm deliberately uses basis orders, not reduced
    // derivative orders. Report the implemented rule without changing it.
    const size_t gradientOrder = Kind == 4 ? massOrder : diffusionOrder;
    const size_t gradientPoints = Kind == 1 || Kind == 2 || Kind == 4
      ? QF::QuadratureFormula::get(gradientOrder, cell.getGeometry()).getSize()
      : 0;
    state.counters["dimension"] = dimension;
    state.counters["field_degree"] = Order;
    state.counters["components"] = fes.getVectorDimension();
    state.counters["local_dofs"] = fes.getFiniteElement(dimension, 0).getCount();
    state.counters["coefficient_degree"] = coefficientDegree;
    state.counters["mass_points_per_cell"] = massPoints;
    state.counters["gradient_points_per_cell"] = gradientPoints;
    state.counters["cold_geometry"] = geometryScope == GeometryScope::Flushed;
    state.counters["paused_timing"] = geometryScope != GeometryScope::Retained;
    state.counters["scatter_slots"] = mesh.getCellCount() *
      fes.getFiniteElement(dimension, 0).getCount() *
      fes.getFiniteElement(dimension, 0).getCount();
  }

  template <size_t Order, size_t Kind>
  void registerEvaluation()
  {
    const char* kinds[] = {
      "Mass", "Diffusion", "Helmholtz", "VectorMass", "LinearElasticity"};
    for (size_t coefficient = 0; coefficient < 2; ++coefficient)
    {
      for (size_t first = 0; first < 2; ++first)
      {
        for (size_t integral = 0; integral < 2; ++integral)
        {
          for (const auto scope :
            {GeometryScope::Retained, GeometryScope::PauseOnly, GeometryScope::Flushed})
          {
            const std::string prefix = scope == GeometryScope::Retained
              ? "NamedEvaluation/P"
              : scope == GeometryScope::PauseOnly ? "NamedEvaluationTimingControl/P"
                                                  : "NamedEvaluationColdGeometry/P";
            const std::string name = prefix + std::to_string(Order) + "/" + kinds[Kind] +
              (coefficient ? "/FieldCoefficient" : "/Constant") +
              (integral ? "/Integral" : "/Named") + (first ? "/First" : "/Warm");
            auto* registration = benchmark::RegisterBenchmark(
              name.c_str(), [scope, coefficient](benchmark::State& state) {
                if (coefficient)
                  namedEvaluation<Order, Kind, true>(state, scope);
                else
                  namedEvaluation<Order, Kind, false>(state, scope);
              });
            registration
              ->Args({2, static_cast<int64_t>(first), static_cast<int64_t>(integral)})
              ->Args({3, static_cast<int64_t>(first), static_cast<int64_t>(integral)})
              ->ArgNames({"dimension", "first", "integral"})
              ->Unit(benchmark::kMillisecond);
          }
        }
      }
    }
  }

  const bool registeredEvaluation = [] {
    registerEvaluation<1, 0>();
    registerEvaluation<2, 0>();
    registerEvaluation<1, 1>();
    registerEvaluation<2, 1>();
    registerEvaluation<1, 2>();
    registerEvaluation<2, 2>();
    registerEvaluation<1, 3>();
    registerEvaluation<2, 3>();
    registerEvaluation<1, 4>();
    registerEvaluation<2, 4>();
    return true;
  }();
}

int main(int argc, char** argv)
{
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Rodin::Tests::Benchmarks::namedEvaluationFailed ? 1 : 0;
}
