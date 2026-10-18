/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Profiles P1 mass and transport binding with warm and refreshed scalar data.
 */
#include <benchmark/benchmark.h>
#include <array>
#include <cmath>
#include <string>
#include <type_traits>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "QuadratureReference.h"
#include "convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  bool p1BindingFailed = false;

  /**
   * @brief Measures the local binding stages consumed by P1 form assembly.
   *
   * ## Mathematical foundation
   *
   * On one mapped cell, the forms are @f$(u,v)_K@f$, @f$(2u,v)_K@f$ and
   * @f$(J(u)b,v)_K@f$, with real @f$b_j=j+1@f$. The vector field has zero
   * through three components, independently of the spatial dimension. Real and
   * complex fields use the same real transport coefficient. Empty ranges have
   * no matrix entries and are not assigned a fictitious entry read.
   *
   * ## Architecture
   *
   * Warm scopes retain scalar reference data at one integration order. Refresh
   * scopes alternate that order and the order four degrees higher, forcing the
   * existing owner's preparation stage without constructing a second cache.
   * Both include owned mapped quadrature and local contraction; neither measures
   * mesh traversal, global insertion or solving. Three base orders study local
   * sample-count scaling, not mesh-refinement convergence. Point formulas
   * always retain one sample and are not assigned artificial refinement rates.
   *
   * Independent entry-loop qualification and perturbation rejection precede
   * timing. Repeated binding must reproduce every retained entry exactly at
   * the final order. Empty-range qualification uses logical DOF counts instead.
   * Benchmarks are opt-in and are not CTest or publication timing gates.
   *
   * @tparam Scalar Real or complex field scalar.
   * @tparam Term Zero for mass, one for weighted mass, two for convection.
   */
  template <class Scalar, size_t Term>
  class P1QuadratureBinding
  {
    public:
      /// @brief Roundoff budget for the independently associated entry loop.
      static constexpr Real ReferenceTolerance = 1e-12;
      /// @brief Degree increment distinguishing the two refresh formulas.
      static constexpr size_t OrderIncrement = 4;

      /**
       * @brief Measures one-cell binding after numerical qualification.
       * @param state Component count and base integration order.
       * @param geometry Reference-cell family of the mapped cell.
       * @param refresh Whether integration orders alternate inside timing.
       */
      static void run(
        benchmark::State& state, Geometry::Polytope::Type geometry, bool refresh)
      {
        using namespace Geometry;
        using namespace Variational;
        const size_t components = state.range(0);
        const size_t baseOrder = state.range(1);
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        auto mesh = [&]() {
          if (dimension == 0)
          {
            auto builder = LocalMesh::Builder();
            builder.initialize(3).nodes(1).vertex({0, 0, 0});
            return builder.finalize();
          }
          const Array<size_t> shape = Array<size_t>::Constant(dimension, 2);
          return LocalMesh::UniformGrid(geometry, shape);
        }();
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        const auto cell = *mesh.getCell();
        P1<Math::SpatialVector<Scalar>> space(mesh, components);
        TrialFunction u(space);
        TestFunction v(space);
        const size_t ambient = mesh.getSpaceDimension();
        const VectorFunction direction(ambient, [ambient](const Point&) {
          Math::SpatialVector<Real> value(ambient);
          for (size_t axis = 0; axis < ambient; ++axis)
            value[axis] = Real(axis + 1);
          return value;
        });
        auto integral = [&]() {
          if constexpr (Term == 0)
            return Integral(u, v);
          else if constexpr (Term == 1)
            return Integral(Real(2) * u, v);
          else
            return Integral(Jacobian(u) * direction, v);
        }();
        static_assert(decltype(integral)::Specialized);
        const size_t ntr = u.getDOFs(cell), nte = v.getDOFs(cell);
        const std::array<size_t, 2> orders{
          baseOrder, refresh ? baseOrder + OrderIncrement : baseOrder};
        std::array<Math::Matrix<Scalar>, 2> initial;
        Real maxReferenceError = 0;
        for (size_t index = 0; index < orders.size(); ++index)
        {
          integral.setOrder(orders[index]);
          integral.setPolytope(cell);
          initial[index].resize(nte, ntr);
          for (size_t te = 0; te < nte; ++te)
          {
            for (size_t tr = 0; tr < ntr; ++tr)
              initial[index](te, tr) = integral.integrate(tr, te);
          }
          if (components == 0)
          {
            if (ntr != 0 || nte != 0)
            {
              p1BindingFailed = true;
              state.SkipWithError("Empty P1 range has nonempty DOFs");
              return;
            }
            continue;
          }
          QuadratureReference reference(integral);
          reference.assemble(cell, orders[index]);
          const auto& expected = reference.getOperator();
          const Real scale = std::max(Real(1), expected.norm());
          const auto matches = [&](const auto& actual) {
            return actual.rows() == expected.rows() && actual.cols() == expected.cols() &&
              (actual - expected).norm() <= ReferenceTolerance * scale;
          };
          auto wrong = expected;
          wrong(0, 0) += scale;
          if (!matches(initial[index]) || matches(wrong))
          {
            p1BindingFailed = true;
            state.SkipWithError("P1 local matrix or sensitivity check failed");
            return;
          }
          maxReferenceError =
            std::max(maxReferenceError, (initial[index] - expected).norm() / scale);
        }
        size_t last = 1;
        for (auto _ : state)
        {
          if (refresh)
            last = 1 - last;
          integral.setOrder(orders[last]);
          integral.setPolytope(cell);
          if (ntr != 0 && nte != 0)
            benchmark::DoNotOptimize(integral.integrate(0, 0));
          benchmark::ClobberMemory();
        }
        for (size_t te = 0; te < nte; ++te)
        {
          for (size_t tr = 0; tr < ntr; ++tr)
          {
            if (integral.integrate(tr, te) != initial[last](te, tr))
            {
              p1BindingFailed = true;
              state.SkipWithError("Repeated P1 binding changed an entry");
              return;
            }
          }
        }
        state.counters["components"] = components;
        state.counters["dimension"] = dimension;
        state.counters["ambient_dimension"] = ambient;
        state.counters["rows"] = nte;
        state.counters["columns"] = ntr;
        state.counters["entries"] = ntr * nte;
        state.counters["base_order"] = orders[0];
        state.counters["alternate_order"] = orders[1];
        state.counters["base_samples"] =
          QF::PolytopeQuadratureFormula::get(orders[0], geometry).getSize();
        state.counters["alternate_samples"] =
          QF::PolytopeQuadratureFormula::get(orders[1], geometry).getSize();
        state.counters["reference_error"] = maxReferenceError;
        state.counters["specialized"] = decltype(integral)::Specialized;
      }

      /// @brief Registers component and sample-count sweeps on every cell family.
      static void registerCases()
      {
        using Geometry::Polytope;
        for (const auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          const std::string geometryName = geometry == Polytope::Type::Point
            ? "Point"
            : std::string(Convergence::UniformGrid::getGeometryName(geometry));
          for (const bool refresh : {false, true})
          {
            const std::string name = std::string("P1/Binding/") +
              (Term == 0      ? "Mass/"
                  : Term == 1 ? "WeightedMass/"
                              : "Convection/") +
              (std::is_same_v<Scalar, Real> ? "Real/" : "Complex/") + geometryName +
              (refresh ? "/Refresh" : "/Warm");
            auto* registered = benchmark::RegisterBenchmark(name.c_str(),
              [geometry, refresh](auto& state) { run(state, geometry, refresh); });
            for (const int components : {0, 1, 2, 3})
            {
              for (const int order : {4, 8, 16})
                registered->Args({components, order});
            }
            registered->UseRealTime();
          }
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  P1QuadratureBinding<Real, 0>::registerCases();
  P1QuadratureBinding<Real, 1>::registerCases();
  P1QuadratureBinding<Real, 2>::registerCases();
  P1QuadratureBinding<Complex, 0>::registerCases();
  P1QuadratureBinding<Complex, 1>::registerCases();
  P1QuadratureBinding<Complex, 2>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext(
    "scope", "One-cell P1 binding; warm and order-refresh preparation");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || p1BindingFailed ? 1 : 0;
}
