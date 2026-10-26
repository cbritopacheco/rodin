/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>
#include <array>
#include <type_traits>
#include "Rodin/Variational.h"
#include "../../../../QuadratureReference.h"
#include "../../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  bool convectionFailed = false;

  /**
   * @brief Measures H1 directional-gradient preparation and local contraction.
   *
   * ## Mathematical foundation
   *
   * The measured form is @f$a_b(u,v)=\int_K(Ju)b\cdot\overline v\,dx@f$.
   * Fields are real or complex; complex fields admit either real or genuinely
   * complex directions. Real fields use real directions, as required by the
   * current shape-Dot range constraint. An
   * independent expression-entry quadrature loop and a deliberate perturbation
   * qualify every matrix before timing. The final matrix must repeat exactly.
   * A defective baseline is not timed for the unsupported coefficient semantics.
   *
   * ## Architecture
   *
   * Binding includes owned mapped quadrature, scalar tabulation lookup, physical
   * gradients and local contraction. Traversal, insertion and solving are absent.
   * Warm scopes reuse one order; Refresh alternates two orders on the same owner.
   * Degree, components, quadrature samples and both physical dimensions are
   * reported separately. This component executable is opt-in and outside CTest.
   *
   * @tparam Scalar Field scalar type.
   * @tparam Degree Scalar polynomial degree.
   * @tparam ComplexDirection Whether the direction has a nonzero imaginary part.
   */
  template <class Scalar, size_t Degree, bool ComplexDirection>
    requires(!ComplexDirection || std::is_same_v<Scalar, Complex>)
  class H1Convection
  {
    public:
      /// @brief Budget for differently associated generic entry-loop quadrature.
      static constexpr Real ReferenceBudget = 1e-12;

      /// @brief Integration-order increment for a reference-data refresh.
      static constexpr size_t OrderIncrement = 4;

      /**
       * @brief Measures one-cell binding after complete local qualification.
       * @param state Component count and base quadrature order.
       * @param geometry Native reference-cell family.
       * @param refresh Whether to alternate quadrature orders inside timing.
       */
      static void run(benchmark::State& state, Polytope::Type geometry, bool refresh)
      {
        const size_t components = state.range(0), order = state.range(1);
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        auto mesh = [&]() {
          if (dimension == 0)
          {
            return LocalMesh::Builder()
              .initialize(3)
              .nodes(1)
              .vertex({0, 0, 0})
              .finalize();
          }
          return LocalMesh::UniformGrid(geometry, Array<size_t>::Constant(dimension, 2));
        }();
        for (size_t d = 1; d <= dimension; ++d)
        {
          for (size_t lower = 0; lower < d; ++lower)
            mesh.getConnectivity().compute(d, lower);
        }
        const auto cell = *mesh.getCell();
        H1<Degree, Math::SpatialVector<Scalar>> space(
          std::integral_constant<size_t, Degree>{}, mesh, components);
        TrialFunction u(space);
        TestFunction v(space);
        const size_t ambient = mesh.getSpaceDimension();
        using DirectionScalar = std::conditional_t<ComplexDirection, Complex, Real>;
        const VectorFunction direction(ambient, [ambient](const Point&) {
          Math::SpatialVector<DirectionScalar> result(ambient);
          for (size_t j = 0; j < ambient; ++j)
          {
            if constexpr (ComplexDirection)
              result[j] = Complex(Real(j + 1), Real(j + 2));
            else
              result[j] = Real(j + 1);
          }
          return result;
        });
        const auto integrand = Dot(Jacobian(u) * direction, v);
        Integral integral(integrand);
        static_assert(decltype(integral)::Specialized);
        const size_t n = u.getDOFs(cell);
        const std::array<size_t, 2> orders{
          order, refresh ? order + OrderIncrement : order};
        std::array<Math::Matrix<Scalar>, 2> initial;
        Real error = 0;
        for (size_t index = 0; index < orders.size(); ++index)
        {
          integral.setOrder(orders[index]);
          integral.setPolytope(cell);
          initial[index].resize(n, n);
          for (size_t te = 0; te < n; ++te)
          {
            for (size_t tr = 0; tr < n; ++tr)
              initial[index](te, tr) = integral.integrate(tr, te);
          }
          if (n == 0)
            continue;
          QuadratureReference<decltype(integral)> reference(integrand);
          reference.assemble(cell, orders[index]);
          const auto& expected = reference.getOperator();
          const Real scale = std::max(Real(1), expected.norm());
          const Real deviation = (initial[index] - expected).norm() / scale;
          auto wrong = initial[index];
          wrong(0, 0) += scale;
          if (!(deviation <= ReferenceBudget) ||
            !((wrong - expected).norm() / scale > ReferenceBudget))
          {
            convectionFailed = true;
            state.SkipWithError("Directional matrix qualification failed");
            return;
          }
          error = std::max(error, deviation);
        }
        size_t last = 1;
        for (auto _ : state)
        {
          if (refresh)
            last = 1 - last;
          integral.setOrder(orders[last]);
          integral.setPolytope(cell);
          if (n != 0)
            benchmark::DoNotOptimize(integral.integrate(0, 0));
          benchmark::ClobberMemory();
        }
        for (size_t te = 0; te < n; ++te)
        {
          for (size_t tr = 0; tr < n; ++tr)
          {
            if (integral.integrate(tr, te) != initial[last](te, tr))
            {
              convectionFailed = true;
              state.SkipWithError("Repeated directional binding changed an entry");
              return;
            }
          }
        }
        state.counters["degree"] = Degree;
        state.counters["field_complex"] = std::is_same_v<Scalar, Complex>;
        state.counters["components"] = components;
        state.counters["dimension"] = dimension;
        state.counters["ambient_dimension"] = ambient;
        state.counters["local_dofs"] = n;
        state.counters["rows"] = n;
        state.counters["columns"] = n;
        state.counters["entries"] = n * n;
        state.counters["base_order"] = orders[0];
        state.counters["alternate_order"] = orders[1];
        state.counters["base_samples"] =
          QF::QuadratureFormula::get(orders[0], geometry).getSize();
        state.counters["alternate_samples"] =
          QF::QuadratureFormula::get(orders[1], geometry).getSize();
        state.counters["reference_error"] = error;
      }

      /// @brief Registers native geometry, component and sample-count sweeps.
      static void registerCases()
      {
        for (const auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (const bool refresh : {false, true})
          {
            const std::string name = "H1/Convection/" +
              std::string(std::is_same_v<Scalar, Real> ? "RealField/" : "") + "P" +
              std::to_string(Degree) +
              (ComplexDirection ? "/ComplexDirection/" : "/RealDirection/") +
              (geometry == Polytope::Type::Point
                  ? std::string("Point")
                  : std::string(Convergence::UniformGrid::getGeometryName(geometry))) +
              (refresh ? "/Refresh" : "/Warm");
            auto* registered = benchmark::RegisterBenchmark(name.c_str(),
              [geometry, refresh](auto& state) { run(state, geometry, refresh); });
            for (const int components : {0, 1, 2, 3})
            {
              for (const int order : {4, 8})
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
  using namespace Rodin::Tests::Benchmarks;
  H1Convection<Real, 1, false>::registerCases();
  H1Convection<Real, 2, false>::registerCases();
  H1Convection<Real, 3, false>::registerCases();
  H1Convection<Complex, 1, false>::registerCases();
  H1Convection<Complex, 1, true>::registerCases();
  H1Convection<Complex, 2, false>::registerCases();
  H1Convection<Complex, 2, true>::registerCases();
  H1Convection<Complex, 3, false>::registerCases();
  H1Convection<Complex, 3, true>::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("scope",
    "H1 directional-gradient binding; real/complex fields and compatible directions");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || convectionFailed ? 1 : 0;
}
