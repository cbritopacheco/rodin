/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Public grid-function differential evaluation across binding protocols.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <string>
#include <type_traits>
#include <vector>
#include "convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  bool fieldDerivativeFailed = false;

  enum class FieldDerivativeProtocol
  {
    Point,
    Provenance,
    Associated
  };

  /**
   * @brief Measures field derivatives rather than individual shape derivatives.
   *
   * ## Mathematical foundation
   *
   * Degree-one scopes interpolate the affine field
   * @f$f(x)=\eta(1+\sum_j(j+1)x_j)@f$. Higher-order scopes use
   * @f$f(x)=\eta(1+\sum_j(j+1)(x_j+x_j^2))@f$, whose varying derivative
   * checks sample dependence. Vector fields have components @f$u_c=(c+1)f@f$.
   * Here @f$\eta=1@f$ or @f$1+i/2@f$. Interpolation applies the space's DOF
   * functionals; no coefficient layout or nodal-value interpretation is assumed.
   * Analytic gradients/Jacobians, perturbation rejection and coefficient scaling
   * @f$1\to2\to1@f$ are checked outside timing. Repeating the timed protocol must
   * reproduce every retained entry exactly; different protocols are not required
   * to have identical floating-point association.
   *
   * ## Architecture
   *
   * The same retained cell and mapped samples are used by direct geometric-point,
   * provenance-only integration-point and owner/index-associated protocols.
   * Geometry, scalar reference tables, interpolation and output buffers are
   * prepared before timing. Public differential evaluation and output assignment
   * remain timed. This isolates the source fallback consumed by error integration
   * without changing convergence quadrature, meshes, norms or solver requirements.
   * Native P1 remains a separate control: an associated binding does not promise
   * that every differential implementation consumes reference tables.
   * These are local, opt-in component benchmarks, not CTest or CI timing gates.
   *
   * @tparam K H1 degree, or one for native P1.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Vector Whether to evaluate a vector-field Jacobian instead of a gradient.
   * @tparam NativeP1 Whether to use the distinct native P1 implementation.
   */
  template <size_t K, class Scalar, bool Vector, bool NativeP1 = false>
  class FieldDerivative
  {
    public:
      /// @brief Relative analytic derivative budget allowing basis roundoff.
      static constexpr Real Tolerance = 1e-10;

      /**
       * @brief Sweeps a prepared cell through the selected public differential path.
       * @param state Quadrature order and iteration protocol.
       * @param geometry Positive-dimensional reference-cell family.
       * @param protocol Point or integration-point binding contract.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        FieldDerivativeProtocol protocol)
      {
        using namespace Variational;
        using Range = std::conditional_t<Vector, Math::SpatialVector<Scalar>, Scalar>;
        using Space = std::conditional_t<NativeP1, P1<Range>, H1<K, Range>>;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
        const auto cell = *mesh.getCell();
        const size_t dim = mesh.getDimension();
        auto space = [&]() {
          if constexpr (NativeP1 && Vector)
            return Space(mesh, dim);
          else if constexpr (NativeP1)
            return Space(mesh);
          else if constexpr (Vector)
            return Space(std::integral_constant<size_t, K>{}, mesh, dim);
          else
            return Space(std::integral_constant<size_t, K>{}, mesh);
        }();
        const Scalar eta = []() {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 0.5);
          else
            return Scalar(1);
        }();
        Real scale = 1;
        const auto scalarValue = [&](const Geometry::Point& point) {
          Scalar value = eta;
          const auto& x = point.getPhysicalCoordinates();
          for (size_t j = 0; j < dim; ++j)
          {
            const Real coordinate = x(j);
            value += eta * Real(j + 1) *
              (coordinate + (K > 1 ? coordinate * coordinate : Real(0)));
          }
          return scale * value;
        };
        const auto source = [&]() {
          if constexpr (Vector)
          {
            return VectorFunction(dim, [&](const Geometry::Point& point) {
              Range value(dim);
              const Scalar scalar = scalarValue(point);
              for (size_t c = 0; c < dim; ++c)
                value(c) = Scalar(c + 1) * scalar;
              return value;
            });
          }
          else if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(scalarValue);
          else
            return RealFunction(scalarValue);
        }();
        GridFunction field(space);
        field = source;
        auto differential = [&]() {
          if constexpr (Vector)
            return Jacobian(field);
          else
            return Grad(field);
        }();
        const auto& formula = QF::QuadratureFormula::get(state.range(0), geometry);
        const Geometry::PolytopeQuadrature quadrature(cell, formula);
        using Value =
          std::remove_cvref_t<decltype(differential.getValue(quadrature.getPoint(0)))>;
        using ScalarElement =
          std::conditional_t<NativeP1, P1Element<Scalar>, H1Element<K, Scalar>>;
        const ScalarElement scalarElement(geometry);
        scalarElement.getTabulation(formula);
        std::vector<IntegrationPoint> samples;
        samples.reserve(quadrature.getSize());
        for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
        {
          const auto& point = quadrature.getPoint(qp);
          point.getPhysicalCoordinates();
          point.getJacobianInverse();
          if (protocol == FieldDerivativeProtocol::Provenance)
            samples.emplace_back(point, &formula, qp);
          else
            samples.emplace_back(quadrature, qp);
          const bool associated =
            samples.back().getReferenceQuadratureFormula() != nullptr;
          if (associated != (protocol != FieldDerivativeProtocol::Provenance))
          {
            fieldDerivativeFailed = true;
            state.SkipWithError("Integration-point association contract changed");
            return;
          }
        }
        std::vector<Value> measured(samples.size());
        const auto evaluate = [&]() {
          for (size_t qp = 0; qp < samples.size(); ++qp)
          {
            measured[qp] = protocol == FieldDerivativeProtocol::Point
              ? differential.getValue(samples[qp].getPoint())
              : differential.getValue(samples[qp]);
          }
        };
        Real maximumRelativeError = 0;
        const auto valid = [&]() {
          maximumRelativeError = 0;
          for (size_t qp = 0; qp < samples.size(); ++qp)
          {
            const auto& x = samples[qp].getPoint().getPhysicalCoordinates();
            for (size_t c = 0; c < (Vector ? dim : 1); ++c)
            {
              for (size_t j = 0; j < dim; ++j)
              {
                const Scalar expected = scale * eta * Real(c + 1) * Real(j + 1) *
                  (K > 1 ? Real(1) + Real(2) * x(j) : Real(1));
                const Scalar actual = [&]() {
                  if constexpr (Vector)
                    return measured[qp](c, j);
                  else
                    return measured[qp](j);
                }();
                const Real error =
                  std::abs(actual - expected) / std::max(Real(1), std::abs(expected));
                maximumRelativeError = std::max(maximumRelativeError, error);
                if (!std::isfinite(std::abs(actual)) || !std::isfinite(error) ||
                  error > Tolerance)
                  return false;
              }
            }
          }
          return true;
        };
        evaluate();
        const auto initial = measured;
        bool qualified = valid();
        if constexpr (Vector)
          measured.front()(0, 0) += Scalar(1);
        else
          measured.front()(0) += Scalar(1);
        qualified = qualified && !valid();
        scale = 2;
        field = source;
        evaluate();
        qualified = qualified && valid();
        scale = 1;
        field = source;
        evaluate();
        qualified = qualified && valid();
        if (!qualified)
        {
          fieldDerivativeFailed = true;
          state.SkipWithError(
            "Analytic differential or coefficient-refresh oracle failed");
          return;
        }
        for (auto _ : state)
        {
          evaluate();
          benchmark::DoNotOptimize(measured.data());
          benchmark::ClobberMemory();
        }
        bool unchanged = valid();
        Scalar checksum = 0;
        for (size_t qp = 0; qp < samples.size(); ++qp)
        {
          for (size_t c = 0; c < (Vector ? dim : 1); ++c)
          {
            for (size_t j = 0; j < dim; ++j)
            {
              if constexpr (Vector)
              {
                unchanged = unchanged && measured[qp](c, j) == initial[qp](c, j);
                checksum += measured[qp](c, j);
              }
              else
              {
                unchanged = unchanged && measured[qp](j) == initial[qp](j);
                checksum += measured[qp](j);
              }
            }
          }
        }
        if (!unchanged)
        {
          fieldDerivativeFailed = true;
          state.SkipWithError(
            "Repeated differential evaluation changed a retained entry");
        }
        state.counters["cells"] = 1;
        state.counters["dimension"] = dim;
        state.counters["degree"] = K;
        state.counters["global_dofs"] = space.getSize();
        state.counters["scalar_basis"] = scalarElement.getCount();
        state.counters["quadrature_order"] = state.range(0);
        state.counters["point_evaluations"] = samples.size();
        state.counters["scalar_output_entries"] =
          samples.size() * dim * (Vector ? dim : 1);
        state.counters["relative_derivative_error"] = maximumRelativeError;
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          state.counters["checksum_real"] = checksum.real();
          state.counters["checksum_imag"] = checksum.imag();
        }
        else
        {
          state.counters["checksum_real"] = checksum;
          state.counters["checksum_imag"] = 0;
        }
      }

      /// @brief Registers each geometry and binding independently at three quadrature orders.
      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral,
               G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
        {
          for (const auto protocol :
            {FieldDerivativeProtocol::Point, FieldDerivativeProtocol::Provenance,
              FieldDerivativeProtocol::Associated})
          {
            const std::string binding = protocol == FieldDerivativeProtocol::Point
              ? "Point"
              : protocol == FieldDerivativeProtocol::Provenance ? "Provenance"
                                                                : "Associated";
            const std::string name = std::string("FieldDerivative/") +
              (Vector ? "Jacobian/" : "Gradient/") + binding +
              (NativeP1 ? "/P1/" : "/H1/P" + std::to_string(K) + "/") +
              (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
              std::string(Convergence::UniformGrid::getGeometryName(geometry));
            benchmark::RegisterBenchmark(name.c_str(),
              [geometry, protocol](auto& state) { run(state, geometry, protocol); })
              ->Arg(4)
              ->Arg(8)
              ->Arg(16)
              ->UseRealTime();
          }
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  const auto registerScalar = []<class Scalar>() {
    FieldDerivative<1, Scalar, false, true>::registerCases();
    FieldDerivative<1, Scalar, true, true>::registerCases();
    FieldDerivative<1, Scalar, false>::registerCases();
    FieldDerivative<1, Scalar, true>::registerCases();
    FieldDerivative<2, Scalar, false>::registerCases();
    FieldDerivative<2, Scalar, true>::registerCases();
    FieldDerivative<3, Scalar, false>::registerCases();
    FieldDerivative<3, Scalar, true>::registerCases();
  };
  registerScalar.template operator()<Real>();
  registerScalar.template operator()<Complex>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "One-cell public field differential; warmed maps/tables");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || fieldDerivativeFailed ? 1 : 0;
}
