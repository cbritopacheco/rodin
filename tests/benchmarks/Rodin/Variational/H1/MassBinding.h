/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODIN_VARIATIONAL_H1_MASSBINDING_H
#define RODIN_RODIN_VARIATIONAL_H1_MASSBINDING_H

/**
 * @file
 * @brief Shared fixture for H1 mass and scalar operand-order measurements.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <string>
#include <type_traits>
#include "Rodin/Variational.h"
#include "../../../../QuadratureReference.h"
#include "../../../../convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  enum class MassRange
  {
    Scalar,
    Vector,
    Matrix
  };

  /// @brief Reports untimed qualification failures to the executable entry point.
  inline bool massBindingFailed = false;

  /**
   * @brief Measures scalar-weighted H1 mass preparation and component contraction.
   *
   * ## Mathematical foundation
   *
   * For scalar weight @f$\rho@f$, the discrete component entries satisfy
   * @f[
   * A_{bm+e,am+c}=\delta_{ec}S_{ba},\qquad
   * S_{ba}=\sum_q w_q D_q\rho(x_q)\phi_a(x_q)\overline{\phi_b(x_q)}.
   * @f]
   * Matrix fields use their local row-major component convention and the
   * Frobenius product. Unit mass uses @f$\rho=1@f$. Weighted cases use
   * @f$\rho=s(2+x_0)@f$ for real fields and
   * @f$\rho=s(2+x_0+i/2)@f$ for complex fields. This is not a Hermitian
   * symmetry claim when the coefficient is complex.
   *
   * ## Measurement boundary
   *
   * One-cell public binding includes mapped quadrature, scalar reference-data
   * lookup and complete local matrix preparation. Construction, independent
   * qualification and final repeated-entry checks are untimed. Scatter, load
   * assembly and projection solves belong to their owning consumer modules.
   *
   * The original expression-entry loop qualifies only the scalar block. Every
   * component entry is then checked exactly against the retained first block,
   * including structural zeros. A perturbed scalar block must be rejected.
   * The reference storage is therefore scalar-sized, not a second full
   * component matrix. Doubling the borrowed coefficient scale after timing
   * checks replacement rather than accumulation of stale operator values.
   * Explicit orders define discrete operators; underintegrated samples are
   * not claims about polynomial integration exactness or PDE convergence.
   *
   * Scalar weights are measured in both operand orders. Although
   * @f$\rho u=u\rho@f$ algebraically, their expression types and local
   * quadrature dispatch may differ. The dispatch counter records this
   * implementation property without requiring one particular implementation.
   * General matrix coefficients are not commuted.
   *
   * @tparam K Scalar polynomial degree.
   * @tparam Scalar Real or complex field scalar.
   * @tparam Range Scalar, vector or matrix field range.
   * @tparam Weighted Whether the scalar coefficient is nonconstant.
   * @tparam Right Whether a scalar weight multiplies the trial field on the right.
   */
  template <size_t K, class Scalar, MassRange Range, bool Weighted, bool Right = false>
  class H1MassBinding
  {
    public:
      static_assert(!Right || Weighted, "A right operand requires a scalar weight");

      /// @brief Budget for independently associated scalar expression quadrature.
      static constexpr Real ReferenceBudget = 1e-12;

      /**
       * @brief Measures repeated binding of one selected reference-cell family.
       * @param state Field rows, columns and requested quadrature order.
       * @param geometry Native reference-cell family, including points.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
      {
        using namespace Geometry;
        using namespace Variational;
        const size_t rows = state.range(0), cols = state.range(1);
        const size_t components = rows * cols, order = state.range(2);
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
        auto space = [&]() {
          if constexpr (Range == MassRange::Scalar)
            return H1<K, Scalar>(std::integral_constant<size_t, K>{}, mesh);
          else if constexpr (Range == MassRange::Vector)
          {
            return H1<K, Math::SpatialVector<Scalar>>(
              std::integral_constant<size_t, K>{}, mesh, rows);
          }
          else
          {
            return H1<K, Math::SpatialMatrix<Scalar>>(
              std::integral_constant<size_t, K>{}, mesh, rows, cols);
          }
        }();
        H1<K, Scalar> scalarSpace(std::integral_constant<size_t, K>{}, mesh);
        TrialFunction u(space);
        TestFunction v(space);
        TrialFunction scalarU(scalarSpace);
        TestFunction scalarV(scalarSpace);
        Real scale = 1;
        const auto coefficient = [&]() {
          if constexpr (std::is_same_v<Scalar, Complex>)
          {
            return ComplexFunction([&](const Point& point) {
              const Real x = point.getDimension() == 0 ? Real(0) : point.x();
              return scale * Complex(Real(2) + x, Real(0.5));
            });
          }
          else
          {
            return RealFunction([&](const Point& point) {
              const Real x = point.getDimension() == 0 ? Real(0) : point.x();
              return scale * (Real(2) + x);
            });
          }
        }();
        auto integral = [&]() {
          if constexpr (Weighted && Right)
            return Integral(u * coefficient, v);
          else if constexpr (Weighted)
            return Integral(coefficient * u, v);
          else
            return Integral(u, v);
        }();
        auto scalarIntegral = [&]() {
          if constexpr (Weighted && Right)
            return Integral(scalarU * coefficient, scalarV);
          else if constexpr (Weighted)
            return Integral(coefficient * scalarU, scalarV);
          else
            return Integral(scalarU, scalarV);
        }();
        if constexpr (!Right)
        {
          static_assert(decltype(integral)::Specialized);
        }
        integral.setOrder(order);
        integral.setPolytope(cell);
        const size_t nb = scalarU.getDOFs(cell), n = u.getDOFs(cell);
        if (n != components * nb || v.getDOFs(cell) != n || n == 0)
        {
          massBindingFailed = true;
          state.SkipWithError("Local component counts do not match the scalar basis");
          return;
        }
        Math::Matrix<Scalar> initial(nb, nb);
        for (size_t b = 0; b < nb; ++b)
        {
          for (size_t a = 0; a < nb; ++a)
            initial(b, a) = integral.integrate(a * components, b * components);
        }
        QuadratureReference reference(scalarIntegral);
        reference.assemble(cell, order);
        const auto& expected = reference.getOperator();
        const Real normalization = std::max(Real(1), expected.norm());
        const Real deviation = (initial - expected).norm() / normalization;
        auto wrong = initial;
        wrong(0, 0) += normalization;
        if (!(deviation <= ReferenceBudget) ||
          !((wrong - expected).norm() / normalization > ReferenceBudget))
        {
          massBindingFailed = true;
          state.SkipWithError("Scalar mass block or perturbation rejection failed");
          return;
        }
        const auto matches = [&](Scalar factor) {
          for (size_t b = 0; b < nb; ++b)
          {
            for (size_t a = 0; a < nb; ++a)
            {
              for (size_t e = 0; e < components; ++e)
              {
                for (size_t c = 0; c < components; ++c)
                {
                  const Scalar value =
                    integral.integrate(a * components + c, b * components + e);
                  const Scalar target = c == e ? factor * initial(b, a) : Scalar(0);
                  if (!std::isfinite(std::abs(value)) || value != target)
                    return false;
                }
              }
            }
          }
          return true;
        };
        if (!matches(Scalar(1)))
        {
          massBindingFailed = true;
          state.SkipWithError("Initial component blocks or structural zeros differ");
          return;
        }
        for (auto _ : state)
        {
          integral.setPolytope(cell);
          benchmark::DoNotOptimize(integral.integrate(0, 0));
          benchmark::ClobberMemory();
        }
        if (!matches(Scalar(1)))
        {
          massBindingFailed = true;
          state.SkipWithError("Repeated binding changed a component entry");
          return;
        }
        if constexpr (Weighted)
        {
          scale = 2;
          integral.setPolytope(cell);
          if (!matches(Scalar(2)))
          {
            massBindingFailed = true;
            state.SkipWithError("Coefficient update did not replace component blocks");
            return;
          }
        }
        const size_t nq = QF::PolytopeQuadratureFormula::get(order, geometry).getSize();
        state.counters["degree"] = K;
        state.counters["dimension"] = dimension;
        state.counters["ambient_dimension"] = mesh.getSpaceDimension();
        state.counters["components"] = components;
        state.counters["field_rows"] = rows;
        state.counters["field_columns"] = cols;
        state.counters["cells"] = 1;
        state.counters["fixture_cells"] = mesh.getPolytopeCount(dimension);
        state.counters["global_dofs"] = space.getSize();
        state.counters["geometry_factor_order"] =
          cell.getTransformation().getFactorOrder();
        state.counters["scalar_basis"] = nb;
        state.counters["local_basis"] = n;
        state.counters["local_entries"] = n * n;
        state.counters["quadrature_order"] = order;
        state.counters["quadrature_points"] = nq;
        state.counters["scalar_pair_sample_proxy"] = nq * nb * nb;
        state.counters["component_pair_sample_proxy"] = nq * nb * nb * components;
        state.counters["dense_entry_sample_proxy"] = nq * n * n;
        state.counters["reference_error"] = deviation;
        state.counters["specialized_kernel"] = decltype(integral)::Specialized;
      }

      /**
       * @brief Registers native geometries at fixed and degree-scaled orders.
       *
       * The original orders @f$r\in\{0,4,8\}@f$ are retained. The additional
       * orders @f$r\in\{2K+1,2K+3,2K+5\}@f$ scale with the basis degree.
       * On affine polynomial elements, @f$2K+1@f$ accounts for the product
       * of two degree-@f$K@f$ bases and the linear weight; the bound is
       * coordinatewise for tensor-product bases and conservative for unit
       * mass. These orders do not assert exactness for rational pyramidal
       * bases. Point cases remain structural checks, not refinement studies.
       *
       * @param prefix Module-owned benchmark name prefix.
       */
      static void registerCases(const std::string& prefix = "H1/MassBinding/")
      {
        using Geometry::Polytope;
        constexpr int DegreeOrder = 2 * static_cast<int>(K) + 1;
        for (const auto geometry :
          {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge})
        {
          const std::string name = prefix +
            std::string(Right ? "RightScalarWeighted/"
                : Weighted    ? "ScalarWeighted/"
                              : "Unit/") +
            (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
            (Range == MassRange::Scalar      ? "Scalar/"
                : Range == MassRange::Vector ? "Vector/"
                                             : "Matrix/") +
            "P" + std::to_string(K) + "/" +
            (geometry == Polytope::Type::Point
                ? std::string("Point")
                : std::string(Convergence::UniformGrid::getGeometryName(geometry))) +
            "/Warm";
          auto* registered = benchmark::RegisterBenchmark(
            name.c_str(), [geometry](auto& state) { run(state, geometry); });
          for (const int order : {0, 4, 8, DegreeOrder, DegreeOrder + 2, DegreeOrder + 4})
          {
            if constexpr (Range == MassRange::Scalar)
              registered->Args({1, 1, order});
            else if constexpr (Range == MassRange::Vector)
            {
              for (const int components : {1, 2, 3})
                registered->Args({components, 1, order});
            }
            else
            {
              registered->Args({1, 2, order});
              registered->Args({2, 3, order});
            }
          }
          registered->UseRealTime();
        }
      }
  };
}

#endif // RODIN_RODIN_VARIATIONAL_H1_MASSBINDING_H
