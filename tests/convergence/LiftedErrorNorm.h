/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_CONVERGENCE_LIFTED_ERROR_NORM_H
#define RODIN_TESTS_CONVERGENCE_LIFTED_ERROR_NORM_H
#include <array>
#include <type_traits>
#include "Convergence.h"
namespace Rodin::Tests::Convergence
{
  /** @brief Real/complex scalar/vector exact-domain error decomposition.
   * @pre Meshes are full-dimensional and share ordered logical cell vertices;
   * the exact map is regular and orientation-preserving on the original box.
   * @par Architecture
   * The original unit-box mesh retains charts and exact logical cell indices.
   * The represented mesh is its geometry-only copy. Cached quadrature and
   * Geometry::Point provide both cell charts; the analytic map supplies the
   * exact-domain metric. For @f$x=\Phi(\xi)@f$ and
   * @f$x_h=\Phi_h(\xi)@f$, the lift is @f$u_h^\ell(x)=u_h(x_h)@f$.
   * Its physical gradient is
   * @f$D\Phi^{-T}D\Phi_h^T\nabla u_h(x_h)@f$.
   * For vector fields with component-row Jacobians, the lifted derivative is
   * @f$Ju_h(x_h)D\Phi_hD\Phi^{-1}@f$. An optional observer receives
   * the three derivative defects and the exact-domain quadrature weight.
   * Field and geometry defects add pointwise, but their norms do not.
   * Complex scalar and gradient errors use the full modulus squared; the
   * real geometric differential acts on both real and imaginary components.
   * MPI integrates owned original cells, reduces squared contributions, then
   * takes square roots. No inverse point location or coordinate matching is used.
   */
  class LiftedErrorNorm
  {
    public:
      struct Result
      {
          ErrorNorms field{0, 0}, geometry{0, 0}, total{0, 0};
      };
      template <class Mesh, class GF, class Data, class Map>
      static Result compute(const Mesh& reference, const Mesh& represented, const GF& uh,
        const Data& data, const Map& exactMap, size_t order)
      {
        return compute(
          reference, represented, uh, data, exactMap, order, [](const auto&, Real) {});
      }

      template <class Mesh, class GF, class Data, class Map, class Observer>
      static Result compute(const Mesh& reference, const Mesh& represented, const GF& uh,
        const Data& data, const Map& exactMap, size_t order, Observer&& observe)
      {
        using Space = typename FormLanguage::Traits<GF>::FESType;
        using Range = typename FormLanguage::Traits<Space>::RangeType;
        using Scalar = typename FormLanguage::Traits<Range>::ScalarType;
        constexpr bool IsVector = FormLanguage::IsVectorRange<Range>::Value;
        static_assert(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>);
        static_assert(IsVector || std::is_same_v<Range, Scalar>);
        using Derivative = std::conditional_t<IsVector, Math::SpatialMatrix<Scalar>,
          Math::SpatialVector<Scalar>>;
        std::array<Real, 6> squared{};
        const auto gradient = [&] {
          if constexpr (IsVector)
            return Variational::Jacobian(uh);
          else
            return Variational::Grad(uh);
        }();
        const auto exactGradient = [&data](const Math::SpatialPoint& x) {
          if constexpr (IsVector)
            return data.getJacobian(x);
          else
            return data.getGradient(x);
        };
        const size_t dim = reference.getSpaceDimension();
        for (auto cell = reference.getCell(); cell; ++cell)
        {
          if constexpr (requires { reference.getShard(); })
            if (!reference.getShard().isOwned(dim, cell->getIndex()))
              continue;
          const auto mapped = represented.getCell(cell->getIndex());
          assert(mapped->getGeometry() == cell->getGeometry());
          [[maybe_unused]] const auto vertices = cell->getVertices(),
                                      mappedVertices = mapped->getVertices();
          assert(vertices.size() == mappedVertices.size());
          for (size_t i = 0; i < vertices.size(); ++i)
            assert(vertices[i] == mappedVertices[i]);
          const auto& qf = QF::PolytopeQuadratureFormula::get(order, cell->getGeometry());
          const auto& originalQuadrature = cell->getQuadrature(qf);
          const auto& mappedQuadrature = mapped->getQuadrature(qf);
          for (size_t qp = 0; qp < originalQuadrature.getSize(); ++qp)
          {
            const auto& original = originalQuadrature.getPoint(qp);
            const auto& point = mappedQuadrature.getPoint(qp);
            const Variational::IntegrationPoint ip(point, &qf, qp);
            const auto exactPosition = exactMap(original.getPhysicalCoordinates());
            const auto exactJacobian =
              exactMap.getJacobian(original.getPhysicalCoordinates());
            const Math::SpatialMatrix<Real> mappedJacobian =
              point.getJacobian() * original.getJacobianInverse();
            const Math::SpatialMatrix<Real> lift =
              exactJacobian.inverse().transpose() * mappedJacobian.transpose();
            const Real weight =
              qf.getWeight(qp) * original.getDistortion() * exactJacobian.determinant();
            assert(original.getDistortion() > 0 && exactJacobian.determinant() > 0);
            const Range value = uh(ip);
            const Range representedValue =
              data.getSolution(point.getPhysicalCoordinates());
            const Range exactValue = data.getSolution(exactPosition);
            const Derivative derivative = gradient(ip);
            const auto representedDerivative =
              exactGradient(point.getPhysicalCoordinates());
            const auto exactDerivative = exactGradient(exactPosition);
            const auto transform = [&lift](const Derivative& derivative) -> Derivative {
              if constexpr (IsVector)
                return derivative * lift.transpose();
              else
                return lift * derivative;
            };
            const std::array<Range, 3> values{value - representedValue,
              representedValue - exactValue, value - exactValue};
            const std::array<Derivative, 3> derivatives{
              transform(derivative - representedDerivative),
              transform(representedDerivative) - exactDerivative,
              transform(derivative) - exactDerivative};
            for (size_t i = 0; i < values.size(); ++i)
            {
              squared[2 * i] += weight * ErrorNorm::squaredMagnitude(values[i]);
              squared[2 * i + 1] += weight * ErrorNorm::squaredMagnitude(derivatives[i]);
            }
            observe(derivatives, weight);
          }
        }
#ifdef RODIN_USE_MPI
        if constexpr (requires { reference.getShard(); })
          for (Real& value : squared)
            value = boost::mpi::all_reduce(
              reference.getContext().getCommunicator(), value, std::plus<Real>());
#endif
        return {{std::sqrt(squared[0]), std::sqrt(squared[1])},
          {std::sqrt(squared[2]), std::sqrt(squared[3])},
          {std::sqrt(squared[4]), std::sqrt(squared[5])}};
      }
  };
}
#endif
