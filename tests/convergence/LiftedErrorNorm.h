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
  /** @brief Real/complex scalar field, geometry and total exact-domain error norms.
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
        std::array<Real, 6> squared{};
        const auto gradient = Variational::Grad(uh);
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
            using Space = typename FormLanguage::Traits<GF>::FESType;
            using Scalar = typename FormLanguage::Traits<Space>::RangeType;
            static_assert(
              std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>);
            const Scalar value = uh(ip);
            const Scalar representedValue =
              data.getSolution(point.getPhysicalCoordinates());
            const Scalar exactValue = data.getSolution(exactPosition);
            const Math::SpatialVector<Scalar> derivative = gradient(ip);
            const auto representedDerivative =
              data.getGradient(point.getPhysicalCoordinates());
            const auto exactDerivative = data.getGradient(exactPosition);
            const std::array<Scalar, 3> values{value - representedValue,
              representedValue - exactValue, value - exactValue};
            const std::array<Math::SpatialVector<Scalar>, 3> derivatives{
              lift * (derivative - representedDerivative),
              lift * representedDerivative - exactDerivative,
              lift * derivative - exactDerivative};
            for (size_t i = 0; i < values.size(); ++i)
            {
              squared[2 * i] += weight * ErrorNorm::squaredMagnitude(values[i]);
              squared[2 * i + 1] += weight * ErrorNorm::squaredMagnitude(derivatives[i]);
            }
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
