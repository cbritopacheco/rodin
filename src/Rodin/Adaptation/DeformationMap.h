/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_DEFORMATIONMAP_H
#define RODIN_ADAPTATION_DEFORMATIONMAP_H

#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Adaptation
{
  /**
   * @brief The deformation @f$\Phi(x)=x+u(x)@f$ induced by a displacement,
   * evaluated with memoisation.
   *
   * Locating @f$\Phi(x)@f$ in the mesh requires inverting a polytope
   * transformation and, on a face, searching the adjacent cells. That cost is
   * invariant with respect to the basis index, but the form language evaluates
   * coefficients through @c getBasis(local) and so may request the same moved
   * point once per basis function. This class makes the request idempotent: the
   * last computed point is cached and returned directly when the integration
   * point has not changed.
   *
   * A single-entry cache suffices because repeated requests are always
   * consecutive — they come from the innermost basis loop of a single
   * quadrature point — so the memory cost is @f$O(1)@f$ regardless of mesh size
   * and the hit rate is the full redundancy.
   *
   * @tparam Displacement Displacement grid function type.
   */
  template <class Displacement, class LocatorType>
  class DeformationMap
  {
    public:
      /**
       * @brief Constructs the deformation induced by @p u.
       * @param u Displacement field.
       * @param locator Bounding-volume locator bound to the displaced mesh.
       */
      DeformationMap(const Displacement& u, const LocatorType& locator)
        : m_u(u),
          m_locator(locator)
      {}

      /**
       * @brief The underlying displacement field.
       * @returns The displacement.
       */
      const Displacement& getDisplacement() const
      {
        return m_u.get();
      }

      /**
       * @brief The displacement @f$u(x)@f$ at an integration point, as a
       * spatial vector.
       * @param pt Point at which the operation is evaluated.
       * @param ip Integration point at which the expression is evaluated.
       * @returns The displacement value.
       */
      Math::SpatialVector<Real> getDisplacementValue(
        const Geometry::Point& pt, const Variational::IntegrationPoint& ip) const
      {
        const auto uq =
          ip.getQuadratureFormula() ? m_u.get().getValue(ip) : m_u.get().getValue(pt);
        Math::SpatialVector<Real> displacement(pt.getPolytope().getMesh().getDimension());
        displacement.setZero();
        for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
          displacement(static_cast<Eigen::Index>(r)) = uq(r);
        return displacement;
      }

      /**
       * @brief The moved point @f$\Phi(x)=x+u(x)@f$, located in the mesh.
       *
       * Consecutive calls may reuse a located point only for a certified
       * reference sample with unchanged mesh, entity, formula version and
       * geometry state. Sampling provenance alone does not authorize reuse.
       * Displacement or locator updates still require @ref invalidate.
       *
       * @param ip Integration point on the reference configuration; its
       * reference association and local geometry identify the source sample.
       * @returns Located point in the deformed configuration.
       */
      const Geometry::Point& getMovedPoint(const Variational::IntegrationPoint& ip) const
      {
        const auto& pt = ip.getPoint();
        const auto& polytope = pt.getPolytope();
        const auto* formula = ip.getReferenceQuadratureFormula();
        const auto identity = polytope.getTransformation().getCacheIdentity();
        const Key key{&polytope.getMesh(), polytope.getDimension(), polytope.getIndex(),
          formula, formula ? formula->getCacheIdentity() : 0, ip.getIndex(), identity};
        if (!formula || !identity || !m_moved || !m_key || !(*m_key == key))
        {
            // Geometry::Point is constructible but not assignable, so the slot
            // is re-emplaced rather than overwritten.
          m_moved.emplace(makeMovedPoint(pt,
            pt.getPhysicalCoordinates() + getDisplacementValue(pt, ip), m_locator.get()));
          m_key = key;
        }
        return *m_moved;
      }

      /// @brief Drops the cached point, forcing the next request to relocate.
      void invalidate() const
      {
        m_moved.reset();
        m_key.reset();
      }

    private:
      /**
       * @brief Resolves a point at moved physical coordinates.
       * @param source Source point retaining the reference cell and coordinates.
       * @param physical Moved physical coordinates.
       * @param locator Locator used to resolve the moved point.
       * @returns Located moved point, or an inverse-mapped point in the source volume cell or an adjacent volume cell when location fails.
       */
      static Geometry::Point makeMovedPoint(const Geometry::Point& source,
        const Math::SpatialVector<Real>& physical, const LocatorType& locator)
      {
        const auto& sourcePolytope = source.getPolytope();
        const auto& mesh = sourcePolytope.getMesh();
        const std::size_t dimension = mesh.getDimension();
        if (auto located = locator.locate(dimension, physical))
          return std::move(*located);

        Index cell = Index(-1);
        if (sourcePolytope.getDimension() == dimension)
          cell = sourcePolytope.getIndex();
        else
        {
          assert(sourcePolytope.getDimension() + 1 == dimension);
          const auto& adjacent = mesh.getConnectivity()
                                   .getIncidence(dimension - 1, dimension)
                                   .at(sourcePolytope.getIndex());
          if (adjacent.empty())
            throw std::runtime_error("Deformation map source face has no adjacent cell.");
          cell = *adjacent.begin();
        }

        Math::SpatialPoint reference;
        mesh.getPolytopeTransformation(dimension, cell).inverse(reference, physical);
        return Geometry::Point(*mesh.getPolytope(dimension, cell), reference, physical);
      }

      struct Key
      {
        const Geometry::MeshBase* mesh;
        std::size_t dimension;
        Index index;
        const QF::QuadratureFormulaBase* quadrature;
        size_t quadratureIdentity;
        std::size_t qp;
        Optional<size_t> transformationIdentity;
        /**
           * @brief Compares the operands for equality.
           * @param other Object compared with this object.
           * @returns True when the operands compare equal.
           */

        bool operator==(const Key& other) const
        {
          return mesh == other.mesh && dimension == other.dimension &&
            index == other.index && quadrature == other.quadrature &&
            quadratureIdentity == other.quadratureIdentity && qp == other.qp &&
            transformationIdentity == other.transformationIdentity;
        }
      };

      std::reference_wrapper<const Displacement> m_u;
      std::reference_wrapper<const LocatorType> m_locator;
      mutable Optional<Key> m_key;
      mutable Optional<Geometry::Point> m_moved;
  };
}

#endif
