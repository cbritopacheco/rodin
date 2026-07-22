/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_DEFORMATIONMAP_H
#define RODIN_ADAPTATION_DEFORMATIONMAP_H

#include "Rodin/Geometry/Point.h"
#include "Rodin/Types.h"

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
        : m_u(u), m_locator(locator)
      {}

      /// @brief The underlying displacement field.
      const Displacement& getDisplacement() const { return m_u.get(); }

      /**
       * @brief The displacement @f$u(x)@f$ at an integration point, as a
       * spatial vector.
       */
      Math::SpatialVector<Real> getDisplacementValue(
        const Geometry::Point& pt, const Variational::IntegrationPoint& ip) const
      {
        const auto uq = m_u.get().getValue(ip);
        Math::SpatialVector<Real> displacement(pt.getPolytope().getMesh().getDimension());
        displacement.setZero();
        for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
          displacement(static_cast<Eigen::Index>(r)) = uq(r);
        return displacement;
      }

      /**
       * @brief The moved point @f$\Phi(x)=x+u(x)@f$, located in the mesh.
       *
       * Cached: consecutive calls with the same polytope and quadrature index
       * return the previously located point without repeating the search.
       *
       * @param ip Integration point on the reference configuration; its
       * polytope and quadrature index form the cache key.
       */
      const Geometry::Point& getMovedPoint(
        const Variational::IntegrationPoint& ip) const
      {
        const auto& pt = ip.getPoint();
        const auto& polytope = pt.getPolytope();
        const Key key{polytope.getDimension(), polytope.getIndex(), ip.getIndex()};
        if (!m_moved || !m_key || !(*m_key == key))
        {
            // Geometry::Point is constructible but not assignable, so the slot
            // is re-emplaced rather than overwritten.
          m_moved.emplace(makeMovedPoint(pt,
            pt.getPhysicalCoordinates() + getDisplacementValue(pt, ip),
            m_locator.get()));
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
          const auto& adjacent = mesh.getConnectivity().getIncidence(
            dimension - 1, dimension).at(sourcePolytope.getIndex());
          if (adjacent.empty())
            throw std::runtime_error(
              "Deformation map source face has no adjacent cell.");
          cell = *adjacent.begin();
        }

        Math::SpatialPoint reference;
        mesh.getPolytopeTransformation(dimension, cell).inverse(reference, physical);
        return Geometry::Point(*mesh.getPolytope(dimension, cell), reference, physical);
      }

      struct Key
      {
        std::size_t dimension;
        Index index;
        std::size_t qp;

        bool operator==(const Key& other) const
        {
          return dimension == other.dimension && index == other.index &&
            qp == other.qp;
        }
      };

      std::reference_wrapper<const Displacement> m_u;
      std::reference_wrapper<const LocatorType> m_locator;
      mutable Optional<Key> m_key;
      mutable Optional<Geometry::Point> m_moved;
  };
}

#endif
