/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_ADAPTATION_DEFORMATIONMAP_H
#define RODIN_ADAPTATION_DEFORMATIONMAP_H

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/IntegrationPoint.h"

namespace Rodin::Adaptation
{
  /**
   * @brief The deformation @f$\Phi(x)=x+u(x)@f$ induced by a displacement.
   *
   * ## Architecture
   *
   * A point binding locates the moved point once and owns that result for all
   * subsequent basis requests in the binding. Cell binding or invalidation
   * releases it. Direct point evaluation always locates afresh. Neither the
   * integration point nor a transformation-state token is retained. Geometry,
   * displacement and locator data must remain unchanged during a binding.
   * Copies retain the displacement and locator references but start unbound.
   * Evaluation objects are local to a pass and must not be shared concurrently.
   *
   * @tparam Displacement Displacement grid function type.
   * @tparam LocatorType Point locator used in the displaced configuration.
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
       * @brief Copies dependencies into an unbound evaluation object.
       * @param other Evaluation whose displacement and locator are retained.
       */
      DeformationMap(const DeformationMap& other)
        : m_u(other.m_u),
          m_locator(other.m_locator)
      {}

      /**
       * @brief Starts a fresh cell binding, releasing the previous moved point.
       * @param polytope Source entity for the next evaluation pass.
       * @returns This unbound point evaluator.
       */
      DeformationMap& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope)
      {
        invalidate();
        return *this;
      }

      /**
       * @brief Locates and owns the moved point for a fresh point binding.
       * @param ip Source integration point; no reference to it is retained.
       * @returns This evaluator with an owned moved-point result.
       */
      DeformationMap& setIntegrationPoint(const Variational::IntegrationPoint& ip)
      {
        getMovedPoint(ip);
        return *this;
      }

      /**
       * @brief Gets the moved point prepared by the current point binding.
       * @returns The owned point, valid until rebinding or invalidation.
       * @throws Alert::Exception If no point is bound.
       */
      const Geometry::Point& getMovedPoint() const
      {
        if (!m_moved)
          Alert::MemberFunctionException(*this, __func__)
            << "No moved point is bound." << Alert::Raise;
        return *m_moved;
      }

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
       * @brief Evaluates and owns a fresh moved point outside a binding.
       * @param ip Source integration point.
       * @returns Located moved point, replaced by the next evaluation.
       */
      const Geometry::Point& getMovedPoint(const Variational::IntegrationPoint& ip)
      {
        invalidate();
        const auto& point = ip.getPoint();
        m_moved.emplace(makeMovedPoint(
          point, point.getPhysicalCoordinates() + getDisplacementValue(point, ip)));
        return *m_moved;
      }

      /// @brief Drops the cached point, forcing the next request to relocate.
      void invalidate()
      {
        m_moved.reset();
      }

    private:
      /**
       * @brief Resolves a point at moved physical coordinates.
       * @param source Source point retaining the reference cell and coordinates.
       * @param physical Moved physical coordinates.
       * @returns Located moved point, or an inverse-mapped point in the source volume cell or an adjacent volume cell when location fails.
       */
      Geometry::Point makeMovedPoint(
        const Geometry::Point& source, const Math::SpatialVector<Real>& physical) const
      {
        const auto& sourcePolytope = source.getPolytope();
        const auto& mesh = sourcePolytope.getMesh();
        const std::size_t dimension = mesh.getDimension();
        if (auto located = m_locator.get().locate(dimension, physical))
          return std::move(*located);

        Index cell;
        if (sourcePolytope.getDimension() == dimension)
          cell = sourcePolytope.getIndex();
        else
        {
          if (sourcePolytope.getDimension() + 1 != dimension)
            Alert::MemberFunctionException(*this, __func__)
              << "Deformation fallback requires a source cell or facet." << Alert::Raise;
          const auto& adjacent = mesh.getConnectivity()
                                   .getIncidence(dimension - 1, dimension)
                                   .at(sourcePolytope.getIndex());
          if (adjacent.empty())
            Alert::MemberFunctionException(*this, __func__)
              << "Deformation map source face has no adjacent cell." << Alert::Raise;
          cell = *adjacent.begin();
        }

        Math::SpatialPoint reference;
        mesh.getPolytopeTransformation(dimension, cell).inverse(reference, physical);
        return Geometry::Point(*mesh.getPolytope(dimension, cell), reference, physical);
      }

      std::reference_wrapper<const Displacement> m_u;
      std::reference_wrapper<const LocatorType> m_locator;
      Optional<Geometry::Point> m_moved;
  };
}

#endif
