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
   * A point binding locates the deformed point once and owns that result for all
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
       * @brief Starts a fresh cell binding, releasing the previous deformed point.
       * @param polytope Source entity for the next evaluation pass.
       * @returns This unbound point evaluator.
       */
      DeformationMap& setPolytope([[maybe_unused]] const Geometry::Polytope& polytope)
      {
        invalidate();
        return *this;
      }

      /**
       * @brief Locates and owns the deformed point for a fresh point binding.
       *
       * The locator resolves the displaced physical coordinates. If location
       * fails, the point is inverse-mapped in the source cell or an adjacent
       * volume cell for a source facet. Every bind evaluates afresh, including
       * repeated binds at an unchanged integration point.
       * @param ip Source integration point; no reference to it is retained.
       * @returns This evaluator with an owned deformed-point result.
       */
      DeformationMap& setIntegrationPoint(const Variational::IntegrationPoint& ip)
      {
        invalidate();
        const auto& point = ip.getPoint();
        const Math::SpatialPoint physical =
          point.getPhysicalCoordinates() + getDisplacementValue(point, ip);
        const auto& sourcePolytope = point.getPolytope();
        const auto& mesh = sourcePolytope.getMesh();
        const std::size_t dimension = mesh.getDimension();
        if (auto located = m_locator.get().locate(dimension, physical))
        {
          m_moved.emplace(std::move(*located));
          return *this;
        }

        Index cell;
        if (sourcePolytope.getDimension() == dimension)
          cell = sourcePolytope.getIndex();
        else
        {
          if (sourcePolytope.getDimension() + 1 != dimension)
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Deformation fallback requires a source cell or facet." << Alert::Raise;
          }
          const auto& adjacent = mesh.getConnectivity()
                                   .getIncidence(dimension - 1, dimension)
                                   .at(sourcePolytope.getIndex());
          if (adjacent.empty())
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Deformation map source face has no adjacent cell." << Alert::Raise;
          }
          cell = *adjacent.begin();
        }

        Math::SpatialPoint reference;
        mesh.getPolytopeTransformation(dimension, cell).inverse(reference, physical);
        m_moved.emplace(*mesh.getPolytope(dimension, cell), reference, physical);
        return *this;
      }

      /**
       * @brief Gets the deformed point prepared by the current point binding.
       * @returns The owned point, valid until rebinding or invalidation.
       * @throws Alert::Exception If no point is bound.
       */
      const Geometry::Point& getPoint() const
      {
        if (!m_moved)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No deformed point is bound." << Alert::Raise;
        }
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
        Math::SpatialVector<Real> displacement(
          pt.getPolytope().getMesh().getSpaceDimension());
        if (static_cast<size_t>(uq.size()) != static_cast<size_t>(displacement.size()))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Displacement dimension must match the mesh ambient dimension."
            << Alert::Raise;
        }
        displacement.setZero();
        for (std::size_t r = 0; r < static_cast<std::size_t>(displacement.size()); ++r)
          displacement(static_cast<Eigen::Index>(r)) = uq(r);
        return displacement;
      }

      /// @brief Drops the cached point, forcing the next request to relocate.
      void invalidate()
      {
        m_moved.reset();
      }

    private:
      std::reference_wrapper<const Displacement> m_u;
      std::reference_wrapper<const LocatorType> m_locator;
      Optional<Geometry::Point> m_moved;
  };
}

#endif
