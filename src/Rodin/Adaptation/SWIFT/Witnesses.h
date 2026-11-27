/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_WITNESSES_H
#define RODIN_ADAPTATION_SWIFT_WITNESSES_H

#include "Rodin/Geometry/PolytopeTransformation.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Parameters.h"
#include "QualityCovering.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Geometric witnesses on a background cell or interface face.
   *
   * The covering need not include the vertices. The mapped points are borrowed
   * from the mesh until its geometry changes. They are not an integration rule
   * of prescribed polynomial exactness.
   */
  class Witnesses
  {
    public:
      /**
       * @brief Selects a covering from the element order and geometric transformation.
       * @param entity Background cell or face.
       * @param order Order used by the sampling policy.
       * @param sampling Witness subdivision policy.
       */
      Witnesses(const Geometry::Polytope& entity, size_t order,
        const Parameters::Sampling& sampling)
        : m_points(entity.getQuadrature(QualityCovering::get(entity.getGeometry(),
            sampling.getSubdivision(order,
              entity.getTransformation().getDegree().value_or(2),
              Geometry::Polytope::Traits(entity.getGeometry()).getVertexCount() ==
                entity.getDimension() + 1))))
      {}

      /// @returns The mesh-owned mapped covering, shared until a geometry flush.
      const Geometry::PolytopeQuadrature& getQuadrature() const
      {
        return m_points;
      }

      /// @returns Number of witnesses.
      size_t getSize() const
      {
        return m_points.getSize();
      }

      /**
       * @param q Witness index in covering order.
       * @returns A point retaining its mapped quadrature and logical reference index.
       */
      Variational::IntegrationPoint getPoint(size_t q) const
      {
        return {m_points, q};
      }

      /**
       * @param q Witness index in covering order.
       * @returns The equal reference mass mapped to the background entity at this point.
       */
      Real getWeight(size_t q) const
      {
        return m_points.getQuadratureFormula().getWeight(0) *
          m_points.getPoint(q).getDistortion();
      }

      /**
       * @brief Visits points without evaluating weights.
       * @param evaluate Callable accepting a mapped integration point.
       */
      template <class Evaluate>
      void forEach(Evaluate&& evaluate) const
      {
        for (size_t q = 0; q < getSize(); ++q)
          evaluate(getPoint(q));
      }

    private:
      const Geometry::PolytopeQuadrature& m_points;
  };
}

#endif
