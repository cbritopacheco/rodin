/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_QFVertexP1_H
#define RODIN_VARIATIONAL_QF_QFVertexP1_H

#include "Rodin/Geometry/GeometryIndexed.h"
#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Vertex-based quadrature formula with degree of exactness 1.
   *
   * QFVertexP1 implements a quadrature rule that places quadrature points
   * at the vertices of the reference polytope. This is useful for certain
   * finite element applications, particularly when working with vertex-based
   * function evaluations or when a simple low-order rule is sufficient.
   *
   * ## Characteristics
   *
   * - Quadrature points are located at polytope vertices
   * - Exact for linear (degree 1) polynomials
   * - Number of points equals number of vertices for the geometry
   * - Weights are chosen to ensure exactness for linear functions
   *
   * ## Supported Geometries
   *
   * Works with all standard polytope geometries (segments, triangles,
   * quadrilaterals, tetrahedra, etc.) that have well-defined vertices.
   *
   * @note This quadrature rule is primarily used for special applications
   * and may not provide optimal accuracy for general integration tasks.
   * Consider GaussLegendre or other higher-order rules for better accuracy.
   */
  class QFVertexP1 final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs a QFVertexP1 quadrature formula for the given geometry.
       * @param g Geometry type (e.g., triangle, quadrilateral, etc.)
       *
       * Initializes the vertex-based quadrature rule with points at the
       * vertices of the reference element for geometry @p g.
       */
      constexpr QFVertexP1(Geometry::Polytope::Type g)
        : Parent(g)
      {}

      /**
       * @brief Gets the number of quadrature points (vertices).
       * @return Number of vertices in the reference polytope
       */
      size_t getSize() const override
      {
        return s_size[getGeometry()];
      }

      /**
       * @brief Gets the coordinates of the i-th quadrature point (vertex).
       * @param i Index of the quadrature point (vertex index)
       * @return Reference to the vertex coordinates in reference space
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        return s_points[getGeometry()][i];
      }

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight @f$ w_i @f$ associated with the i-th vertex
       */
      Real getWeight(size_t i) const override
      {
        return s_weights[getGeometry()][i];
      }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new QFVertexP1 instance (ownership transferred to caller)
       */
      QFVertexP1* copy() const noexcept override
      {
        return new QFVertexP1(*this);
      }

    private:
      /// Static data: number of vertices for each geometry type
      static const Geometry::GeometryIndexed<size_t> s_size;
      
      /// Static data: vertex coordinates for each geometry type
      static const Geometry::GeometryIndexed<std::vector<Math::SpatialVector<Real>>> s_points;
      
      /// Static data: vertex weights for each geometry type
      static const Geometry::GeometryIndexed<Math::Vector<Real>> s_weights;
  };
}

#endif
