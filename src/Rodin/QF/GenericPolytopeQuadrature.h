/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_H
#define RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_H

/**
 * @file
 * @brief Defines the GenericPolytopeQuadrature dispatcher class.
 */

/**
 * @ingroup RodinDirectives
 * @brief Default order for GenericPolytopeQuadrature.
 *
 * This macro defines the default polynomial degree of exactness when
 * constructing a GenericPolytopeQuadrature without specifying an order.
 */
#define RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_DEFAULT_ORDER 1


#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Generic quadrature formula dispatcher for all polytope types.
   *
   * GenericPolytopeQuadrature provides a convenient interface for creating
   * quadrature formulas without needing to know which specific quadrature
   * rule to use for each geometry type. It automatically selects an
   * appropriate quadrature formula based on:
   * - The polytope geometry (triangle, tetrahedron, etc.)
   * - The desired polynomial degree of exactness
   *
   * ## Automatic Rule Selection
   *
   * The class internally selects the most appropriate quadrature rule:
   * - For simplices (triangles, tetrahedra): Grundmann-Möller rules
   * - For tensor-product elements (quads, hexahedra, wedges): Gauss-Legendre rules
   * - For segments: Gauss-Legendre rules
   * - For points: Single point with unit weight
   *
   * ## Usage
   *
   * @code{.cpp}
   * // Create a quadrature rule exact for degree 3 polynomials on a triangle
   * GenericPolytopeQuadrature qf(3, Geometry::Polytope::Type::Triangle);
   *
   * // Use default order (degree 1)
   * GenericPolytopeQuadrature qf_default(Geometry::Polytope::Type::Tetrahedron);
   * @endcode
   *
   * @see QuadratureFormulaBase
   * @see GaussLegendre
   * @see GrundmannMoller
   */
  class GenericPolytopeQuadrature : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      static const QuadratureFormulaBase& get(size_t order, Geometry::Polytope::Type g);

      static std::unique_ptr<const QuadratureFormulaBase> build(size_t order, Geometry::Polytope::Type g);

      /**
       * @brief Constructs quadrature with default order.
       * @param g Polytope geometry type
       *
       * Creates a quadrature formula with the default polynomial degree
       * of exactness (defined by 
       * RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_DEFAULT_ORDER).
       */
      GenericPolytopeQuadrature(Geometry::Polytope::Type g)
        : GenericPolytopeQuadrature(RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_DEFAULT_ORDER, g)
      {}

      /**
       * @brief Constructs quadrature with specified order.
       * @param order Desired polynomial degree of exactness
       * @param g Polytope geometry type
       *
       * Automatically selects and constructs an appropriate quadrature
       * formula that is exact for all polynomials of degree up to @p order.
       */
      GenericPolytopeQuadrature(size_t order, Geometry::Polytope::Type g);

      /**
       * @brief Copy constructor.
       * @param other GenericPolytopeQuadrature to copy from
       */
      GenericPolytopeQuadrature(const GenericPolytopeQuadrature& other)
        : Parent(other),
          m_qf(other.m_qf->copy()),
          m_order(other.m_order)
      {}

      /**
       * @brief Move constructor.
       * @param other GenericPolytopeQuadrature to move from
       */
      GenericPolytopeQuadrature(GenericPolytopeQuadrature&& other)
        : Parent(std::move(other)),
          m_qf(std::move(other.m_qf)),
          m_order(std::move(other.m_order))
      {}

      /**
       * @brief Gets the number of quadrature points.
       * @return Total number of quadrature points in the underlying formula
       */
      inline
      size_t getSize() const override
      {
        return m_qf->getSize();
      }

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight @f$ w_i @f$ from the underlying formula
       */
      inline
      Real getWeight(size_t i) const override
      {
        return m_qf->getWeight(i);
      }

      /**
       * @brief Gets the coordinates of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Reference to the point coordinates from the underlying formula
       */
      inline
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        return m_qf->getPoint(i);
      }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new GenericPolytopeQuadrature instance
       */
      inline
      GenericPolytopeQuadrature* copy() const noexcept override
      {
        return new GenericPolytopeQuadrature(*this);
      }

    private:
      std::unique_ptr<const QuadratureFormulaBase> m_qf; ///< Underlying quadrature formula
      const size_t m_order; ///< Polynomial degree of exactness
  };
}

#endif
