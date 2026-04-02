/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QUADRATUREFORMULA_H
#define RODIN_VARIATIONAL_QUADRATUREFORMULA_H

/**
 * @file
 * @brief Defines the QuadratureFormulaBase abstract base class.
 */

#include <vector>
#include <utility>

#include "Rodin/Math.h"
#include "Rodin/Copyable.h"
#include "Rodin/Identifiable.h"
#include "Rodin/Geometry.h"


namespace Rodin::QF
{
  /**
   * @defgroup RodinQuadrature Quadrature formulae
   * @brief Quadrature formulae utilized by Rodin
   * @see QuadratureFormulaBase
   */

  /**
   * @brief Abstract base class for quadrature formulas.
   * @ingroup RodinQuadrature
   *
   * QuadratureFormulaBase defines the interface for all quadrature formulas
   * in Rodin. A quadrature formula consists of a set of points (nodes) and
   * associated weights that approximate integrals over a reference polytope.
   *
   * Each quadrature formula is associated with a specific polytope geometry
   * and provides methods to query:
   * - The number of quadrature points
   * - The coordinates of each point in the reference element
   * - The weight associated with each point
   *
   * Concrete implementations include Gauss-Legendre, Gauss-Lobatto,
   * Grundmann-Möller, and other specialized rules for various polytope types.
   */
  class QuadratureFormulaBase : public Copyable
  {
    public:
      constexpr QuadratureFormulaBase() = default;

      /**
       * @brief Copy constructor.
       * @param other Another quadrature formula to copy from
       */
      constexpr
      QuadratureFormulaBase(const QuadratureFormulaBase& other) = default;

      /**
       * @brief Virtual destructor.
       */
      virtual ~QuadratureFormulaBase() = default;

      /**
       * @brief Gets the number of quadrature points.
       * @return The total number of quadrature points in this formula
       */
      virtual size_t getSize() const = 0;

      /**
       * @brief Gets the weight associated with the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight @f$ w_i @f$ for the i-th point
       */
      virtual Real getWeight(size_t i) const = 0;

      /**
       * @brief Gets the coordinates of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Reference to the coordinates @f$ x_i @f$ in reference space
       *
       * @note The returned reference must remain valid for the lifetime
       * of the quadrature formula object.
       */
      virtual const Math::SpatialVector<Real>& getPoint(size_t i) const = 0;

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new dynamically allocated copy
       *
       * @note The caller is responsible for managing the lifetime of the
       * returned pointer.
       */
      virtual QuadratureFormulaBase* copy() const noexcept override = 0;
  };
}

#endif
