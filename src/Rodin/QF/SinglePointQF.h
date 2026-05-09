/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_SINGLEPOINTQF_H
#define RODIN_QF_SINGLEPOINTQF_H

/**
 * @file
 * @brief Defines a single-point quadrature formula with a user-supplied
 *        reference coordinate.
 */

#include <limits>

#include "QuadratureFormula.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Single-point quadrature formula at an arbitrary reference point.
   *
   * Holds one reference-element point (and unit weight) so that consumers of
   * @ref QuadratureFormulaBase that fetch the formula's point — for instance
   * the @c ShapeFunction caches built around @c qf.getPoint(qp) — sample at
   * the supplied location instead of a fixed canonical point. Used to
   * evaluate a @ref Variational::ShapeFunctionBase at an arbitrary geometric
   * point inside the identification Dirichlet BC assembler.
   *
   * Stores the point by value to keep the lifetime self-contained.
   */
  class SinglePointQF final : public QuadratureFormulaBase
  {
    public:
      SinglePointQF(const Math::SpatialVector<Real>& point)
        : m_point(point)
      {}

      SinglePointQF(Math::SpatialVector<Real>&& point)
        : m_point(std::move(point))
      {}

      SinglePointQF(const SinglePointQF& other) = default;

      // Reports a large nominal size so consumers that key caches on
      // (qf, qp) pairs — e.g. shape-function basis caches — can pass a
      // monotonically increasing index to invalidate stale cached values
      // when the stored point is changed across calls. Regardless of the
      // index, getPoint and getWeight return the single stored entry.
      size_t getSize() const override
      {
        return std::numeric_limits<size_t>::max();
      }

      Real getWeight(size_t i) const override
      {
        (void) i;
        return Real(1);
      }

      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        (void) i;
        return m_point;
      }

      SinglePointQF* copy() const noexcept override
      {
        return new SinglePointQF(*this);
      }

    private:
      Math::SpatialVector<Real> m_point;
  };
}

#endif
