/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_TENSORPRODUCT_H
#define RODIN_QF_TENSORPRODUCT_H

/**
 * @file
 * @brief Owns Cartesian products of reference quadrature formulas.
 */

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Cartesian product of two quadrature formulae.
   *
   * The coordinates of the right factor are appended to those of the left
   * factor and the corresponding weights are multiplied. This covers both
   * tensor-product cells and products of unlike factors, such as a triangular
   * rule times a segment rule on a wedge.
   */
  class TensorProduct final : public QuadratureFormulaBase
  {
    public:
      /**
       * @brief Constructs the Cartesian product of two quadrature formulae.
       * @param g Reference geometry of the Cartesian product.
       * @param left Quadrature formula for the first tensor-product factor.
       * @param right Quadrature formula for the second tensor-product factor.
       */
      TensorProduct(Geometry::Polytope::Type g, const QuadratureFormulaBase& left,
        const QuadratureFormulaBase& right)
        : m_geometry(g)
      {
        append(left, right);
      }

      /**
       * @brief Constructs the Cartesian product of three quadrature formulae.
       * @param g Reference geometry of the Cartesian product.
       * @param first Quadrature formula for the first tensor-product factor.
       * @param second Quadrature formula for the second tensor-product factor.
       * @param third Quadrature formula for the third tensor-product factor.
       */
      TensorProduct(Geometry::Polytope::Type g, const QuadratureFormulaBase& first,
        const QuadratureFormulaBase& second, const QuadratureFormulaBase& third)
        : m_geometry(g)
      {
        TensorProduct pair(g, first, second);
        append(pair, third);
      }

      /**
       * @brief Copies a tensor-product quadrature formula.
       * @param other Object to copy from.
       */
      TensorProduct(const TensorProduct& other) = default;

      size_t getSize() const override
      {
        return m_points.size();
      }

      Real getWeight(size_t i) const override
      {
        assert(i < m_weights.size());
        return m_weights[i];
      }

      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i < m_points.size());
        return m_points[i];
      }

      /**
       * @brief Returns the reference element associated with the product rule.
       * @returns The reference element associated with the product rule.
       */
      Geometry::Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      /**
       * @brief Copies the concrete formula without regenerating its samples.
       * @returns A new formula owned by the caller.
       */
      TensorProduct* copy() const noexcept
      {
        return new TensorProduct(*this);
      }

    private:
      /**
       * @brief Appends product quadrature points and weights.
       *
       * Factor sizes and left samples are invariant over their inner loops.
       * References are borrowed only during construction; the product owns all
       * coordinates and weights. Concatenation order and floating-point weight
       * multiplication are unchanged, including for aliased factors.
       * Each point is zero-constructed in its final owning slot before filling
       * active coordinates, avoiding a temporary fixed-capacity vector move.
       *
       * @param left Left quadrature formula.
       * @param right Right quadrature formula.
       */
      void append(const QuadratureFormulaBase& left, const QuadratureFormulaBase& right)
      {
        const size_t leftCount = left.getSize(), rightCount = right.getSize();
        if (leftCount == 0 || rightCount == 0)
          return;
        m_points.reserve(leftCount * rightCount);
        m_weights.reserve(leftCount * rightCount);
        for (size_t i = 0; i < leftCount; ++i)
        {
          const auto& leftPoint = left.getPoint(i);
          const auto leftDimension = leftPoint.size();
          const Real leftWeight = left.getWeight(i);
          for (size_t j = 0; j < rightCount; ++j)
          {
            const auto& rightPoint = right.getPoint(j);
            const auto rightDimension = rightPoint.size();
            auto& p = m_points.emplace_back(leftDimension + rightDimension);
            for (Eigen::Index k = 0; k < leftDimension; ++k)
              p[k] = leftPoint[k];
            for (Eigen::Index k = 0; k < rightDimension; ++k)
              p[leftDimension + k] = rightPoint[k];
            m_weights.push_back(leftWeight * right.getWeight(j));
        }
      }
      }

      Geometry::Polytope::Type m_geometry;
      std::vector<Math::SpatialVector<Real>> m_points;
      std::vector<Real> m_weights;
  };
}

#endif
