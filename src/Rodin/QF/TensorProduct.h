/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_TENSORPRODUCT_H
#define RODIN_QF_TENSORPRODUCT_H

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
      /// @brief Constructs the Cartesian product of two quadrature formulae.
      TensorProduct(Geometry::Polytope::Type g, const QuadratureFormulaBase& left,
        const QuadratureFormulaBase& right)
        : m_geometry(g)
      {
        append(left, right);
      }

      /// @brief Constructs the Cartesian product of three quadrature formulae.
      TensorProduct(Geometry::Polytope::Type g, const QuadratureFormulaBase& first,
        const QuadratureFormulaBase& second, const QuadratureFormulaBase& third)
        : m_geometry(g)
      {
        TensorProduct pair(g, first, second);
        append(pair, third);
      }

      /// @brief Copies a tensor-product quadrature formula.
      TensorProduct(const TensorProduct&) = default;

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

      /// @brief Returns the reference element associated with the product rule.
      Geometry::Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      TensorProduct* copy() const noexcept override
      {
        return new TensorProduct(*this);
      }

    private:
      void append(const QuadratureFormulaBase& left, const QuadratureFormulaBase& right)
      {
        m_points.reserve(left.getSize() * right.getSize());
        m_weights.reserve(left.getSize() * right.getSize());
        for (size_t i = 0; i < left.getSize(); ++i)
        {
          for (size_t j = 0; j < right.getSize(); ++j)
          {
            Math::SpatialVector<Real> p(
              left.getPoint(i).size() + right.getPoint(j).size());
            for (Eigen::Index k = 0; k < left.getPoint(i).size(); ++k)
              p[k] = left.getPoint(i)[k];
            for (Eigen::Index k = 0; k < right.getPoint(j).size(); ++k)
              p[left.getPoint(i).size() + k] = right.getPoint(j)[k];
            m_points.push_back(std::move(p));
            m_weights.push_back(left.getWeight(i) * right.getWeight(j));
          }
        }
      }

      Geometry::Polytope::Type m_geometry;
      std::vector<Math::SpatialVector<Real>> m_points;
      std::vector<Real> m_weights;
  };
}

#endif
