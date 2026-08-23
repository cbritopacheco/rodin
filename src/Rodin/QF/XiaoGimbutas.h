/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_XIAOGIMBUTAS_H
#define RODIN_QF_XIAOGIMBUTAS_H

/**
 * @file
 * @brief Defines the XiaoGimbutas simplex quadrature formula.
 */

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Near-minimal positive interior quadrature on simplices.
   *
   * Implements rules of the family described by Xiao and Gimbutas
   * @cite xiao2010numerical: begin from a rule that is already exact, remove
   * the node contributing least to the moments, and solve the moment
   * equations again from the reduced rule, repeating until nothing more can
   * go. The result has positive weights and interior nodes, and uses close to
   * the fewest points the moment count allows.
   *
   * The coefficients shipped here were computed by Rodin, by
   * NodeElimination::reduce, and not transcribed from any published table.
   * They are therefore not guaranteed to coincide with the rules of
   * @cite xiao2010numerical even where the point counts agree: the moment
   * equations have many solutions and the search finds one of them. What is
   * guaranteed, and tested, is the property that matters --- exactness to the
   * stated degree, positive weights, and nodes inside the element.
   *
   * Node elimination does not preserve symmetry, so these rules are not
   * symmetric, unlike those of Witherden and Vincent
   * @cite witherden2015identification. At tetrahedron degree 4 that lets the
   * count fall below the published symmetric minimum.
   *
   * Grundmann-Möller remains the fallback beyond the tabulated degrees and for
   * simplices of dimension above three, where it is the only rule available.
   *
   * @see NodeElimination
   * @see GrundmannMoller
   */
  class XiaoGimbutas final : public QuadratureFormulaBase
  {
    public:
      /// @brief Highest tabulated degree for @p g, or zero if unsupported.
      static size_t getMaxDegree(Geometry::Polytope::Type g);

      /// @brief Whether a rule of strength @p degree is tabulated for @p g.
      static bool isAvailable(size_t degree, Geometry::Polytope::Type g)
      {
        return degree >= 1 && degree <= getMaxDegree(g);
      }

      /**
       * @brief Constructs the rule of strength @p degree on @p g.
       * @pre isAvailable(degree, g)
       */
      XiaoGimbutas(size_t degree, Geometry::Polytope::Type g);

      XiaoGimbutas(const XiaoGimbutas&) = default;

      size_t getSize() const override
      {
        return m_count;
      }

      Real getWeight(size_t i) const override
      {
        assert(i < m_count);
        return m_data[i * (m_dimension + 1) + m_dimension];
      }

      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i < m_count);
        return m_points[i];
      }

      /// @brief The element this rule is defined on.
      Geometry::Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      XiaoGimbutas* copy() const noexcept override
      {
        return new XiaoGimbutas(*this);
      }

    private:
      Geometry::Polytope::Type m_geometry;
      size_t m_dimension;
      size_t m_count;
      const Real* m_data;
      std::vector<Math::SpatialVector<Real>> m_points;
  };
}

#endif
