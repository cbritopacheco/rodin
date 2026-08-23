/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_WITHERDENVINCENT_H
#define RODIN_QF_WITHERDENVINCENT_H

/**
 * @file
 * @brief Defines the WitherdenVincent quadrature formula.
 */

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Fully symmetric positive interior quadrature.
   *
   * Implements rules of the family described by Witherden and Vincent
   * @cite witherden2015identification. Where NodeElimination discovers a rule
   * by removing nodes from a larger one, and in doing so destroys any
   * symmetry it had, these rules are symmetric by construction: the nodes are
   * organised into orbits of the element's symmetry group and the moment
   * equations are solved for the orbit parameters directly.
   *
   * Symmetry costs points and buys structure. A symmetric rule commutes with
   * vertex renumbering, so neighbouring elements of a mesh integrate a shared
   * face identically however their vertices happen to be ordered, which an
   * asymmetric rule does not guarantee.
   *
   * It also reaches counts an eliminating search cannot. At tetrahedron degree
   * 5 the minimum attains the counting bound exactly --- 56 moments against 14
   * nodes carrying 4 unknowns each --- which leaves the moment system square
   * and its solutions isolated. NodeElimination relies throughout on the null
   * space of an underdetermined system and stops at 15 points there; imposing
   * symmetry reduces the unknowns to a handful of orbit parameters and finds
   * the 14-point rule directly.
   *
   * The coefficients shipped here were computed by Rodin, by
   * SymmetricRuleSolver::search, and not transcribed from any published table.
   * They are therefore not guaranteed to coincide with the rules of
   * @cite witherden2015identification even where the point counts agree. What
   * is guaranteed, and tested, is exactness to the stated degree, positive
   * weights, and nodes inside the element.
   *
   * @see SymmetricRuleSolver
   * @see XiaoGimbutas
   * @see GrundmannMoller
   */
  class WitherdenVincent final : public QuadratureFormulaBase
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
      WitherdenVincent(size_t degree, Geometry::Polytope::Type g);

      WitherdenVincent(const WitherdenVincent&) = default;

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

      WitherdenVincent* copy() const noexcept override
      {
        return new WitherdenVincent(*this);
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
