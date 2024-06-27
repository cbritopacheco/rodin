/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_H
#define RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_H

/**
 * @ingroup RodinDirectives
 * @brief Default order for GenericPolytopeQuadrature.
 */
#define RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_DEFAULT_ORDER 1

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Quadrature on a polytope with any of @ref
   * Geometry::Polytope::Geometry "Rodin's supported geometries".
   * @see @ref Geometry::Polytope::Geometry "Polytope::Geometry"
   */
  class GenericPolytopeQuadrature : public QuadratureFormulaBase
  {
    public:
      /// Parent class
      using Parent = QuadratureFormulaBase;

      GenericPolytopeQuadrature(Geometry::Polytope::Type g)
        : GenericPolytopeQuadrature(RODIN_VARIATIONAL_QF_GENERIC_POLYTOPE_QUADRATURE_DEFAULT_ORDER, g)
      {}

      GenericPolytopeQuadrature(size_t order, Geometry::Polytope::Type g);

      GenericPolytopeQuadrature(const GenericPolytopeQuadrature& other)
        : Parent(other),
          m_qf(other.m_qf->copy()),
          m_order(other.m_order)
      {}

      GenericPolytopeQuadrature(GenericPolytopeQuadrature&& other)
        : Parent(std::move(other)),
          m_qf(std::move(other.m_qf)),
          m_order(std::move(other.m_order))
      {}

      inline
      size_t getSize() const override
      {
        return m_qf->getSize();
      }

      inline
      Scalar getWeight(size_t i) const override
      {
        return m_qf->getWeight(i);
      }

      inline
      const Math::SpatialVector<Scalar>& getPoint(size_t i) const override
      {
        return m_qf->getPoint(i);
      }

      inline
      GenericPolytopeQuadrature* copy() const noexcept override
      {
        return new GenericPolytopeQuadrature(*this);
      }

    private:
      std::unique_ptr<const QuadratureFormulaBase> m_qf;
      const size_t m_order;
  };
}

#endif
