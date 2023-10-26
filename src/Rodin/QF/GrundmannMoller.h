/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_GRUNDMANNMOLLER_H
#define RODIN_VARIATIONAL_QF_GRUNDMANNMOLLER_H

/**
 * @ingroup RodinDirectives
 * @brief Maximum value permitted for the @f$ s \geq 0 @f$ parameter in
 * Grundmann-Moller quadrature.
 * @see Rodin::QF::GrundmannMoller
 * @see RODIN_QF_GRUNDMANNMOLLER_MAX_ORDER
 */
#define RODIN_QF_GRUNDMANNMOLLER_MAX_S 16

/**
 * @ingroup RodinDirectives
 * @brief Maximum degree permitted for Grundmann-Moller quadrature.
 *
 * The degree is given by @f$ d = 2s + 1 @f$ where @f$ s \geq 0 @f$.
 * @see Rodin::QF::GrundmannMoller
 * @see RODIN_QF_GRUNDMANNMOLLER_MAX_S
 */
#define RODIN_QF_GRUNDMANNMOLLER_MAX_ORDER 2 * RODIN_QF_GRUNDMANNMOLLER_MAX_S + 1

#include <boost/multi_array.hpp>

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Grundmann-Möller Simplex Quadrature
   *
   * Implements the Grundmann-Möller integration rules as described in
   * @cite grundmann1978invariant. These quadrature formulae are defined only
   * on simplices and have an odd degree @f$ d = 2s + 1 @f$ parametrized by @f$
   * s \geq 0 @f$.
   */
  class GrundmannMoller : public QuadratureFormulaBase
  {
    public:
      /// Parent class
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs the quadrature formula.
       * @param[in] s Parameter which will indicate the degree.
       * @param[in] geom Simplex where the quadrature is defined.
       */
      GrundmannMoller(size_t s, Geometry::Polytope::Type geom);

      constexpr
      GrundmannMoller(const GrundmannMoller& other)
        : Parent(other),
          m_s(other.m_s),
          m_n(other.m_n),
          m_order(other.m_order)
      {}

      constexpr
      GrundmannMoller(GrundmannMoller&& other)
        : Parent(std::move(other)),
          m_s(std::move(other.m_s)),
          m_n(std::move(other.m_n)),
          m_order(std::move(other.m_order))
      {}

      inline
      size_t getOrder() const
      {
        return m_order;
      }

      size_t getSize() const override;

      Scalar getWeight(size_t i) const override;

      const Math::SpatialVector& getPoint(size_t i) const override;

      inline
      GrundmannMoller* copy() const noexcept override
      {
        return new GrundmannMoller(*this);
      }

      static boost::multi_array<size_t, 2> initSizes();
      static boost::multi_array<Math::Vector, 2> initWeights();
      static boost::multi_array<std::vector<Math::SpatialVector>, 2> initPoints();
    private:

      static boost::multi_array<size_t, 2> s_sizes;
      static boost::multi_array<Math::Vector, 2> s_weights;
      static boost::multi_array<std::vector<Math::SpatialVector>, 2> s_points;

      const size_t m_s;
      const size_t m_n;
      const size_t m_order;
  };
}

#endif
