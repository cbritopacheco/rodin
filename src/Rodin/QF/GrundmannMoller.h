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
   * @brief Grundmann-Möller simplex quadrature formula.
   *
   * Implements the Grundmann-Möller integration rules as described in
   * @cite grundmann1978invariant. These quadrature formulae are specifically
   * designed for simplicial elements (segments, triangles, tetrahedra) and
   * provide exact integration for polynomials of odd degree @f$ d = 2s + 1 @f$
   * where @f$ s \geq 0 @f$ is a parameter.
   *
   * ## Mathematical Foundation
   *
   * For a simplex @f$ K @f$ of dimension @f$ n @f$, the Grundmann-Möller
   * rule approximates:
   * @f[
   *   \int_K f(x) \, dx \approx \sum_{i=1}^N w_i f(x_i)
   * @f]
   * where the points @f$ x_i @f$ and weights @f$ w_i @f$ are determined by
   * the parameter @f$ s @f$ and the simplex dimension @f$ n @f$.
   *
   * ## Degree of Exactness
   *
   * The quadrature rule is exact for all polynomials of degree up to
   * @f$ d = 2s + 1 @f$. For example:
   * - @f$ s = 0 @f$: degree 1 (linear exactness)
   * - @f$ s = 1 @f$: degree 3 (cubic exactness)
   * - @f$ s = 2 @f$: degree 5 (quintic exactness)
   *
   * ## Supported Geometries
   *
   * - Segment (1D simplex)
   * - Triangle (2D simplex)
   * - Tetrahedron (3D simplex)
   *
   * @note This class only supports simplicial geometries. For tensor-product
   * elements (quadrilaterals, hexahedra), use GaussLegendre or other
   * appropriate rules.
   *
   * @note The parameter @f$ s @f$ is limited by RODIN_QF_GRUNDMANNMOLLER_MAX_S.
   */
  class GrundmannMoller : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs the Grundmann-Möller quadrature formula.
       * @param s Parameter determining the degree @f$ d = 2s + 1 @f$
       * @param geom Simplex geometry (Segment, Triangle, or Tetrahedron)
       *
       * Creates a Grundmann-Möller quadrature rule for the specified simplex
       * geometry with degree of exactness @f$ 2s + 1 @f$.
       *
       * @note @p s must satisfy @f$ 0 \leq s \leq @f$ 
       * RODIN_QF_GRUNDMANNMOLLER_MAX_S
       * @note @p geom must be a simplicial geometry type
       */
      GrundmannMoller(size_t s, Geometry::Polytope::Type geom);

      /**
       * @brief Copy constructor.
       * @param other GrundmannMoller instance to copy from
       */
      constexpr
      GrundmannMoller(const GrundmannMoller& other)
        : Parent(other),
          m_s(other.m_s),
          m_n(other.m_n),
          m_order(other.m_order)
      {}

      /**
       * @brief Move constructor.
       * @param other GrundmannMoller instance to move from
       */
      constexpr
      GrundmannMoller(GrundmannMoller&& other)
        : Parent(std::move(other)),
          m_s(std::move(other.m_s)),
          m_n(std::move(other.m_n)),
          m_order(std::move(other.m_order))
      {}

      /**
       * @brief Gets the degree of polynomial exactness.
       * @return Degree @f$ d = 2s + 1 @f$ for which this rule is exact
       */
      inline
      size_t getOrder() const
      {
        return m_order;
      }

      /**
       * @brief Gets the number of quadrature points.
       * @return Total number of quadrature points
       */
      size_t getSize() const override;

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight @f$ w_i @f$ associated with the point
       */
      Real getWeight(size_t i) const override;

      /**
       * @brief Gets the coordinates of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Reference to the point coordinates in reference space
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override;

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new GrundmannMoller instance
       */
      inline
      GrundmannMoller* copy() const noexcept override
      {
        return new GrundmannMoller(*this);
      }

      /**
       * @brief Initializes the static size table.
       * @return Multi-dimensional array indexed by [s][n] containing sizes
       */
      static boost::multi_array<size_t, 2> initSizes();

      /**
       * @brief Initializes the static weights table.
       * @return Multi-dimensional array indexed by [s][n] containing weight vectors
       */
      static boost::multi_array<Math::Vector<Real>, 2> initWeights();

      /**
       * @brief Initializes the static points table.
       * @return Multi-dimensional array indexed by [s][n] containing point vectors
       */
      static boost::multi_array<std::vector<Math::SpatialVector<Real>>, 2> initPoints();

    private:
      /// Static table: number of points for each [s][n] combination
      static boost::multi_array<size_t, 2> s_sizes;
      
      /// Static table: quadrature weights for each [s][n] combination
      static boost::multi_array<Math::Vector<Real>, 2> s_weights;
      
      /// Static table: quadrature points for each [s][n] combination
      static boost::multi_array<std::vector<Math::SpatialVector<Real>>, 2> s_points;

      const size_t m_s;     ///< Parameter @f$ s @f$ determining degree @f$ 2s+1 @f$
      const size_t m_n;     ///< Simplex dimension
      const size_t m_order; ///< Degree of exactness @f$ d = 2s + 1 @f$
  };
}

#endif
