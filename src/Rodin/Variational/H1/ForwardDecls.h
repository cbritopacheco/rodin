/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_FORWARDDECLS_H
#define RODIN_VARIATIONAL_H1_FORWARDDECLS_H

/**
 * @file
 * @brief Forward declarations and type aliases for H1Element classes.
 */

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Variational
{
  /**
   * @brief Degree k H1-conforming Lagrange element.
   * @tparam K Polynomial degree
   * @tparam Range Range value type (Scalar or Math::Vector<Scalar>)
   *
   * H1Element provides a high-order H1-conforming finite element implementation
   * with arbitrary polynomial degree K. The element uses Lagrange nodal basis
   * functions satisfying @f$ \phi_i(x_j) = \delta_{ij} @f$.
   *
   * @note For an overview of all the possible specializations of the
   * H1 class, please see @ref H1ElementSpecializations.
   *
   * @see H1ElementSpecializations
   */
  template <size_t K, class Range>
  class H1Element;

  /**
   * @brief Convenience alias for real-valued H1 element.
   * @tparam K Polynomial degree
   *
   * Equivalent to `H1Element<K, Real>`. Used for scalar-valued problems
   * such as the Poisson equation or heat conduction.
   */
  template <size_t K>
  using RealH1Element = H1Element<K, Real>;

  /**
   * @brief Convenience alias for complex-valued H1 element.
   * @tparam K Polynomial degree
   *
   * Equivalent to `H1Element<K, Complex>`. Used for complex-valued problems
   * such as wave equations in frequency domain.
   */
  template <size_t K>
  using ComplexH1Element = H1Element<K, Complex>;

  /**
   * @brief Convenience alias for real vector-valued H1 element.
   * @tparam K Polynomial degree
   *
   * Equivalent to `H1Element<K, Math::Vector<Real>>`. Used for vector-valued
   * problems such as linear elasticity or fluid mechanics.
   */
  template <size_t K>
  using VectorH1Element = H1Element<K, Math::Vector<Real>>;

  /**
   * @brief Convenience alias for complex vector-valued H1 element.
   * @tparam K Polynomial degree
   *
   * Equivalent to `H1Element<K, Math::Vector<Complex>>`. Used for complex
   * vector-valued problems such as electromagnetic wave propagation.
   */
  template <size_t K>
  using ComplexVectorH1Element = H1Element<K, Math::Vector<Complex>>;
}

#endif
