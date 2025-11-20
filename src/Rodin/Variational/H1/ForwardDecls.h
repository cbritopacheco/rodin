/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_FORWARDDECLS_H
#define RODIN_VARIATIONAL_H1_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @brief Degree k H1-conforming Lagrange element
   * @tparam K Polynomial degree
   * @tparam Range Range value type
   *
   * @note For an overview of all the possible specializations of the
   * H1 class, please see @ref H1ElementSpecializations.
   *
   * @see H1ElementSpecializations
   */
  template <size_t K, class Range>
  class H1Element;

  /**
   * @brief Convenience alias for real-valued H1 element
   */
  template <size_t K>
  using RealH1Element = H1Element<K, Real>;

  /**
   * @brief Convenience alias for complex-valued H1 element
   */
  template <size_t K>
  using ComplexH1Element = H1Element<K, Complex>;

  /**
   * @brief Convenience alias for real vector-valued H1 element
   */
  template <size_t K>
  using VectorH1Element = H1Element<K, Math::Vector<Real>>;

  /**
   * @brief Convenience alias for complex vector-valued H1 element
   */
  template <size_t K>
  using ComplexVectorH1Element = H1Element<K, Math::Vector<Complex>>;
}

#endif
