/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_FORWARDDECLS_H
#define RODIN_VARIATIONAL_P1_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Variational
{

  /**
   * @brief Degree 1 Lagrange element
   * @tparam Range Range value type
   *
   * @note For an overview of all the possible specializations of the
   * P1 class, please see @ref P1Specializations.
   *
   * @see P1ElementSpecializations
   */
  template <class Range>
  class P1Element;

  /**
   * @brief Degree 1 Lagrange finite element space
   * @tparam Range Range value type
   * @tparam Context Context type
   * @tparam Args Additional arguments
   *
   * Represents the finite element space composed of continuous, piecewise
   * linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h)^d = \{ v \in C^0(\mathcal{T}_h)^d \mid v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ ,
   * @f]
   * for a given vector dimension @f$ d \in \mathbb{N} @f$.
   *
   * @note For an overview of all the possible specializations of the
   * P1 class, please see @ref P1Specializations.
   *
   * @see P1Specializations
   */
  template <class Range, class Context, class Mesh>
  class P1;

  /**
   * @ingroup GridFunctionSpecializations
   * @brief GridFunction on the P1 finite element space.
   */
  template <class ... Ts>
  class GridFunction<P1<Ts...>>;

  /**
   * @brief Alias for P1Element<Scalar>
   */
  using ScalarP1Element = P1Element<Scalar>;

  /**
   * @brief Alias for P1Element<Math::Vector>
   */
  using VectorP1Element = P1Element<Math::Vector>;
}

#endif
