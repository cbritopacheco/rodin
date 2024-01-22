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
   * @brief Degree 1 Lagrange element
   * @tparam Range Range value type
   *
   * @note For an overview of all the possible specializations of the
   * H1 class, please see @ref H1Specializations.
   *
   * @see H1ElementSpecializations
   */
  template <class Range>
  class H1Element;

  /**
   * @brief Arbitrary order @f$ H^1(\Omega)^d @f$ continuous finite element
   * space.
   * @tparam Context Type of context for the finite element space
   *
   * Given some triangulation @f$ \mathcal{T}_h @f$ of @f$ \Omega @f$,
   * instances of this class will represent the finite element space:
   * @f[
   *   V_h := \left\{ v : \overline{\Omega} \rightarrow \mathbb{R}^d \mid
   *   \forall \tau \in \mathcal{T}_h , \ v|_\tau \in \mathcal{P}_\tau
   *   \right\}
   * @f]
   * where @f$ \mathcal{P}_\tau \subset H^1(\tau)^d @f$ such that @f$ V_h \subset
   * C^0(\Omega)^d @f$.
   *
   */
  template <class Range, class Context, class Mesh>
  class H1;

  /**
   * @ingroup GridFunctionSpecializations
   * @brief GridFunction on the H1 finite element space.
   */
  template <class ... Ts>
  class GridFunction<H1<Ts...>>;

  /**
   * @brief Alias for H1Element<Scalar>
   */
  using ScalarH1Element = H1Element<Scalar>;

  /**
   * @brief Alias for H1Element<Math::Vector>
   */
  using VectorH1Element = H1Element<Math::Vector>;
}

#endif
