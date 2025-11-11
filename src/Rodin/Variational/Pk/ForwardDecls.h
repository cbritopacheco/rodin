/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_Pk_FORWARDDECLS_H
#define RODIN_VARIATIONAL_Pk_FORWARDDECLS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @brief Degree k Lagrange element
   * @tparam K Polynomial degree
   * @tparam Range Range value type
   *
   * @note For an overview of all the possible specializations of the
   * Pk class, please see @ref PkSpecializations.
   *
   * @see PkElementSpecializations
   */
  template <size_t K, class Range>
  class PkElement;
}

#endif
