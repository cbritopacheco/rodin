/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FORWARDDECLS_H
#define RODIN_FORWARDDECLS_H

/**
 * @file
 * @brief Forward declarations for core Rodin types.
 *
 * This header provides forward declarations for fundamental template classes
 * used throughout the Rodin library, enabling header file organization and
 * reducing compilation dependencies.
 */

namespace Rodin
{
  /**
   * @brief Forward declaration of the Pair template class.
   * @tparam L Type of the first element.
   * @tparam R Type of the second element.
   * @see Pair
   */
  template <class L, class R>
  class Pair;

  /**
   * @brief Forward declaration of the Tuple template class.
   * @tparam Ts Types of the tuple elements.
   * @see Tuple
   */
  template <class ... Ts>
  class Tuple;
}

#endif
