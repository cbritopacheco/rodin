/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CORE_CONSTANTS_H
#define RODIN_CORE_CONSTANTS_H

#include <cmath>
#include <type_traits>

namespace Rodin::Math::Constants
{
    /**
     * @brief Computes the number @f$ \pi @f$ to machine precision.
     */
    template <class T>
    constexpr
    inline
    std::enable_if_t<std::is_arithmetic_v<T>, T>
    pi()
    {
      return static_cast<T>(std::acos(-1));
    }
}

#endif
