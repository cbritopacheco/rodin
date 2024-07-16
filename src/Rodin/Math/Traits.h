/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_TRAITS_H
#define RODIN_MATH_TRAITS_H

#include <type_traits>

#include <Eigen/Core>

namespace Rodin::FormLanguage
{
  template <class T>
  struct IsEigenObject
  {
    static constexpr bool Value =
      std::is_base_of_v<Eigen::EigenBase<typename std::decay<T>::type>, typename std::decay<T>::type>;
  };
}

#endif

