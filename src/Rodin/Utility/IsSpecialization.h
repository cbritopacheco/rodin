/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ISSPECIALIZATION_H
#define RODIN_UTILITY_ISSPECIALIZATION_H

#include <type_traits>

namespace Rodin::Utility
{
  template <class Test, template<class...> class Ref>
  struct IsSpecialization
  {
    static constexpr const bool Value = false;
  };

  template<template<class...> class Ref, class... Args>
  struct IsSpecialization<Ref<Args...>, Ref>
  {
    static constexpr const bool Value = true;
  };
}

#endif
