/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_EXTRACT_H
#define RODIN_UTILITY_EXTRACT_H

/**
 * @file
 * @brief Defines the Extract metafunction for tuple type transformation.
 */

#include "Rodin/Tuple.h"

/**
 * @defgroup UtilityModule Utility
 * @brief Template metaprogramming utilities and type manipulation tools.
 *
 * The Utility module provides a comprehensive set of template metaprogramming
 * utilities for type manipulation, parameter pack handling, and compile-time
 * computations used throughout the Rodin finite element library.
 */

namespace Rodin::Utility
{
  /**
   * @brief Template metafunction for extracting types from tuple structures.
   * @ingroup UtilityModule
   * @tparam Ts Template parameter pack (specialized for specific tuple types).
   *
   * The Extract template provides a mechanism to apply a type transformation
   * (Extractor) to each element type in a tuple, producing a new tuple with
   * the transformed types.
   *
   * This is useful for type-level transformations such as adding const qualifiers,
   * extracting nested types, or applying other template metafunctions to tuple elements.
   */
  template <class ...>
  struct Extract;

  /**
   * @brief Specialization of Extract for Tuple types.
   * @ingroup UtilityModule
   * @tparam T The first type in the tuple.
   * @tparam Ts The remaining types in the tuple.
   *
   * This specialization enables the extraction and transformation of types
   * from a Rodin::Tuple by applying an Extractor metafunction to each element type.
   *
   * Example usage:
   * @code{.cpp}
   * template<typename T>
   * struct AddConst { using Type = const T; };
   * 
   * using OriginalTuple = Tuple<int, double, std::string>;
   * using ConstTuple = Extract<OriginalTuple>::Type<AddConst>;
   * // ConstTuple is Tuple<const int, const double, const std::string>
   * @endcode
   */
  template <class T, class ... Ts>
  struct Extract<Tuple<T, Ts...>>
  {
    /**
     * @brief Applies an extractor metafunction to all tuple element types.
     * @tparam Extractor A template metafunction that transforms a type T to Extractor<T>::Type.
     *
     * This type alias creates a new Tuple where each element type has been
     * transformed by the provided Extractor metafunction.
     */
    template <template <class> class Extractor>
    using Type = Tuple<typename Extractor<T>::Type, typename Extractor<Ts>::Type...>;
  };
}

#endif

