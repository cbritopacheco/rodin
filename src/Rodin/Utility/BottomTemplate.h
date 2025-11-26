/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_BOTTOMTEMPLATE_H
#define RODIN_UTILITY_BOTTOMTEMPLATE_H

/**
 * @file
 * @brief Defines the BottomTemplate metafunction for extracting innermost types.
 */

namespace Rodin::Utility
{
  /**
   * @brief Metafunction to extract the innermost (bottom) type from nested templates.
   * @ingroup UtilityModule
   * @tparam T The type to extract the bottom type from.
   *
   * BottomTemplate recursively unwraps template wrappers to find the
   * innermost type. For a non-template type T, the result is T itself.
   * For a nested template like A<B<C<D>>>, the result is D.
   *
   * Example usage:
   * @code
   * using Nested = std::vector<std::optional<int>>;
   * using Bottom = BottomTemplate<Nested>::Type;
   * // Bottom is int (assuming single-param template specialization)
   * @endcode
   */
  template <class T>
  struct BottomTemplate
  {
    using Type = T;  ///< The type itself when T is not a template instantiation
  };

  /**
   * @brief Specialization for single-parameter template instantiations.
   * @ingroup UtilityModule
   * @tparam T Template with a single type parameter.
   * @tparam S The type parameter of T.
   *
   * Recursively extracts the bottom type from the inner type S.
   */
  template <template <class> class T, class S>
  struct BottomTemplate<T<S>>
  {
    using Type = typename BottomTemplate<S>::Type;  ///< Recursively extract from inner type
  };
}

#endif
