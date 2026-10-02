/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Traits.h
 * @brief Type traits system for form language objects.
 *
 * This file defines the Traits template, which provides a uniform interface
 * for querying type information and properties of form language objects.
 * Specialized versions of Traits are defined throughout the library for
 * specific form language types.
 *
 * ## Design Pattern
 * The Traits system follows the trait-based design pattern common in C++
 * template metaprogramming. Each form language type can specialize the
 * Traits template to expose:
 * - Associated types (e.g., scalar type, function space type)
 * - Compile-time properties (e.g., space type, dimension)
 * - Type transformations
 *
 * ## Usage
 * @code{.cpp}
 * // Query traits of a form language object
 * using MyFES = typename Traits<MyFormObject>::FESType;
 * constexpr auto space = Traits<MyFormObject>::SpaceType;
 * @endcode
 */
#ifndef RODIN_FORMLANGUAGE_TRAITS_H
#define RODIN_FORMLANGUAGE_TRAITS_H

#include <type_traits>
#include <utility>
#include <boost/type_index.hpp>

#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Primary template for type traits of form language objects.
   * @tparam T Form language type to query traits from
   * @tparam Enable SFINAE enabler for conditional specializations
   *
   * The Traits template provides a centralized mechanism for querying
   * compile-time type information about form language objects. This is
   * the primary (unspecialized) template; specific types require explicit
   * specializations to be useful.
   *
   * ## Common Trait Members
   * While the primary template is intentionally incomplete, specializations
   * typically provide some or all of the following:
   *
   * ### Type Members
   * - **FESType**: Associated finite element space type
   * - **ScalarType**: Scalar type for numerical values
   * - **RangeType**: Range type for function evaluations
   *
   * ### Static Constants
   * - **SpaceType**: Indicates trial space, test space, or other
   * - **Dimension**: Spatial or value dimension
   *
   * ## Specialization Example
   * @code{.cpp}
   * template <class Derived, class FES>
   * struct Traits<MyFormType<Derived, FES>>
   * {
   *   using FESType = FES;
   *   using ScalarType = typename FES::ScalarType;
   *   static constexpr auto SpaceType = ShapeFunctionSpaceType::Trial;
   * };
   * @endcode
   *
   * @note This primary template is intentionally left incomplete. Users must
   * specialize Traits for their specific types. Attempting to use the primary
   * template directly will result in a compilation error.
   */
  template <class T, class Enable = void>
  struct Traits;

  /**
   * @brief Derived type @c D of the Variational::FunctionBase<D> that @c T
   * inherits from.
   *
   * Function expressions do not name themselves as their FunctionBase
   * parameter: a RealFunction, for one, derives from a FunctionBase over a
   * ScalarFunctionBase. Anything that stores a function through its
   * FunctionBase, the way the form language's own operators do, needs this
   * type rather than the function's.
   */
  template <class T>
  struct FunctionDerived
  {
    private:
      template <class D>
      static D deduce(const Variational::FunctionBase<D>&);

    public:
      /// @brief The FunctionBase parameter of @c T.
      using Type = decltype(deduce(std::declval<const T&>()));
  };

  /**
   * @brief Type trait: whether @c T is a named bilinear form.
   *
   * A named form names its integrand instead of building it through the form
   * language and brings its own local kernel, which lets it be assembled
   * through Assembly::ScatterMap. Each named form specializes this trait next
   * to its declaration.
   *
   * The trait is matched on the form's template rather than detected from
   * its members because the assembly is selected while the form's own class
   * definition is still being instantiated, when none of its members can be
   * looked up yet.
   */
  template <class T>
  struct IsNamedForm : std::false_type
  {
      /// @brief False for types that are not named forms.
      static constexpr bool Value = false;
  };
}

#endif
