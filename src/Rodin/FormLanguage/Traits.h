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
 * @code
 * // Query traits of a form language object
 * using MyFES = typename Traits<MyFormObject>::FESType;
 * constexpr auto space = Traits<MyFormObject>::SpaceType;
 * @endcode
 */
#ifndef RODIN_FORMLANGUAGE_TRAITS_H
#define RODIN_FORMLANGUAGE_TRAITS_H

#include <type_traits>
#include <boost/type_index.hpp>

#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /**
   * @defgroup TraitsSpecializations Traits Template Specializations
   * @brief Template specializations of the Traits class.
   * @ingroup RodinFormLanguage
   *
   * This group contains all specializations of the Traits template for
   * different form language types. Each specialization provides type
   * information specific to that form language construct.
   *
   * @see Traits
   */

  /**
   * @brief Primary template for type traits of form language objects.
   * @tparam T Form language type to query traits from
   * @tparam Enable SFINAE enabler for conditional specializations
   * @ingroup RodinFormLanguage
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
   * @code
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
}

#endif
