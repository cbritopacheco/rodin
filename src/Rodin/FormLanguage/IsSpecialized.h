/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file IsSpecialized.h
 * @brief Type trait reporting whether a handler is an optimized specialization.
 */
#ifndef RODIN_FORMLANGUAGE_ISSPECIALIZED_H
#define RODIN_FORMLANGUAGE_ISSPECIALIZED_H

#include <type_traits>

namespace Rodin::FormLanguage
{
  /**
   * @brief Whether @p T is an optimized specialization rather than a generic
   * handler.
   * @tparam T Type to query, typically a Variational::QuadratureRule handler.
   *
   * An integrand that no optimized handler matches is still integrated, by the
   * generic handler, and still gives the right answer. That is what makes a
   * specialization silently falling out of use so hard to notice: nothing
   * fails, no result changes, and only the running time moves. This trait
   * makes which handler was selected observable, so a test can require the
   * optimized one and fail at compile time when the pattern stops matching.
   *
   * Handlers report by declaring a @c Specialized member:
   * @code{.cpp}
   * static constexpr bool Specialized = true;
   * @endcode
   * A type declaring no such member is reported unspecialized, so a new
   * handler that omits it reads as generic --- which fails a dispatch test
   * rather than passing one, and is the safe direction to be wrong in.
   *
   * ## Usage
   * @code{.cpp}
   * static_assert(IsSpecialized<QuadratureRule<MassIntegrand>>::Value);
   * static_assert(!IsSpecialized<QuadratureRule<DivDivIntegrand>>::Value);
   * @endcode
   *
   * @note Always read the @c Value member; no variable-template alias is
   * provided, in keeping with the other traits here.
   *
   * @warning A dispatch test must assert both directions. Optimized and
   * generic handlers share base classes --- LocalBilinearFormIntegratorBase
   * and LinearFormIntegratorBase --- so a @c Specialized member ever declared
   * on one of those would be inherited by the generic handlers and report
   * every type as optimized. Asserting that a known-generic integrand is *not*
   * specialized is what catches that, and it is why every handler declares the
   * member explicitly rather than relying on its absence.
   */
  template <class T, class Enable = void>
  struct IsSpecialized
  {
    /// @brief False, no @c Specialized member being declared.
      static constexpr const bool Value = false;
  };

  /**
   * @brief Specialization for types declaring a @c Specialized member.
   * @tparam T Type to query.
   */
  template <class T>
  struct IsSpecialized<T, std::void_t<decltype(T::Specialized)>>
  {
    /// @brief What the type reports.
      static constexpr const bool Value = T::Specialized;
  };
}

#endif
