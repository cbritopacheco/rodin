/**
 * @file Debug.h
 * @brief Debugging utilities for form language objects.
 *
 * This file provides template utilities for debugging and introspecting
 * form language objects during development. The Debug template can be
 * specialized to provide type-specific debugging information.
 */
#ifndef RODIN_FORMLANGUAGE_DEBUG_H
#define RODIN_FORMLANGUAGE_DEBUG_H

namespace Rodin::FormLanguage
{
  /**
   * @brief Template for debugging form language types.
   * @tparam T Type to debug
   * @ingroup RodinFormLanguage
   *
   * The Debug template provides a mechanism for introspecting and debugging
   * form language types during compilation and runtime. Specialize this template
   * for specific types to provide custom debugging information.
   *
   * ## Usage
   * This template is primarily used during development to understand type
   * transformations and check template instantiations. When instantiated,
   * compiler errors will reveal the actual type of T.
   *
   * Example usage:
   * @code
   * // This will trigger a compile error showing the actual type
   * Debug<decltype(myExpression)> debug;
   * @endcode
   */
  template <typename> struct Debug;
}

#endif
