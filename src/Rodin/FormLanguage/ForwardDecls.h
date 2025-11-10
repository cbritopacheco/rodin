/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the FormLanguage module.
 *
 * This file provides forward declarations for the main types in the FormLanguage
 * module, enabling reduced compilation dependencies and supporting recursive type
 * definitions.
 */
#ifndef RODIN_FORMLANGUAGE_FORWARDDECLS_H
#define RODIN_FORMLANGUAGE_FORWARDDECLS_H

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all form language objects.
   * @see Base
   */
  class Base;

  /**
   * @brief Container for polymorphic form language objects.
   * @tparam T Type of elements stored in the list
   * @see List
   */
  template <class T>
  class List;

  /**
   * @brief Type traits for form language objects.
   * @tparam Args Template parameters for trait specializations
   * @see Traits
   */
  template <class ... Args>
  struct Traits;
}

#endif
