/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FORMLANGUAGE_BASE_H
#define RODIN_FORMLANGUAGE_BASE_H

#include <atomic>
#include <vector>
#include <cassert>
#include <variant>
#include <typeinfo>

#include "Rodin/Types.h"
#include "Rodin/Copyable.h"
#include "Rodin/Identifiable.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Traits.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all objects in Rodin's FormLanguage system.
   *
   * This class serves as the foundational base for all form language objects,
   * providing essential services such as unique identification, object lifetime
   * management, and polymorphic copying capabilities. All form language expressions,
   * functions, and operators derive from this base class.
   *
   * @warning The FormLanguage::Base class is not thread safe. Only one thread
   * should access the methods of a FormLanguage object at a time.
   *
   * ## Key Features
   * - **Unique Identification**: Each instance receives a unique UUID for tracking
   * - **Polymorphic Operations**: Support for copying and cloning operations
   */
  class Base : public Copyable, public Identifiable
  {
    public:
      /**
       * @brief Constructor.
       */
      Base() = default;

      /**
       * @brief Copy constructor.
       */
      Base(const Base& other)
        : Copyable(other),
          Identifiable(other)
      {}

      /**
       * @brief Move constructor.
       */
      Base(Base&& other)
        : Copyable(std::move(other)),
          Identifiable(std::move(other))
      {}

      /**
       * @brief Destructor.
       */
      virtual ~Base() = default;

      /**
       * @brief Copy assignment is not allowed.
       */
      Base& operator=(const Base&) = delete;

      /**
       * @brief Move assignment is not allowed.
       */
      Base& operator=(Base&&) = delete;

      /**
       * @brief Gets the human-readable name of this object.
       * 
       * @return const char* Object name string, or nullptr if no name is set
       * 
       * Returns a string representation of the object type or expression.
       * Derived classes should override this method to provide meaningful
       * names for debugging and diagnostic purposes.
       */
      virtual Optional<StringView> getName() const
      {
        return {};
      }

      /**
       * @brief Creates a polymorphic copy of this object.
       * @return Non-owning pointer to the copied object
       *
       * Pure virtual function that must be implemented by derived classes
       * to support polymorphic copying. The returned pointer is non-owning;
       * the caller is responsible for managing its lifetime.
       *
       * @note This is a CRTP function to be overridden in derived classes.
       */
      virtual Base* copy() const noexcept override = 0;

  };
}

#endif
