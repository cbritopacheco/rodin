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
#include <memory>
#include <cassert>
#include <variant>
#include <typeinfo>

#include "Rodin/Types.h"
#include "Rodin/Copyable.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Traits.h"
#include "IsPlaneObject.h"

namespace Rodin::FormLanguage
{
  /**
   * @defgroup RodinFormLanguage Form Language Module
   * @brief Domain-specific language for expressing finite element forms and expressions.
   *
   * The FormLanguage module provides a comprehensive framework for expressing
   * variational forms, mathematical expressions, and finite element operations
   * using a high-level, expressive domain-specific language. This module forms
   * the foundation for Rodin's symbolic computation capabilities.
   */

  /**
   * @ingroup RodinFormLanguage
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
   * - **Object Management**: Automatic lifetime management for temporary objects
   * - **Polymorphic Operations**: Support for copying and cloning operations
   * - **Type Safety**: Template-based object storage with type validation
   */
  class Base : public Copyable
  {
    using ObjectTable = std::vector<std::shared_ptr<const void>>;

    public:
      /**
       * @brief Type alias for unique object identifiers.
       *
       * UUID (Universally Unique Identifier) is used to uniquely identify
       * each FormLanguage::Base instance during its lifetime. The identifier
       * is assigned during construction and remains constant.
       */
      using UUID = size_t;

      /**
       * @brief Constructor.
       */
      Base()
        : m_uuid(s_id++)
      {}

      /**
       * @brief Copy constructor.
       */
      Base(const Base& other) = default;

      /**
       * @brief Move constructor.
       */
      Base(Base&&) = default;

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
       * @brief Gets the unique identifier associated with this instance.
       * 
       * @return UUID Unique identifier for this form language object
       * 
       * Each FormLanguage::Base instance receives a unique identifier during
       * construction that persists for the lifetime of the object. This UUID
       * can be used for object tracking, caching, and debugging purposes.
       */
      const UUID& getUUID() const
      {
        return m_uuid;
      }

      /**
       * @brief Gets the human-readable name of this object.
       * 
       * @return const char* Object name string, or nullptr if no name is set
       * 
       * Returns a string representation of the object type or expression.
       * Derived classes should override this method to provide meaningful
       * names for debugging and diagnostic purposes.
       */
      virtual const char* getName() const
      {
        return nullptr;
      }

      /**
       * @brief Stores an object for automatic lifetime management.
       * 
       * @tparam T Type of object to store (must be a plain object type)
       * @param obj Object to store (rvalue) or reference (lvalue) 
       * @return const T& Reference to the stored object
       * 
       * This method provides automatic lifetime management for temporary objects
       * used in form language expressions. For rvalue references, the object is
       * moved into internal storage and its lifetime is tied to this Base instance.
       * For lvalue references, the original object is returned unchanged.
       * 
       * @note Only plain object types (as defined by IsPlainObject) are accepted
       */
      template <class T, typename =
        std::enable_if_t<FormLanguage::IsPlainObject<std::remove_reference_t<T>>::Value>>
      constexpr
      const T& object(T&& obj) const noexcept
      {
        if constexpr (std::is_lvalue_reference_v<T>)
        {
          return obj;
        }
        else
        {
          using R = typename std::remove_reference_t<T>;
          const R* res = new R(std::forward<T>(obj));
          m_objs.emplace_back(res);
          return *res;
        }
      }

      /**
       * @brief Forwards non-plain objects unchanged.
       * @tparam T Type of object (must not be a plain object type)
       * @param[in] obj Object to forward
       * @return Forwarded object
       *
       * This overload handles non-plain object types (such as scalars, references,
       * or expression templates) by forwarding them directly without storage.
       * It is selected via SFINAE when T is not a plain object.
       */
      template <class T, typename =
        std::enable_if_t<!FormLanguage::IsPlainObject<std::remove_reference_t<T>>::Value>>
      constexpr
      T object(T&& obj) const noexcept
      {
        return std::forward<T>(obj);
      }

      /**
       * @brief Clears all stored objects, releasing their memory.
       *
       * Destroys all objects that were stored via the object() method,
       * freeing the associated memory. This is useful for managing
       * temporary object lifetimes explicitly.
       *
       * @note After calling clear(), any references obtained from previous
       * object() calls become invalid.
       */
      void clear()
      {
        m_objs.clear();
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

    private:
      thread_local static UUID s_id;

      const size_t m_uuid;
      mutable std::vector<std::shared_ptr<const void>> m_objs;
  };
}

#endif
