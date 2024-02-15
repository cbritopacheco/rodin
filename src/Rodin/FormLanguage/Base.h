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
#include "Rodin/Threads/Unsafe.h"
#include "Rodin/Utility/DependentFalse.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "Traits.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all classes which are part of Rodin's FormLanguage.
   *
   * @warning The FormLanguage::Base class is not thread safe, i.e. only one
   * thread should have access to the methods of the FormLanguage object at a
   * time.
   */
  class Base : public Copyable
  {
    using ObjectTable = std::vector<std::shared_ptr<const void>>;

    public:
      using UUID = size_t;

      /**
       * @brief Constructor.
       */
      Base()
        : m_uuid(s_id++)
      {
        m_objs.write(
            [](auto& obj)
            {
              obj.reserve(8);
            });
      }

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
       * @brief Gets the unique identifier associated to the instance.
       */
      inline
      const UUID& getUUID() const
      {
        return m_uuid;
      }

      /**
       * @brief Gets the name of the object which it represents.
       */
      virtual const char* getName() const
      {
        return nullptr;
      }

      /**
       * @brief Keeps the passed object in memory for later use.
       */
      template <class T, typename = std::enable_if_t<FormLanguage::IsPlainObject<T>::Value>>
      inline
      constexpr
      const T& object(T&& obj) const noexcept
      {
        using R = typename std::remove_reference_t<T>;
        const R* res = new R(std::forward<T>(obj));
        m_objs.write([&](auto& obj){ obj.emplace_back(res); });;
        return *res;
      }

      /**
       * @brief Returns the same object.
       */
      template <class T, typename = std::enable_if_t<!FormLanguage::IsPlainObject<T>::Value>>
      inline
      constexpr
      T object(T&& obj) const noexcept
      {
        return std::forward<T>(obj);
      }

      /**
       * @brief Destructs the objects stored.
       */
      void clear()
      {
        m_objs.write([](auto& obj){ obj.clear(); });
      }

      /**
       * @brief Copies the object and returns a non-owning pointer to the
       * copied object.
       * @returns Non-owning pointer to the copied object.
       * @note CRTP function to be overriden in the Derived class.
       *
       */
      virtual Base* copy() const noexcept override = 0;

    private:
      thread_local static UUID s_id;

      const size_t m_uuid;
      mutable Threads::Unsafe<std::vector<std::shared_ptr<const void>>> m_objs;
  };
}

#endif
