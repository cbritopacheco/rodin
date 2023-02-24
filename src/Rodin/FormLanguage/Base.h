/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FORMLANGUAGE_BASE_H
#define RODIN_FORMLANGUAGE_BASE_H

#include <deque>
#include <memory>
#include <cassert>
#include <variant>
#include <typeinfo>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "Rodin/Types.h"
#include "Rodin/Utility/DependentFalse.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/BasisOperator.h"

#include "Traits.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all classes which are part of Variational::FormLanguage.
   */
  class Base
  {
    static boost::uuids::random_generator s_gen;

    public:
      Base()
        : m_uuid(s_gen())
      {}

      Base(const Base& other)
        : m_uuid(other.m_uuid)
      {}

      Base(Base&&) = default;

      Base& operator=(const Base&) = delete;

      Base& operator=(Base&&) = delete;

      inline
      const boost::uuids::uuid& getUUID() const
      {
        return m_uuid;
      }

      /**
       * @brief Virtual destructor.
       */
      virtual ~Base() = default;

      virtual const char* getName() const
      {
        return "Rodin::FormLanguage::Base";
      }

      template <class T, typename = std::enable_if_t<FormLanguage::IsPlainObject<T>::Value>>
      inline
      constexpr
      const T& object(T&& obj) const
      {
        T* res = new T(std::forward<T>(obj));
        m_objs.emplace_back(res);
        return *res;
      }

      template <class T, typename = std::enable_if_t<!FormLanguage::IsPlainObject<T>::Value>>
      inline
      constexpr
      T&& object(T&& obj) const
      {
        return std::forward<T>(obj);
      }

      void clear()
      {
        m_objs.clear();
      }

      /**
       * @internal
       * @brief Copies the object and returns a non-owning pointer to the
       * copied object.
       * @returns Non-owning pointer to the copied object.
       */
      virtual Base* copy() const noexcept = 0;

    private:
      const boost::uuids::uuid m_uuid;
      mutable std::deque<std::shared_ptr<void>> m_objs;
  };
}

#endif
