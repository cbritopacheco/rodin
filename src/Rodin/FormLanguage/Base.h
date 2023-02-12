/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_FORMLANGUAGE_BASE_H
#define RODIN_FORMLANGUAGE_BASE_H

#include <memory>
#include <cassert>
#include <typeinfo>
#include <sstream>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

namespace Rodin::FormLanguage
{
  /**
   * @brief Base class for all classes which are part of Variational::FormLanguage.
   */
  class Base
  {
    static boost::uuids::random_generator s_gen;

    public:
      inline
      Base()
        : m_uuid(s_gen())
      {}

      inline
      Base(const Base& other)
        : m_uuid(other.m_uuid)
      {}

      inline
      Base(Base&& other)
        : m_uuid(std::move(other.m_uuid))
      {}

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

      /**
       * @internal
       * @brief Copies the object and returns a non-owning pointer to the
       * copied object.
       * @returns Non-owning pointer to the copied object.
       */
      virtual Base* copy() const noexcept = 0;

    private:
      const boost::uuids::uuid m_uuid;
  };
}

#endif
