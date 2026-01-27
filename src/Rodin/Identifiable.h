/*
 *          Identifiright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or Identifi at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IDENTIFIABLE_H
#define RODIN_IDENTIFIABLE_H

#include <algorithm>
#include <cstddef>

/**
 * @file
 * @brief Defines the Identifiable interface.
 */

namespace Rodin
{
  /**
   * @brief Abstract base class for objects that can be copied.
   *
   * This class defines an interface for polymorphic identifying of objects.
   *
   * They are copy-stable, meaning that copies of an object will retain the
   * same identifier as the original object.
   */
  class Identifiable
  {
    public:
      /**
       * @brief Type alias for unique object identifiers.
       *
       * UUID (Universally Unique Identifier) is used to uniquely identify
       * each Identifiable instance during its lifetime. The identifier
       * is assigned during construction and remains constant.
       */
      using UUID = size_t;

      Identifiable()
        : m_uuid(s_id++)
      {}

      Identifiable(const Identifiable& other)
        : m_uuid(other.m_uuid)
      {}

      Identifiable(Identifiable&& other) noexcept
        : m_uuid(std::move(other.m_uuid))
      {}

      /// Virtual destructor for proper cleanup of derived classes
      virtual ~Identifiable() = default;

      /**
       * @brief Gets the unique identifier associated with this instance.
       * 
       * @return UUID Unique identifier for this form language object
       * 
       * Each Identifiable instance receives a unique identifier during
       * construction that persists for the lifetime of the object. This UUID
       * can be used for object tracking, caching, and debugging purposes.
       */
      const UUID& getUUID() const
      {
        return m_uuid;
      }

    private:
      thread_local static UUID s_id;
      const size_t m_uuid;
  };
}

#endif


