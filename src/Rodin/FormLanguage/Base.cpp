/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Base.h"

namespace Rodin::FormLanguage
{
   boost::uuids::random_generator Base::s_gen;

   Base::Base()
      : m_uuid(s_gen())
   {}

   Base::Base(const Base& other)
      : m_uuid(other.m_uuid)
   {}

   Base::Base(Base&& other)
      : m_uuid(std::move(other.m_uuid))
   {}

   const boost::uuids::uuid& Base::getUUID() const
   {
      return m_uuid;
   }
}

