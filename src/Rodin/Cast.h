/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CAST_H
#define RODIN_CAST_H

#include <variant>
#include <utility>

#include "Rodin/Alert.h"
#include "Rodin/Utility/Overloaded.h"

namespace Rodin
{
   template <class From, bool Owner = false>
   class Cast;

   template <class From>
   using MutableCast = Cast<From, true>;

   template <class From>
   Cast(From&) -> Cast<From, false>;

   template <class From>
   Cast(const From&) -> Cast<From, false>;

   template <class From>
   class Cast<From, false>
   {
      public:
         Cast(const From& from)
            : m_from(from)
         {}

         const From& from() const
         {
            return m_from;
         }

         template <class To>
         To to() const
         {
            return static_cast<To>(from());
         }

      private:
         const From& m_from;
   };

   template <class From>
   Cast(From&&) -> Cast<From, true>;

   template <class From>
   class Cast<From, true>
   {
      public:
         Cast(From&& from)
            : m_from(std::move(from))
         {}

         const From& from() const
         {
            return m_from;
         }

         From& from()
         {
            return m_from;
         }

         template <class To>
         To to() const
         {
            return static_cast<To>(from());
         }

      private:
         From m_from;
   };
}

#endif
