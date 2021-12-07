/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CAST_H
#define RODIN_CAST_H

#include <utility>

namespace Rodin
{
   template <class From>
   class Cast
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
}

#endif
