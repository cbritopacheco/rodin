/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_DEFAULTCASTS_H
#define RODIN_DEFAULTCASTS_H

#include <vector>
#include <mfem.hpp>

#include "Cast.h"

namespace Rodin
{
   template <class T>
   class Cast<std::vector<T>, false>
   {
      using From = std::vector<T>;
      public:
         Cast(const From& from)
            : m_from(from)
         {}

         const From& from() const
         {
            return m_from;
         }

         template <class To>
         To to() const;

         template <>
         mfem::Array<T> to<mfem::Array<T>>() const
         {
            mfem::Array res(from().size());
            std::copy(from().begin(), from().end(), res.begin());
            return res;
         }
      private:
         const From& m_from;
   };

   template <class T>
   class Cast<std::vector<T>, true>
   {
      using From = std::vector<T>;
      public:
         Cast(From&& from)
            : m_from(std::move(from))
         {}

         From& from()
         {
            return m_from;
         }

         const From& from() const
         {
            return m_from;
         }

         template <class To>
         To to();

         template <>
         mfem::Array<T> to<mfem::Array<T>>()
         {
            int asize = from().size();
            return mfem::Array{from().data(), asize};
         }
      private:
         From m_from;
   };
}

#endif
