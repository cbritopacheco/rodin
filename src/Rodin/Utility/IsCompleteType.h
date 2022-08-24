/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_ISCOMPLETETYPE_H
#define RODIN_UTILITY_ISCOMPLETETYPE_H

#include <cstdlib>
#include <type_traits>

namespace Rodin::Utility
{
   namespace Internal
   {
      template <class T, std::size_t = sizeof(T)>
      std::true_type IsCompleteTypeImpl(T *);

      std::false_type IsCompleteTypeImpl(...);
   }

   template <class T>
   using IsCompleteType = decltype(Internal::IsCompleteTypeImpl(std::declval<T*>()));

}

#endif
