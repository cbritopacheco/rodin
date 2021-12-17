/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_OVERLOADED_H
#define RODIN_UTILITY_OVERLOADED_H

namespace Rodin::Utility
{
   /**
    * @brief Helper type for use with visitor pattern.
    */
   template <class... Ts>
   struct Overloaded : Ts...
   {
      using Ts::operator()...;
   };

   template<class... Ts> Overloaded(Ts...) -> Overloaded<Ts...>;
}

#endif
