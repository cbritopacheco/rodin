/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TRAITS_FORWARDDECLS_H
#define RODIN_TRAITS_FORWARDDECLS_H

namespace Rodin::Traits
{
   /**
    * @brief Empty tag type to indicate a parallel context.
    * @see Serial
    */
   struct Parallel {};

   /**
    * @brief Empty tag type to indicate a serial context.
    * @see Parallel
    */
   struct Serial {};
}

#endif
