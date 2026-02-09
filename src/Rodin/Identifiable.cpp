/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Identifiable.h"

namespace Rodin
{
  thread_local Identifiable::UUID Identifiable::s_id = 0;
}


