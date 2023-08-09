/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <type_traits>

#include "Base.h"

namespace Rodin::FormLanguage
{
  thread_local Base::UUID Base::s_id = 0;
}

