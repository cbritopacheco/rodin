/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_MODEL_HISTORY_H
#define RODIN_HEART_CCMLC2014_MODEL_HISTORY_H

#include <optional>

namespace Rodin::Heart::CCMLC2014::Model
{
  template <class State>
  struct HistoryT
  {
    State n;
    State nm1;
    std::optional<State> nm2;
  };
}

#endif
