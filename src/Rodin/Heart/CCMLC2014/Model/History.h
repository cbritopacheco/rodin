/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file History.h
 * @brief Time-history container used by the CCMLC2014 stepper.
 */
#ifndef RODIN_HEART_CCMLC2014_MODEL_HISTORY_H
#define RODIN_HEART_CCMLC2014_MODEL_HISTORY_H

#include <optional>

namespace Rodin::Heart::CCMLC2014::Model
{
  /**
   * @brief Storage of recent states needed by second-order finite differences.
   *
   * @tparam State State type.
   */
  template <class State>
  struct HistoryT
  {
    State n; ///< State at level @f$ n @f$.
    State nm1; ///< State at level @f$ n-1 @f$.
    std::optional<State> nm2; ///< Optional state at level @f$ n-2 @f$.
  };
}

#endif
