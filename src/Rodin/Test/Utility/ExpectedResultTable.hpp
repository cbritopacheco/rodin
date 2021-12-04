/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TEST_UTILITY_EXPECTEDRESULTTABLE_HPP
#define RODIN_TEST_UTILITY_EXPECTEDRESULTTABLE_HPP

#include <tuple>
#include <vector>
#include <functional>

#include "ExpectedResultTable.h"

namespace Rodin::Test::Utility
{

  template <class Result, class ... Parameters>
  constexpr ExpectedResultTable<Result, Parameters...>::ExpectedResultTable(
      std::function<Result(Parameters&&...)> model,
      std::function<bool(const Result&, const Result&)> compare)
    : m_model(model), m_compare(compare)
  {}

  template <class Result, class ... Parameters>
  constexpr void ExpectedResultTable<Result, Parameters...>::push_back(const ExpectedResult& entry)
  {
    m_table.push_back(entry);
  }

  template <class Result, class ... Parameters>
  constexpr void ExpectedResultTable<Result, Parameters...>::emplace_back(
      Result&& res, Parameters&&... params)
  {
    m_table.push_back(ExpectedResult(
          std::forward<Result>(res), std::forward<Parameters>(params)...));
  }

  template <class Result, class ... Parameters>
  constexpr bool ExpectedResultTable<Result, Parameters...>::evaluate() const
  {
    bool allPassed = true;
    for (const auto& it : m_table)
      allPassed = allPassed && m_compare(
          std::apply(m_model, it.getParameters()), it.getResult());
    return allPassed;
  }
}
#endif
