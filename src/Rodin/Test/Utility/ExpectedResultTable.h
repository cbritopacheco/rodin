/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TEST_UTILITY_EXPECTEDRESULTTABLE_H
#define RODIN_TEST_UTILITY_EXPECTEDRESULTTABLE_H

#include <tuple>
#include <vector>
#include <functional>

namespace Rodin::Test::Utility
{
  template <class Result, class ... Parameters>
  class ExpectedResultTable
  {
   public:
    class ExpectedResult
    {
      public:
       constexpr ExpectedResult(const Result& res, const Parameters&... params)
        : m_res(res), m_params{params...}
       {}

       constexpr Result getResult() const
       {
        return m_res;
       }

       constexpr std::tuple<Parameters...> getParameters() const
       {
        return m_params;
       }

      private:
       Result m_res;
       std::tuple<Parameters...> m_params;
    };
    constexpr ExpectedResultTable(
       std::function<Result(Parameters&&...)> model,
       std::function<bool(const Result&, const Result&)> compare =
       [](const Result& modelResult, const Result& knownResult){
        return modelResult == knownResult;
       });

    constexpr void push_back(const ExpectedResult& entry);
    constexpr void emplace_back(Result&& res, Parameters&&... params);
    constexpr bool evaluate() const;

   private:
    std::vector<ExpectedResult> m_table;
    std::function<Result(Parameters&&...)> m_model;
    std::function<bool(const Result&, const Result&)> m_compare;
  };
}

#include "ExpectedResultTable.hpp"

#endif
