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
  /**
   * @brief Table-driven testing utility for validating computational models.
   *
   * This class provides a framework for table-driven testing where expected 
   * results are stored alongside their corresponding input parameters. It 
   * allows systematic validation of computational models against known results.
   *
   * @tparam Result Type of the expected result
   * @tparam Parameters Types of the input parameters
   */
  template <class Result, class ... Parameters>
  class ExpectedResultTable
  {
   public:
    /**
     * @brief Container for a single expected result and its parameters.
     */
    class ExpectedResult
    {
      public:
       /**
        * @brief Constructs an expected result entry.
        * @param res Expected result value
        * @param params Input parameters that should produce this result
        */
       constexpr ExpectedResult(const Result& res, const Parameters&... params)
        : m_res(res), m_params{params...}
       {}

       /**
        * @brief Gets the expected result value.
        * @return The expected result
        */
       constexpr Result getResult() const
       {
        return m_res;
       }

       /**
        * @brief Gets the input parameters.
        * @return Tuple of input parameters
        */
       constexpr std::tuple<Parameters...> getParameters() const
       {
        return m_params;
       }

      private:
       Result m_res;                        ///< Expected result value
       std::tuple<Parameters...> m_params; ///< Input parameters
    };

    /**
     * @brief Constructs an expected result table.
     * @param model Function/model to test
     * @param compare Comparison function for results (defaults to equality)
     */
    constexpr ExpectedResultTable(
       std::function<Result(Parameters&&...)> model,
       std::function<bool(const Result&, const Result&)> compare =
       [](const Result& modelResult, const Result& knownResult){
        return modelResult == knownResult;
       });

    /**
     * @brief Adds an expected result entry to the table.
     * @param entry Expected result entry to add
     */
    constexpr void push_back(const ExpectedResult& entry);

    /**
     * @brief Constructs and adds an expected result entry to the table.
     * @param res Expected result value
     * @param params Input parameters
     */
    template <class R, class ... Ps>
    constexpr void emplace_back(R&& res, Ps&&... params)
    {
      m_table.push_back(ExpectedResult(std::forward<R>(res), std::forward<Parameters>(params)...));
    }

    /**
     * @brief Evaluates all entries in the table against the model.
     * @return True if all expected results match model outputs, false otherwise
     */
    constexpr bool evaluate() const;

   private:
    std::vector<ExpectedResult> m_table;                             ///< Table of expected results
    std::function<Result(Parameters&&...)> m_model;                 ///< Model function to test
    std::function<bool(const Result&, const Result&)> m_compare;    ///< Result comparison function
  };
}

#include "ExpectedResultTable.hpp"

#endif
