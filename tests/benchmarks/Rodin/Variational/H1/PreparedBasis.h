/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_BENCHMARKS_H1_PREPAREDBASIS_H
#define RODIN_TESTS_BENCHMARKS_H1_PREPAREDBASIS_H

/**
 * @file
 * @brief Registers H1 basis consumption with explicitly prepared reference data.
 *
 * The shared Variational basis-production fixture owns the fields, numerical
 * guards and measured sweep. This component selects its prepared-reference
 * protocol; it does not duplicate basis evaluation or own a reference cache.
 */
namespace Rodin::Tests::Benchmarks
{
  /// @brief Registers degrees one through three in one compilation unit.
  void registerH1PreparedBasisLow();

  /// @brief Registers degrees four through six in one compilation unit.
  void registerH1PreparedBasisHigh();
}

#endif
