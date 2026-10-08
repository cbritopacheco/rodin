/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Forward declarations for the SWIFT fitting API.
 */
#ifndef RODIN_ADAPTATION_SWIFT_FORWARDDECLS_H
#define RODIN_ADAPTATION_SWIFT_FORWARDDECLS_H

#include <type_traits>
#include <utility>

namespace Rodin::Adaptation::SWIFT
{
  struct Parameters;
  struct Report;
  template <class TrialFunctionType, class TestFunctionType>
  class Problem;
  template <class Mesh,
    class ContextType =
      std::remove_cvref_t<decltype(std::declval<const Mesh&>().getContext())>>
  class Adapt;
}

#endif
