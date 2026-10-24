/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_RODIN_KOKKOS_TESTCONTEXT_H
#define RODIN_RODIN_KOKKOS_TESTCONTEXT_H

/**
 * @file
 * @brief Selects the real CUDA backend when enabled and explicit host validation otherwise.
 */

#include <Rodin/Kokkos/Context.h>
#include <Rodin/Types.h>

namespace Rodin::Tests::Unit
{
#ifdef KOKKOS_ENABLE_CUDA
  using TestBackend = Context::CUDA;
#else
  using TestBackend = ::Kokkos::DefaultHostExecutionSpace;
#endif

  using TestContext = Context::Kokkos<Context::Local, TestBackend>;
  using TestExecutionSpace = TestContext::ExecutionSpaceType;
  using TestMemorySpace = TestContext::MemorySpaceType;
}

#endif
