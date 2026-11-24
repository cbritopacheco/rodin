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
 * @brief Selects the real CUDA or HIP backend when enabled and explicit host validation otherwise,
 * and ends the test processes.
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <Rodin/Kokkos/Context.h>
#include <Rodin/Types.h>

namespace Rodin::Tests::Unit
{
#ifdef KOKKOS_ENABLE_CUDA
  using TestBackend = Context::CUDA;
#elif defined(KOKKOS_ENABLE_HIP)
  using TestBackend = Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
  using TestBackend = Context::SYCL;
#else
  using TestBackend = ::Kokkos::DefaultHostExecutionSpace;
#endif

  using TestContext = Context::Kokkos<Context::Local, TestBackend>;
  using TestExecutionSpace = TestContext::ExecutionSpaceType;
  using TestMemorySpace = TestContext::MemorySpaceType;

  /**
   * @brief Ends a test process whose Kokkos, PETSc and MPI runtimes are
   * finalized.
   *
   * When several processes share an AMD GPU, the static destructors of the
   * HIP runtime library sometimes block forever after main has returned,
   * with or without Rodin code on the stack. With HIP, the output is
   * flushed and the process ends without running them.
   * @param status Exit status of the process.
   */
  [[noreturn]] inline void exitTestProcess(int status)
  {
    std::cout.flush();
    std::cerr.flush();
    std::fflush(nullptr);
#ifdef KOKKOS_ENABLE_HIP
    std::_Exit(status);
#else
    std::exit(status);
#endif
  }
}

#endif
