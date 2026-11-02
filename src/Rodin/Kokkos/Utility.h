/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_UTILITY_H
#define RODIN_KOKKOS_UTILITY_H

/**
 * @file
 * @brief Transfers at Kokkos construction and host query boundaries.
 */

#include <algorithm>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>
#include <Kokkos_Core.hpp>
#include "Rodin/Configure.h"
#include "Rodin/Alert/Exception.h"

/**
 * @def RODIN_KOKKOS_UNROLL
 * @brief Unrolls the following loop completely when its trip count is static.
 *
 * Kernel loops over the local basis functions have a static trip count when
 * the polytope and the spaces are known at compile time. Unrolled, their local
 * indices are constants, so the per-component structure of vector bases folds
 * away and the local arrays stay in registers. The pragma is a no-op for a
 * dynamic trip count and on compilers that do not support it.
 */
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#define RODIN_KOKKOS_UNROLL _Pragma("unroll")
#else
#define RODIN_KOKKOS_UNROLL
#endif

namespace Rodin::Kokkos::Internal
{
  /**
   * @brief Runs a host construction loop, on OpenMP threads when enabled.
   * @tparam Body Callable taking the iteration index.
   * @param concurrent Whether @p body may run concurrently on distinct indices.
   * @param count Iteration count.
   * @param body Loop body; it must not throw when @p concurrent is set.
   */
  template <class Body>
  void forEachOnHost(bool concurrent, size_t count, Body&& body)
  {
#ifdef RODIN_USE_OPENMP
    // Below this count the thread team costs more than the loop.
    constexpr size_t grain = 4096;
    if (concurrent && count >= grain)
    {
      const auto n = static_cast<std::int64_t>(count);
#pragma omp parallel for schedule(static)
      for (std::int64_t i = 0; i < n; ++i)
        body(static_cast<size_t>(i));
      return;
    }
#endif
    for (size_t i = 0; i < count; ++i)
      body(i);
  }

  /**
   * @brief Transfers construction data written in place into an execution space's memory.
   * @tparam T Trivially copyable entry type.
   * @tparam Execution Execution-space type.
   * @tparam Fill Callable writing the entries through a host pointer.
   * @param execution Execution instance.
   * @param count Entry count.
   * @param label Allocation label.
   * @param fill Writes the @p count entries; it may throw.
   * @returns Independently owned resident array, after transfer completion.
   */
  template <class T, class Execution, class Fill>
  auto upload(const Execution& execution, size_t count, const char* label, Fill&& fill)
  {
    static_assert(std::is_trivially_copyable_v<T>);
    ::Kokkos::View<T*, typename Execution::memory_space> result(
      ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing, std::string(label)),
      count);
    // On host memory spaces the mirror is the result itself.
    auto host = ::Kokkos::create_mirror_view(::Kokkos::WithoutInitializing, result);
    fill(host.data());
    ::Kokkos::deep_copy(execution, result, host);
    execution.fence("Rodin construction transfer");
    return result;
  }

  /**
   * @brief Transfers construction data into an execution space's memory.
   * @tparam Execution Execution-space type.
   * @tparam T Trivially copyable entry type.
   * @param execution Execution instance.
   * @param values Temporary host construction values.
   * @param label Allocation label.
   * @returns Independently owned resident array, after transfer completion.
   */
  template <class Execution, class T>
  auto upload(const Execution& execution, const std::vector<T>& values, const char* label)
  {
    // Copying on several threads also spreads the first touch of the staging pages.
    return upload<T>(execution, values.size(), label, [&](T* out) {
      forEachOnHost(true, values.size(), [&](size_t i) { out[i] = values[i]; });
    });
  }

  /**
   * @brief Copies a resident array without aliasing its allocation.
   * @tparam Execution Execution-space type.
   * @tparam View Owning view type.
   * @param execution Execution instance.
   * @param source Resident source.
   * @returns Independent resident copy, after transfer completion.
   */
  template <class Execution, class View>
  View clone(const Execution& execution, const View& source)
  {
    View result(::Kokkos::view_alloc(execution, source.label()), source.extent(0));
    ::Kokkos::deep_copy(execution, result, source);
    execution.fence("Rodin resident copy");
    return result;
  }

  /**
   * @brief Reads one resident entry for a host API query.
   * @tparam Execution Execution-space type.
   * @tparam View Resident array type.
   * @param execution Execution instance.
   * @param source Resident source array.
   * @param index Entry index.
   * @returns Host value after execution completion.
   */
  template <class Execution, class View>
  auto read(const Execution& execution, const View& source, size_t index)
  {
    typename View::non_const_value_type value;
    ::Kokkos::deep_copy(execution, value, ::Kokkos::subview(source, index));
    execution.fence("Rodin host query");
    return value;
  }

  /**
   * @brief Clears the pending launch error of the device runtime, so the
   * next checkLaunches() reports only the kernels launched in between.
   */
  inline void clearLaunches()
  {
#ifdef KOKKOS_ENABLE_CUDA
    static_cast<void>(cudaGetLastError());
#endif
  }

  /**
   * @brief Raises the launch error of a kernel launched since
   * clearLaunches().
   *
   * Kokkos checks the launches only with its bounds checks enabled. A kernel
   * refused at launch, for instance for lack of memory for its local arrays
   * on a nearly full device, otherwise leaves its outputs unwritten without
   * an error, and the assembled system silently misses its entries.
   * @param what Name of the launched work, for the message.
   */
  inline void checkLaunches(const char* what)
  {
#ifdef KOKKOS_ENABLE_CUDA
    const cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
      Alert::Exception() << what << " failed to launch: " << cudaGetErrorString(error)
                         << Alert::Raise;
    }
#else
    static_cast<void>(what);
#endif
  }
}

#endif
