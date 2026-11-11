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
#include <cassert>
#include <cstdint>
#include <list>
#include <mutex>
#include <string>
#include <type_traits>
#include <vector>
#include <Kokkos_Core.hpp>
#include "Rodin/Configure.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Types.h"

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
   * @brief Resident copies of host arrays, reused while their content holds.
   *
   * Keeps the arrays last uploaded together with a host snapshot of their
   * content. Uploading the same content again returns the resident copy
   * without any transfer. Arrays whose content lives elsewhere, such as in
   * a PETSc vector, are kept under a key identifying that content instead
   * of a snapshot. Another content of the same size goes to a copy nobody
   * references any more, so that neither an allocation nor a release, both
   * waiting for the device on CUDA, happens between kernel launches. The
   * copies are read only. A copy of the cache starts empty.
   *
   * @tparam MemorySpace Memory space of the copies.
   */
  template <class MemorySpace>
  class TransferCache
  {
    public:
      /// @brief Resident array type.
      using Array = ::Kokkos::View<Real*, MemorySpace>;

      /// @brief Identity of a content kept outside the cache.
      using Key = std::vector<std::uint64_t>;

      /// @brief Largest number of arrays kept.
      static constexpr size_t Capacity = 4;

      TransferCache() = default;

      /// @brief Waits for the transfers still reading the snapshots.
      ~TransferCache()
      {
        if (!m_entries.empty())
          ::Kokkos::fence("Rodin transfer release");
      }

      /**
       * @brief Starts empty, the arrays of the source stay with it.
       * @param other Cache whose arrays are not copied.
       */
      TransferCache([[maybe_unused]] const TransferCache& other) {}

      /**
       * @brief Keeps the arrays of this cache.
       * @param other Cache whose arrays are not copied.
       * @returns This object.
       */
      TransferCache& operator=([[maybe_unused]] const TransferCache& other)
      {
        return *this;
      }

      /**
       * @brief Gets a resident copy of a host array.
       * @param execution Execution instance ordering the transfer.
       * @param values Host values.
       * @param label Label of a new allocation.
       * @returns Resident copy, to be read only.
       */
      template <class Execution>
      Array upload(
        const Execution& execution, std::vector<Real>&& values, const std::string& label)
      {
        std::lock_guard<std::mutex> lock(m_mutex);
        ++m_clock;
        for (auto& entry : m_entries)
        {
          if (entry.key.empty() && entry.snapshot == values)
          {
            entry.used = m_clock;
            return entry.array;
          }
        }
        Entry* target = acquire(execution, values.size(), label);
        if (!target)
        {
          // Every copy is referenced: this one is not kept.
          Array array(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing, label),
            values.size());
          ::Kokkos::deep_copy(execution, array, host(values));
          execution.fence("Rodin transfer");
          return array;
        }
        target->snapshot = std::move(values);
        ::Kokkos::deep_copy(execution, target->array, host(target->snapshot));
        return target->array;
      }

      /**
       * @brief Gets a resident array identified by a key, filled on the
       * device when no array is kept under that key.
       * @param execution Execution instance ordering the fill.
       * @param key Identity of the content, not empty: equal keys must
       * stand for equal contents.
       * @param size Entry count.
       * @param label Label of a new allocation.
       * @param fill Writes the entries of the array it is given, ordered on
       * @p execution.
       * @returns Resident array, to be read only.
       */
      template <class Execution, class Fill>
      Array upload(const Execution& execution, Key key, size_t size,
        const std::string& label, Fill&& fill)
      {
        assert(!key.empty());
        std::lock_guard<std::mutex> lock(m_mutex);
        ++m_clock;
        for (auto& entry : m_entries)
        {
          if (entry.key == key && entry.array.extent(0) == size)
          {
            entry.used = m_clock;
            return entry.array;
          }
        }
        Entry* target = acquire(execution, size, label);
        if (!target)
        {
          // Every copy is referenced: this one is not kept.
          Array array(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing, label), size);
          fill(array);
          return array;
        }
        fill(target->array);
        target->key = std::move(key);
        return target->array;
      }

    private:
      struct Entry
      {
        std::vector<Real> snapshot;
        Key key;
        Array array;
        std::uint64_t used = 0;
      };

      /// @brief Unmanaged host view of host values.
      static auto host(const std::vector<Real>& values)
      {
        return ::Kokkos::View<const Real*, ::Kokkos::HostSpace,
          ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(values.data(), values.size());
      }

      /**
       * @brief Gets the entry receiving a new content, emptied.
       * @param execution Execution instance ordering the transfers.
       * @param size Entry count of the content.
       * @param label Label of a new allocation.
       * @returns Entry holding an array of the size, or null when every
       * copy is referenced.
       *
       * A free copy of the same size is overwritten in place; otherwise the
       * least recently used free copy is released.
       */
      template <class Execution>
      Entry* acquire(const Execution& execution, size_t size, const std::string& label)
      {
        Entry* target = nullptr;
        Entry* oldest = nullptr;
        for (auto& entry : m_entries)
        {
          if (entry.array.use_count() != 1)
            continue;
          if (entry.array.extent(0) == size && (!target || entry.used < target->used))
            target = &entry;
          if (!oldest || entry.used < oldest->used)
            oldest = &entry;
        }
        if (!target)
        {
          if (m_entries.size() < Capacity)
            target = &m_entries.emplace_back();
          else if (oldest)
            target = oldest;
          else
            return nullptr;
          if (target->array.extent(0) != size || !target->array.data())
          {
            target->array = Array();
            target->array = Array(
              ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing, label), size);
          }
        }
        // A former transfer may still read the snapshot being replaced.
        if (!target->snapshot.empty())
          execution.fence("Rodin transfer reuse");
        target->snapshot.clear();
        target->key.clear();
        target->used = m_clock;
        return target;
      }

      std::mutex m_mutex;
      std::list<Entry> m_entries;
      std::uint64_t m_clock = 0;
  };

  /**
   * @brief Clears the pending launch error of the device runtime, so the
   * next checkLaunches() reports only the kernels launched in between.
   */
  inline void clearLaunches()
  {
#if defined(KOKKOS_ENABLE_CUDA)
    static_cast<void>(cudaGetLastError());
#elif defined(KOKKOS_ENABLE_HIP)
    static_cast<void>(hipGetLastError());
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
#if defined(KOKKOS_ENABLE_CUDA)
    const cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
      Alert::Exception() << what << " failed to launch: " << cudaGetErrorString(error)
                         << Alert::Raise;
    }
#elif defined(KOKKOS_ENABLE_HIP)
    const hipError_t error = hipGetLastError();
    if (error != hipSuccess)
    {
      Alert::Exception() << what << " failed to launch: " << hipGetErrorString(error)
                         << Alert::Raise;
    }
#else
    static_cast<void>(what);
#endif
  }
}

#endif
