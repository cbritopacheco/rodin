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
#include <string>
#include <vector>
#include <Kokkos_Core.hpp>

namespace Rodin::Kokkos::Internal
{
  /**
   * @brief Transfers construction data into an execution space's memory.
   * @tparam Execution Execution-space type.
   * @tparam T Entry type.
   * @param execution Execution instance.
   * @param values Temporary host construction values.
   * @param label Allocation label.
   * @returns Independently owned resident array, after transfer completion.
   */
  template <class Execution, class T>
  auto upload(const Execution& execution, const std::vector<T>& values, const char* label)
  {
    ::Kokkos::View<T*, typename Execution::memory_space> result(
      ::Kokkos::view_alloc(execution, std::string(label)), values.size());
    auto host = ::Kokkos::create_mirror(result);
    std::copy(values.begin(), values.end(), host.data());
    ::Kokkos::deep_copy(execution, result, host);
    execution.fence("Rodin construction transfer");
    return result;
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
}

#endif
