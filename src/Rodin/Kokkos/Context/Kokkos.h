/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_CONTEXT_KOKKOS_H
#define RODIN_KOKKOS_CONTEXT_KOKKOS_H

#include <type_traits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "Rodin/Context/Local.h"

#include "CUDA.h"

namespace Rodin::Context
{
  namespace Internal
  {
    /**
     * @brief Resolves a backend tag to its Kokkos execution space.
     * @tparam Backend Backend tag or native Kokkos execution-space type.
     */
    template <class Backend>
    struct KokkosBackend
    {
      /// Resolved native execution space.
      using ExecutionSpace = Backend;
    };

#ifdef KOKKOS_ENABLE_CUDA
    template <>
    struct KokkosBackend<CUDA>
    {
      /// CUDA execution space selected by the Rodin tag.
      using ExecutionSpace = ::Kokkos::Cuda;
    };
#endif
  }

  /**
   * @brief Combines a distribution context with a Kokkos execution space.
   * @tparam Distribution Rodin distribution context, initially @ref Local.
   * @tparam Backend @ref CUDA or a native Kokkos execution-space type.
   *
   * ## Architecture
   *
   * The distribution is the base subobject so existing distribution-facing
   * Rodin interfaces remain available. The execution-space instance controls
   * local kernel execution and transfers. Context copies share that execution
   * association through Kokkos execution-space value semantics.
   *
   * The application must initialize Kokkos before constructing this context
   * and finalize it only after the context and its dependent allocations have
   * been destroyed. The context does not initialize or finalize either
   * runtime. CUDA selection requires a CUDA-enabled Kokkos installation;
   * execution never falls back to a host backend implicitly.
   */
  template <class Distribution, class Backend>
  class Kokkos final : public Distribution
  {
    static_assert(std::is_base_of_v<Base, Distribution>,
      "A Kokkos context requires a Rodin distribution context.");

    public:
      /// Distribution and ownership context.
      using DistributionType = Distribution;

      /// Backend tag or native Kokkos execution-space type.
      using BackendType = Backend;

      /// Concrete Kokkos execution-space type.
      using ExecutionSpaceType =
        typename Internal::KokkosBackend<Backend>::ExecutionSpace;

      /// Memory space associated with the execution backend.
      using MemorySpaceType = typename ExecutionSpaceType::memory_space;

      static_assert(::Kokkos::is_execution_space<ExecutionSpaceType>::value,
        "A Kokkos backend must select a Kokkos execution space.");

      /// Constructs the default distribution and execution instances.
      Kokkos() = default;

      /**
       * @brief Constructs a context from distribution and execution instances.
       * @param distribution Distribution association retained by value.
       * @param execution Local execution instance retained by value.
       */
      explicit Kokkos(
        Distribution distribution, ExecutionSpaceType execution = ExecutionSpaceType())
        : Distribution(std::move(distribution)),
          m_execution(std::move(execution))
      {}

      /**
       * @brief Copies a context and its execution association.
       * @param other Context to copy.
       */
      Kokkos(const Kokkos& other) = default;

      /**
       * @brief Moves a context and its execution association.
       * @param other Context to move.
       */
      Kokkos(Kokkos&& other) = default;

      /**
       * @brief Copies a distribution and execution association.
       * @param other Context to copy.
       * @return This context.
       */
      Kokkos& operator=(const Kokkos& other) = default;

      /**
       * @brief Moves a distribution and execution association.
       * @param other Context to move.
       * @return This context.
       */
      Kokkos& operator=(Kokkos&& other) = default;

      /**
       * @brief Gets the distribution context.
       * @return Distribution base subobject.
       */
      const Distribution& getDistribution() const
      {
        return *this;
      }

      /**
       * @brief Gets the local execution instance.
       * @return Kokkos execution instance used by this context.
       */
      const ExecutionSpaceType& getExecutionSpace() const
      {
        return m_execution;
      }

    private:
      ExecutionSpaceType m_execution;
  };

#ifndef KOKKOS_ENABLE_CUDA
  /**
   * @brief Diagnoses CUDA selection with a Kokkos build lacking CUDA support.
   * @tparam Distribution Requested distribution context.
   */
  template <class Distribution>
  class Kokkos<Distribution, CUDA>
  {
    static_assert(!std::is_same_v<Distribution, Distribution>,
      "Context::CUDA requires Kokkos built with CUDA support.");
  };
#endif
}

#endif
