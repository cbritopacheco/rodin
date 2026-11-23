/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_VARIATIONAL_GRIDFUNCTIONDATA_H
#define RODIN_PETSC_KOKKOS_VARIATIONAL_GRIDFUNCTIONDATA_H

/**
 * @file
 * @brief Assembly-scoped PETSc coefficient views and bounded gather reuse.
 */

#include <list>
#include <utility>
#include <vector>
#include "GridFunctionAdapter.h"

namespace Rodin::PETSc::Internal
{
  /**
   * @brief PETSc assembler storage for reusable distributed coefficient gathers.
   * @tparam Execution Kokkos execution space of the assembler.
   *
   * ## Architecture
   * An assembler owns the idle gather allocations, in most-recently-used order.
   * A scope owns borrows and temporarily prepared arrays until the kernels finish.
   * The thread-local pointer connects host-side expression lowering to that scope;
   * its storage, including all retained device allocations, belongs to the assembler.
   */
  template <class Execution>
  class CoefficientCache
  {
    using Read = GridFunctionRead<Execution>;
    using Array = typename Read::Array;
    using View = typename Read::View;

    struct Key
    {
      typename Read::State state; ///< PETSc vector and local-form state.
      std::uint64_t revision; ///< Rodin numbering revision.
      size_t size; ///< Rodin coefficient count.

      /// @param other Key to compare. @returns Whether the snapshots match.
      bool operator==(const Key& other) const = default;
    };

    struct Entry
    {
      Key key; ///< State of the retained values.
      Array values; ///< Owned gather.
    };

    public:
      /// @brief Default per-assembler retention allowance in bytes.
      static constexpr size_t DefaultBudget = 32 * 1024 * 1024;

      /// @brief Uses the default retention allowance.
      CoefficientCache() = default;

      /// @param budget Maximum retained coefficient bytes, excluding active temporary buffers.
      explicit CoefficientCache(size_t budget)
        : m_budget(budget)
      {}

      /// @returns Bytes held for reuse between assembly calls.
      size_t getRetainedBytes() const
      {
        size_t bytes = 0;
        for (const auto& entry : m_entries)
          bytes += entry.key.size * sizeof(Real);
        return bytes;
      }

      /**
       * @brief Owns coefficient accesses for a single assembly.
       *
       * Expressions copy lightweight read-only views. Scope destruction restores
       * PETSc arrays and discards temporary buffers after execution completes.
       */
      class Scope
      {
        public:
          /**
           * @param cache Assembler retention storage.
           * @param execution Execution instance running the assembly kernels.
           */
          Scope(CoefficientCache& cache, const Execution& execution)
            : m_cache(cache),
              m_execution(execution),
              m_previous(s_current)
          {
            if (!cache.m_execution || *cache.m_execution != execution)
            {
              cache.m_entries.clear();
              cache.m_execution = execution;
            }
            s_current = this;
          }

          Scope(const Scope&) = delete;
          Scope& operator=(const Scope&) = delete;

          /// @brief Restores PETSc access on both normal and exceptional exits.
          ~Scope()
          {
            if (s_current == this)
            {
              if (!m_local.empty() || !m_gathers.empty())
                m_execution.fence("Rodin PETSc coefficient scope");
              release();
            }
          }

          /**
           * @param execution Execution instance of the coefficient's mesh.
           * @returns Active scope for that instance, or null for standalone lowering.
           */
          static Scope* current(const Execution& execution)
          {
            return s_current && s_current->m_execution == execution ? s_current : nullptr;
          }

          /**
           * @param data Local PETSc vector with Rodin coefficient ordering.
           * @returns A view valid until the scope is released.
           */
          View prepare(::Vec data)
          {
            for (const auto& entry : m_local)
            {
              if (entry.vector == data)
                return entry.values;
            }
            auto read = std::make_unique<Read>(m_execution, data);
            if (const auto borrowed = read->borrow())
            {
              m_local.push_back({data, *borrowed, std::move(read)});
            }
            else
            {
              Array staged(
                ::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
                  "Rodin PETSc staged coefficients"),
                read->getSize());
              read->copy(staged);
              m_local.push_back({data, staged, nullptr});
            }
            return m_local.back().values;
          }

          /**
           * @param read Access to the ghosted local form.
           * @param revision Immutable Rodin numbering revision, when available.
           * @returns A gathered view valid until the scope is released.
           */
          View prepare(Read& read, Optional<std::uint64_t> revision)
          {
            Optional<Key> key;
            if (revision)
            {
              key = Key{read.getState(), *revision, read.getSize()};
              for (const auto& entry : m_gathers)
              {
                if (entry.key == key)
                  return entry.values;
              }
            }
            m_gathers.push_back({key, m_cache.prepare(m_execution, read, key)});
            return m_gathers.back().values;
          }

          /**
           * @brief Restores array access and resumes the enclosing scope.
           * @pre All coefficient readers on the execution instance have completed.
           */
          void release()
          {
            assert(s_current == this);
            s_current = m_previous;
            m_local.clear();
            m_gathers.clear();
          }

        private:
          struct Local
          {
            ::Vec vector; ///< Vector consumed by this assembly.
            View values; ///< Borrowed or staged values.
            std::unique_ptr<Read> access; ///< Read interval for a borrowed array.
          };

          struct Gather
          {
            Optional<Key> key; ///< Reusable state, when the numbering has a revision.
            View values; ///< Gather owned until the assembly finishes.
          };

          inline static thread_local Scope* s_current = nullptr;
          CoefficientCache& m_cache;
          const Execution& m_execution;
          Scope* m_previous;
          std::vector<Local> m_local;
          std::vector<Gather> m_gathers;
      };

    private:
      /**
       * @param execution Execution instance consuming the gather.
       * @param read Local-form access for a cache miss.
       * @param key State of the requested gather, when reuse is valid.
       * @returns An owned, prepared gather.
       */
      Array prepare(const Execution& execution, Read& read, const Optional<Key>& key)
      {
        Array result;
        const size_t size = read.getSize();
        const bool eligible = key && size > 0 && size <= m_budget / sizeof(Real);
        if (eligible)
        {
          for (auto it = m_entries.begin(); it != m_entries.end(); ++it)
          {
            if (it->key == *key)
            {
              m_entries.splice(m_entries.begin(), m_entries, it);
              return m_entries.front().values;
            }
            const auto& old = it->key;
            if (old.state.vector == key->state.vector &&
              old.state.form == key->state.form && old.revision == key->revision &&
              old.size == size && it->values.use_count() == 1)
            {
              result = std::move(it->values);
              m_entries.erase(it);
              break;
            }
          }
        }
        bool retain = eligible;
        if (retain)
        {
          const size_t bytes = size * sizeof(Real);
          size_t retained = getRetainedBytes();
          for (auto it = m_entries.end();
               retained > m_budget - bytes && it != m_entries.begin();)
          {
            --it;
            if (it->values.use_count() == 1)
            {
              retained -= it->key.size * sizeof(Real);
              if (!result.data() && it->key.size == size)
                result = std::move(it->values);
              it = m_entries.erase(it);
            }
          }
          retain = retained <= m_budget - bytes;
        }
        if (!result.data())
        {
          result = Array(::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                           "Rodin PETSc gathered coefficients"),
            size);
        }
        read.copy(result);
        if (retain)
          m_entries.push_front({*key, result});
        return result;
      }

      size_t m_budget = DefaultBudget;
      Optional<Execution> m_execution;
      std::list<Entry> m_entries;
  };
}

namespace Rodin::Kokkos::Variational
{
  /// @brief Prepares PETSc coefficients in the local numbering of a Kokkos space.
  template <>
  struct GridFunctionData<::Vec>
  {
    /**
     * @tparam FES Finite element space.
     * @param fes Device space of the grid function.
     * @param data Coefficient vector.
     * @returns Independent device coefficients for standalone evaluation.
     */
    template <class FES>
    static ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> prepare(
      const FES& fes, ::Vec data)
    {
      using Execution = typename FES::MeshType::ContextType::ExecutionSpaceType;
      using Read = ::Rodin::PETSc::Internal::GridFunctionRead<Execution>;
      using Adapter = ::Rodin::PETSc::Internal::GridFunctionAdapter<Execution>;
      const auto& execution = fes.getMesh().getContext().getExecutionSpace();
      Optional<typename Read::Table> table;
      if constexpr (::Rodin::Context::IsMPI<typename FES::MeshType::ContextType>::value)
        table = Adapter::of(data).setSpace(fes).getGhostedPositions();
      Read read(execution, data, std::move(table));
      typename Read::Array values(
        ::Kokkos::view_alloc(
          execution, ::Kokkos::WithoutInitializing, "Rodin PETSc coefficients"),
        read.getSize());
      read.copy(values);
      return values;
    }

    /**
     * @tparam FES Finite element space.
     * @param fes Device space of the grid function.
     * @param data Coefficient vector.
     * @returns Scoped read-only coefficients, or an owned copy for standalone lowering.
     */
    template <class FES>
    static ::Kokkos::View<const Real*, typename FES::MeshType::MemorySpace> prepareView(
      const FES& fes, ::Vec data)
    {
      using Execution = typename FES::MeshType::ContextType::ExecutionSpaceType;
      using Scope =
        typename ::Rodin::PETSc::Internal::CoefficientCache<Execution>::Scope;
      const auto& execution = fes.getMesh().getContext().getExecutionSpace();
      if (auto* scope = Scope::current(execution))
      {
        if constexpr (::Rodin::Context::IsMPI<
                        typename FES::MeshType::ContextType>::value)
        {
          using Read = ::Rodin::PETSc::Internal::GridFunctionRead<Execution>;
          using Adapter = ::Rodin::PETSc::Internal::GridFunctionAdapter<Execution>;
          const auto& adapter = Adapter::of(data).setSpace(fes);
          Read read(execution, data, adapter.getGhostedPositions());
          return scope->prepare(read, adapter.getRevision());
        }
        else
          return scope->prepare(data);
      }
      return prepare(fes, data);
    }
  };
}

#endif
