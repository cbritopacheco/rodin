/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_VARIATIONAL_GHOSTEDPOSITIONS_H
#define RODIN_PETSC_KOKKOS_VARIATIONAL_GHOSTEDPOSITIONS_H

/**
 * @file
 * @brief PETSc-owned local-to-ghosted numbering maps.
 */

#include <memory>
#include <petscvec.h>
#include "Rodin/Kokkos/Utility.h"

namespace Rodin::PETSc::Internal
{
  /**
   * @brief Retains a distributed space's numbering map with its PETSc vector.
   * @tparam MemorySpace Memory space of the cached map.
   *
   * ## Architecture
   *
   * A PETSc container attached to the vector owns one position table per memory
   * space. A space's immutable numbering revision validates reuse, independently
   * of coefficient writes. Copies, moves and reconstruction at the same address
   * have distinct revisions. Replacing a map leaves already returned views valid.
   * Destroying the vector releases the cache; no static device allocations remain.
   * Access is serialized by the caller, as for the PETSc vector itself. Creation
   * of the container uses PETSC_COMM_SELF and adds no inter-rank collective.
   */
  template <class MemorySpace>
  class GhostedPositions
  {
    public:
      /// @brief Positions and size of the local form.
      using Table =
        typename ::Rodin::Kokkos::Internal::GhostedPositions<MemorySpace>::Table;

      /**
       * @brief Gets the adapter state owned by a PETSc vector.
       * @param data Vector owning the adapter cache.
       * @returns Non-const adapter state for explicit space binding.
       */
      static GhostedPositions& of(::Vec data)
      {
        const auto object = reinterpret_cast<PetscObject>(data);
        const std::string key =
          std::string("Rodin ghosted positions ") + MemorySpace::name();
        PetscObject composed = nullptr;
        [[maybe_unused]] PetscErrorCode ierr =
          PetscObjectQuery(object, key.c_str(), &composed);
        assert(ierr == PETSC_SUCCESS);
        if (!composed)
        {
          PetscContainer container = nullptr;
          auto owner = std::make_unique<GhostedPositions>();
          ierr = PetscContainerCreate(PETSC_COMM_SELF, &container);
          assert(ierr == PETSC_SUCCESS);
          ierr = PetscContainerSetPointer(container, owner.get());
          assert(ierr == PETSC_SUCCESS);
#if PETSC_VERSION_GE(3, 23, 0)
          ierr = PetscContainerSetCtxDestroy(container, destroy);
#else
          ierr = PetscContainerSetUserDestroy(container, destroy);
#endif
          assert(ierr == PETSC_SUCCESS);
          composed = reinterpret_cast<PetscObject>(container);
          ierr = PetscObjectCompose(object, key.c_str(), composed);
          assert(ierr == PETSC_SUCCESS);
          owner.release();
          ierr = PetscContainerDestroy(&container);
          assert(ierr == PETSC_SUCCESS);
        }
        void* pointer = nullptr;
        ierr =
          PetscContainerGetPointer(reinterpret_cast<PetscContainer>(composed), &pointer);
        assert(ierr == PETSC_SUCCESS);
        return *static_cast<GhostedPositions*>(pointer);
      }

      /**
       * @brief Binds the immutable numbering and prepares its device map.
       * @tparam FES Distributed finite-element space.
       * @param fes Space supplying owned and ghost numbering.
       * @returns This adapter.
       *
       * Rebinding the same numbering reuses its table. Spaces without a
       * numbering revision are prepared anew, preserving the query interface.
       */
      template <class FES>
      GhostedPositions& setSpace(const FES& fes)
      {
        if constexpr (requires { fes.getRevision(); })
        {
          const auto revision = fes.getRevision();
          if (!m_entry || m_entry->revision != revision)
            m_entry = Entry{revision, fes.getGhostedPositions()};
        }
        else
          m_entry = Entry{Optional<std::uint64_t>{}, fes.getGhostedPositions()};
        return *this;
      }

      /**
       * @brief Reads the prepared map without allocation or state changes.
       * @pre setSpace() has completed.
       * @returns Read-only positions retaining their allocation independently.
       */
      Table getView() const
      {
        assert(m_entry);
        return m_entry->table;
      }

    private:
      struct Entry
      {
        /// @brief Revision of the numbering that produced the table.
        Optional<std::uint64_t> revision;

        /// @brief Retained device position table.
        Table table;
      };

#if PETSC_VERSION_GE(3, 23, 0)
      /**
       * @tparam Pointer Context argument of this PETSc version's destroy callback.
       * @param pointer Address of the owned cache context, cleared after deletion.
       * @returns PETSc success.
       */
      template <class Pointer>
      static PetscErrorCode destroy(Pointer pointer)
      {
        auto** context = static_cast<void**>(pointer);
        delete static_cast<GhostedPositions*>(*context);
        *context = nullptr;
        return PETSC_SUCCESS;
      }
#else
      /// @param pointer Owned cache context. @returns PETSc success.
      static PetscErrorCode destroy(void* pointer)
      {
        delete static_cast<GhostedPositions*>(pointer);
        return PETSC_SUCCESS;
      }
#endif

      Optional<Entry> m_entry;
  };
}

#endif
