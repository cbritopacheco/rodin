/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_VARIATIONAL_GRIDFUNCTIONADAPTER_H
#define RODIN_PETSC_KOKKOS_VARIATIONAL_GRIDFUNCTIONADAPTER_H

/**
 * @file
 * @brief PETSc coefficient reads and distributed numbering.
 */

#include <memory>
#include <string>
#include <petscvec.h>
// PETSc declares its object-state accessors in this implementation header.
#include <petsc/private/petscimpl.h>
#if defined(PETSC_HAVE_KOKKOS) && __has_include(<petsc_kokkos.hpp>)
#include <petsc_kokkos.hpp>
#endif
#include "Rodin/Kokkos/Variational/GridFunction.h"
#include "Rodin/Kokkos/Utility.h"

namespace Rodin::PETSc::Internal
{
  /**
   * @brief Gathers the local form of a ghosted PETSc vector in the local
   * numbering of the device degrees of freedom.
   * @tparam MemorySpace Memory space of the input and output views.
   *
   * Kernel functor at namespace scope, as NVCC requires.
   */
  template <class MemorySpace>
  struct GhostedCoefficientsGather
  {
    /// @brief Coefficients in the local numbering.
    ::Kokkos::View<Real*, MemorySpace> res;

    /// @brief Local form of the vector, owned entries then ghosts.
    ::Kokkos::View<const Real*, MemorySpace> form;

    /// @brief Position in the local form of each local degree of freedom.
    ::Kokkos::View<const Index*, MemorySpace> position;

    /// @param i Local degree of freedom.
    KOKKOS_INLINE_FUNCTION void operator()(Index i) const
    {
      res(i) = form(position(i));
    }
  };

  /**
   * @brief Retains the local-to-ghosted numbering of a PETSc vector.
   * @tparam Execution Kokkos execution space preparing the coefficients.
   *
   * ## Architecture
   * A PETSc container retains one adapter per vector and execution-space type.
   * Space binding prepares the numbering independently of vector writes.
   * The immutable numbering revision validates reuse.
   * Access is serialized by the caller and adds no inter-rank communication.
   */
  template <class Execution>
  class GridFunctionAdapter
  {
    public:
      /// @brief Memory space of the prepared coefficients and positions.
      using MemorySpace = typename Execution::memory_space;

      /// @brief Positions and size of the local form.
      using Table =
        typename ::Rodin::Kokkos::Internal::GhostedPositions<MemorySpace>::Table;

      /// @param owner Identity of the vector that owns the adapter.
      explicit GridFunctionAdapter(PetscObjectId owner)
        : m_owner(owner)
      {}

      /**
       * @brief Gets the adapter state owned by a PETSc vector.
       * @param data Vector owning the adapter.
       * @returns Non-const adapter state for preparation.
       */
      static GridFunctionAdapter& of(::Vec data)
      {
        const auto object = reinterpret_cast<PetscObject>(data);
        static const std::string key =
          std::string("Rodin ghosted numbering ") + Execution::name();
        PetscObjectId id = 0;
        [[maybe_unused]] PetscErrorCode ierr = PetscObjectGetId(object, &id);
        assert(ierr == PETSC_SUCCESS);
        PetscObject composed = nullptr;
        ierr = PetscObjectQuery(object, key.c_str(), &composed);
        assert(ierr == PETSC_SUCCESS);
        if (composed)
        {
          void* pointer = nullptr;
          ierr = PetscContainerGetPointer(
            reinterpret_cast<PetscContainer>(composed), &pointer);
          assert(ierr == PETSC_SUCCESS);
          auto& adapter = *static_cast<GridFunctionAdapter*>(pointer);
          // VecDuplicate can inherit the source vector's composed container.
          if (adapter.m_owner == id)
            return adapter;
        }

        PetscContainer container = nullptr;
        auto owner = std::make_unique<GridFunctionAdapter>(id);
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
        auto& adapter = *owner.release();
        ierr = PetscContainerDestroy(&container);
        assert(ierr == PETSC_SUCCESS);
        return adapter;
      }

      /**
       * @brief Binds the execution instance and prepares the ghosted numbering if needed.
       * @tparam FES Finite element space.
       * @param fes Space supplying execution and numbering.
       * @returns This adapter.
       */
      template <class FES>
      GridFunctionAdapter& setSpace(const FES& fes)
      {
        const auto& execution = fes.getMesh().getContext().getExecutionSpace();
        if (!m_execution || *m_execution != execution)
        {
          m_positions.reset();
          m_execution = execution;
        }
        if constexpr (::Rodin::Context::IsMPI<typename FES::MeshType::ContextType>::value)
        {
          Optional<std::uint64_t> revision;
          if constexpr (requires { fes.getRevision(); })
            revision = static_cast<std::uint64_t>(fes.getRevision());
          if (!revision || !m_positions || m_positions->revision != revision)
          {
            m_positions = PositionEntry{revision, fes.getGhostedPositions()};
          }
        }
        else if (m_positions)
        {
          m_positions.reset();
        }
        return *this;
      }

      /**
       * @brief Reads the prepared map of a distributed space.
       * @pre setSpace() has completed for the distributed space.
       * @returns Read-only positions retaining their allocation independently.
       */
      Table getGhostedPositions() const
      {
        assert(m_positions);
        return m_positions->table;
      }

      /// @returns Immutable numbering revision, if the space supplies one.
      Optional<std::uint64_t> getRevision() const
      {
        assert(m_positions);
        return m_positions->revision;
      }

    private:
      struct PositionEntry
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
        delete static_cast<GridFunctionAdapter*>(*context);
        *context = nullptr;
        return PETSC_SUCCESS;
      }
#else
      /// @param pointer Owned cache context. @returns PETSc success.
      static PetscErrorCode destroy(void* pointer)
      {
        delete static_cast<GridFunctionAdapter*>(pointer);
        return PETSC_SUCCESS;
      }
#endif

      /// @brief Owner of this adapter, even when inherited by VecDuplicate.
      PetscObjectId m_owner;

      /// @brief Bound execution instance ordering reads and refreshes.
      Optional<Execution> m_execution;

      /// @brief Numbering for the bound distributed space.
      Optional<PositionEntry> m_positions;
  };

  /**
   * @brief Reads a PETSc vector or its ghosted local form.
   * @tparam Execution Kokkos execution space consuming the values.
   *
   * Array access starts lazily after cache validation. Borrowed arrays remain
   * valid until destruction. The caller completes all reads before destruction.
   */
  template <class Execution>
  class GridFunctionRead
  {
    public:
      /// @brief Memory space of the consumer.
      using MemorySpace = typename Execution::memory_space;

      /// @brief Owned coefficients in Rodin numbering.
      using Array = ::Kokkos::View<Real*, MemorySpace>;

      /// @brief Read-only coefficients in Rodin numbering.
      using View = ::Kokkos::View<const Real*, MemorySpace>;

      /// @brief Map from Rodin numbering to the ghosted local form.
      using Table = typename GridFunctionAdapter<Execution>::Table;

      /// @brief Identity and state of both PETSc vector handles.
      struct State
      {
        PetscObjectId vector = 0; ///< Owning vector identity.
        PetscObjectId form = 0; ///< Local form identity.
        PetscObjectState values = 0; ///< Owning vector state.
        PetscObjectState ghosts = 0; ///< Local form state.

        /// @param other State to compare. @returns Whether the states match.
        bool operator==(const State& other) const = default;
      };

      /**
       * @param execution Execution instance consuming the data.
       * @param data Coefficient vector.
       * @param positions Numbering map for distributed vectors.
       */
      GridFunctionRead(
        const Execution& execution, ::Vec data, Optional<Table> positions = {})
        : m_execution(execution),
          m_data(data),
          m_form(data),
          m_positions(std::move(positions))
      {
        [[maybe_unused]] PetscErrorCode ierr;
        if (m_positions)
        {
          ierr = VecGhostGetLocalForm(data, &m_form);
          assert(ierr == PETSC_SUCCESS);
          if (!m_form)
          {
            Alert::Exception()
              << "The PETSc vector of a distributed grid function has no ghosts."
              << Alert::Raise;
          }
        }
        ierr = VecGetLocalSize(m_form, &m_count);
        assert(ierr == PETSC_SUCCESS);
        if (m_positions && m_positions->size != static_cast<size_t>(m_count))
        {
          ierr = VecGhostRestoreLocalForm(data, &m_form);
          assert(ierr == PETSC_SUCCESS);
          Alert::Exception()
            << "The PETSc vector of a distributed grid function does not have "
               "the layout of its space."
            << Alert::Raise;
        }
      }

      /// @param other Read access owned by another guard.
      GridFunctionRead(const GridFunctionRead& other) = delete;

      /**
       * @param other Read access owned by another guard.
       * @returns This guard.
       */
      GridFunctionRead& operator=(const GridFunctionRead& other) = delete;

      /// @brief Ends array and local-form access in reverse order.
      ~GridFunctionRead()
      {
        [[maybe_unused]] PetscErrorCode ierr;
        if (m_access == Access::MemoryType)
        {
          ierr = VecRestoreArrayReadAndMemType(m_form, &m_array);
          assert(ierr == PETSC_SUCCESS);
        }
        else if (m_access == Access::Host)
        {
          ierr = VecRestoreArrayRead(m_form, &m_array);
          assert(ierr == PETSC_SUCCESS);
        }
        if (m_positions)
        {
          ierr = VecGhostRestoreLocalForm(m_data, &m_form);
          assert(ierr == PETSC_SUCCESS);
        }
      }

      /// @returns Coefficient count in Rodin numbering.
      size_t getSize() const
      {
        return m_positions ? m_positions->positions.extent(0)
                           : static_cast<size_t>(m_count);
      }

      /// @returns Current vector and local-form identities and states.
      State getState() const
      {
        State state;
        const auto data = reinterpret_cast<PetscObject>(m_data);
        const auto form = reinterpret_cast<PetscObject>(m_form);
        [[maybe_unused]] PetscErrorCode ierr = PetscObjectGetId(data, &state.vector);
        assert(ierr == PETSC_SUCCESS);
        ierr = PetscObjectGetId(form, &state.form);
        assert(ierr == PETSC_SUCCESS);
        ierr = PetscObjectStateGet(data, &state.values);
        assert(ierr == PETSC_SUCCESS);
        ierr = PetscObjectStateGet(form, &state.ghosts);
        assert(ierr == PETSC_SUCCESS);
        return state;
      }

      /// @returns A direct local view when the execution space can read the array.
      Optional<View> borrow()
      {
        assert(!m_positions);
        acquire();
        if (m_accessible)
          return View(m_array, m_count);
        return {};
      }

      /// @param result Output in Rodin numbering. The copy completes before returning.
      void copy(const Array& result)
      {
        acquire();
        using Unmanaged = ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>;
        using Device = ::Kokkos::DefaultExecutionSpace::memory_space;
        View values;
        const auto read = [&](const auto& source) {
          if (!m_positions)
            ::Kokkos::deep_copy(m_execution, result, source);
          else if constexpr (std::is_same_v<
                               typename std::decay_t<decltype(source)>::memory_space,
                               MemorySpace>)
            values = source;
          else
          {
            Array staging(::Kokkos::view_alloc(m_execution, ::Kokkos::WithoutInitializing,
                            "Rodin PETSc local form"),
              m_count);
            ::Kokkos::deep_copy(m_execution, staging, source);
            values = staging;
          }
        };
        if (PetscMemTypeHost(m_type))
        {
          read(::Kokkos::View<const PetscScalar*, ::Kokkos::HostSpace, Unmanaged>(
            m_array, m_count));
        }
        else
          read(::Kokkos::View<const PetscScalar*, Device, Unmanaged>(m_array, m_count));
        if (m_positions)
        {
          ::Kokkos::parallel_for("Rodin PETSc grid function gather",
            ::Kokkos::RangePolicy<Execution, ::Kokkos::IndexType<Index>>(
              m_execution, 0, getSize()),
            GhostedCoefficientsGather<MemorySpace>{
              result, values, m_positions->positions});
        }
        m_execution.fence("Rodin PETSc grid function transfer");
      }

    private:
      /// @brief PETSc access to restore on destruction.
      enum class Access
      {
        None,
        MemoryType,
        Host
      };

      /// @brief Acquires the freshest array and checks whether it can be used in place.
      void acquire()
      {
        if (m_access != Access::None)
          return;
        [[maybe_unused]] PetscErrorCode ierr =
          VecGetArrayReadAndMemType(m_form, &m_array, &m_type);
        assert(ierr == PETSC_SUCCESS);
        m_access = Access::MemoryType;
        if (PetscMemTypeHost(m_type))
        {
          m_accessible = std::is_same_v<MemorySpace, ::Kokkos::HostSpace>;
          return;
        }
#if defined(PETSC_HAVE_KOKKOS) && __has_include(<petsc_kokkos.hpp>)
        if constexpr (std::is_same_v<Execution, ::Kokkos::DefaultExecutionSpace>)
        {
          PetscBool kokkos = PETSC_FALSE;
          ierr = PetscObjectTypeCompare(
            reinterpret_cast<PetscObject>(m_form), VECSEQKOKKOS, &kokkos);
          assert(ierr == PETSC_SUCCESS);
          if (kokkos)
          {
            const Execution producer = PetscGetKokkosExecutionSpace();
            if constexpr (requires { producer.cuda_device(); })
              m_accessible = producer.cuda_device() == m_execution.cuda_device();
            else if constexpr (requires { producer.hip_device(); })
              m_accessible = producer.hip_device() == m_execution.hip_device();
            else if constexpr (requires { producer.sycl_queue(); })
            {
              m_accessible = producer.sycl_queue().get_context() ==
                  m_execution.sycl_queue().get_context() &&
                producer.sycl_queue().get_device() ==
                  m_execution.sycl_queue().get_device();
            }
            if (m_accessible)
            {
              if (producer != m_execution)
                producer.fence("Rodin PETSc coefficient producer");
              return;
            }
          }
        }
#endif
        // An unfamiliar device pointer is staged through PETSc's host access.
        ierr = VecRestoreArrayReadAndMemType(m_form, &m_array);
        assert(ierr == PETSC_SUCCESS);
        m_access = Access::None;
        ierr = VecGetArrayRead(m_form, &m_array);
        assert(ierr == PETSC_SUCCESS);
        m_access = Access::Host;
        m_type = PETSC_MEMTYPE_HOST;
      }

      const Execution& m_execution;
      ::Vec m_data;
      ::Vec m_form;
      Optional<Table> m_positions;
      PetscInt m_count = 0;
      const PetscScalar* m_array = nullptr;
      PetscMemType m_type = PETSC_MEMTYPE_HOST;
      Access m_access = Access::None;
      bool m_accessible = false;
  };
}

#endif
