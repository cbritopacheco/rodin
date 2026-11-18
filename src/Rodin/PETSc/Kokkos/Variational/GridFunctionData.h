/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_KOKKOS_VARIATIONAL_GRIDFUNCTIONDATA_H
#define RODIN_PETSC_KOKKOS_VARIATIONAL_GRIDFUNCTIONDATA_H

/**
 * @file
 * @brief State-aware PETSc coefficient access in local device numbering.
 */

#include <petscvec.h>
#include "GhostedPositions.h"
// PETSc declares its object-state accessors in this implementation header.
#include <petsc/private/petscimpl.h>
#include "Rodin/Kokkos/Variational/GridFunction.h"
#include "Rodin/Kokkos/Utility.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Gathers the local form of a ghosted PETSc vector in the local
   * numbering of the device degrees of freedom.
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
   * @brief Copies the coefficients of a PETSc grid function to the device.
   *
   * The values are read where PETSc keeps them, on the host or on the
   * device. On a distributed mesh the ghosted local form is read, without
   * communication: the collective writes of PETSc grid functions keep the
   * ghosts up to date. Its owned entries are numbered from the first owned
   * global index and its ghosts follow in the order of the local degrees of
   * freedom, as PETSc::Variational::GridFunction lays them out.
   */
  template <>
  struct GridFunctionData<::Vec>
  {
    /**
     * @param fes Device space of the grid function.
     * @param data Coefficients of the grid function.
     * @returns Device coefficients in the local numbering, to be read only.
     *
     * The copy is kept in the transfer cache of the mesh under the PETSc
     * states of the vector and of its local form, so that the integrators
     * reading the same grid function share one copy until PETSc writes the
     * vector again.
     */
    template <class FES>
    static ::Kokkos::View<Real*, typename FES::MeshType::MemorySpace> prepare(
      const FES& fes, ::Vec data)
    {
      using MemorySpace = typename FES::MeshType::MemorySpace;
      using Unmanaged = ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>;
      using ContextType = typename FES::MeshType::ContextType;
      using Array = ::Kokkos::View<Real*, MemorySpace>;
      constexpr bool distributed = ::Rodin::Context::IsMPI<ContextType>::value;
      const auto& execution = fes.getMesh().getContext().getExecutionSpace();

      ::Vec form = data;
      PetscErrorCode ierr = PETSC_SUCCESS;
      if constexpr (distributed)
      {
        ierr = VecGhostGetLocalForm(data, &form);
        assert(ierr == PETSC_SUCCESS);
        if (!form)
        {
          Alert::Exception()
            << "The PETSc vector of a distributed grid function has no ghosts."
            << Alert::Raise;
        }
      }
      PetscInt count = 0;
      ierr = VecGetLocalSize(form, &count);
      assert(ierr == PETSC_SUCCESS);
      size_t size = static_cast<size_t>(count);
      if constexpr (distributed)
        size = fes.getLocalSize();

      // Writes through PETSc raise the state of the vector they write,
      // the global vector for its owned entries and the local form for the
      // ghost updates.
      PetscObjectId dataId = 0, formId = 0;
      PetscObjectState dataState = 0, formState = 0;
      ierr = PetscObjectGetId(reinterpret_cast<::PetscObject>(data), &dataId);
      assert(ierr == PETSC_SUCCESS);
      ierr = PetscObjectStateGet(reinterpret_cast<::PetscObject>(data), &dataState);
      assert(ierr == PETSC_SUCCESS);
      ierr = PetscObjectGetId(reinterpret_cast<::PetscObject>(form), &formId);
      assert(ierr == PETSC_SUCCESS);
      ierr = PetscObjectStateGet(reinterpret_cast<::PetscObject>(form), &formState);
      assert(ierr == PETSC_SUCCESS);
      typename ::Rodin::Kokkos::Internal::TransferCache<MemorySpace>::Key key{
        static_cast<std::uint64_t>(dataId), static_cast<std::uint64_t>(dataState),
        static_cast<std::uint64_t>(formId), static_cast<std::uint64_t>(formState),
        static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(&fes))};
      if constexpr (distributed && requires { fes.getRevision(); })
        key.push_back(fes.getRevision());

      const auto fill = [&](const Array& res) {
        const PetscScalar* array = nullptr;
        PetscMemType type;
        PetscErrorCode error = VecGetArrayReadAndMemType(form, &array, &type);
        assert(error == PETSC_SUCCESS);
        // The local form, resident: read in place when PETSc keeps it in the
        // memory space of the mesh, copied otherwise.
        using Device = ::Kokkos::DefaultExecutionSpace::memory_space;
        ::Kokkos::View<const Real*, MemorySpace> values;
        const auto read = [&](const auto& source) {
          if constexpr (!distributed)
          {
            ::Kokkos::deep_copy(execution, res, source);
          }
          else if constexpr (std::is_same_v<
                               typename std::decay_t<decltype(source)>::memory_space,
                               MemorySpace>)
          {
            values = source;
          }
          else
          {
            ::Kokkos::View<Real*, MemorySpace> copy(
              ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
                std::string("Rodin PETSc grid function")),
              count);
            ::Kokkos::deep_copy(execution, copy, source);
            values = copy;
          }
        };
        if (PetscMemTypeHost(type))
        {
          read(::Kokkos::View<const PetscScalar*, ::Kokkos::HostSpace, Unmanaged>(
            array, count));
        }
        else
          read(::Kokkos::View<const PetscScalar*, Device, Unmanaged>(array, count));
        if constexpr (distributed)
        {
          // PETSc owns the cached map independently of coefficient state.
          auto& positions =
            ::Rodin::PETSc::Internal::GhostedPositions<MemorySpace>::of(data);
          positions.setSpace(fes);
          const auto table = positions.getView();
          if (table.size != static_cast<size_t>(count))
          {
            error = VecRestoreArrayReadAndMemType(form, &array);
            assert(error == PETSC_SUCCESS);
            Alert::Exception()
              << "The PETSc vector of a distributed grid function does not have "
              << "the layout of its space." << Alert::Raise;
          }
          ::Kokkos::parallel_for("Rodin PETSc grid function gather",
            ::Kokkos::RangePolicy<std::decay_t<decltype(execution)>>(
              execution, 0, size),
            GhostedCoefficientsGather<MemorySpace>{res, values, table.positions});
        }
        // PETSc may release the array once it is restored.
        execution.fence("Rodin PETSc grid function transfer");
        error = VecRestoreArrayReadAndMemType(form, &array);
        assert(error == PETSC_SUCCESS);
        (void)error;
      };
      const Array res = fes.getMesh().getTransferCache().upload(
        execution, std::move(key), size, "Rodin grid function", fill);

      if constexpr (distributed)
      {
        ierr = VecGhostRestoreLocalForm(data, &form);
        assert(ierr == PETSC_SUCCESS);
      }
      (void)ierr;
      return res;
    }
  };
}

#endif
