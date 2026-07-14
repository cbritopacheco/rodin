/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_OBJECT_H
#define RODIN_PETSC_OBJECT_H

/**
 * @file
 * @brief RAII base class for PETSc opaque handle wrappers.
 *
 * Provides a common interface for accessing the underlying PETSc handle
 * and querying the associated MPI communicator.
 *
 * @see Rodin::Solver::KSP, Rodin::Solver::SNES
 */

#include <cassert>
#include <petsc.h>

namespace Rodin::PETSc
{
  /**
   * @brief RAII base class for wrappers around PETSc opaque handles.
   *
   * Derived classes such as the KSP and SNES solver wrappers must implement
   * getHandle() to expose
   * their PETSc object, which is used by getComm() to retrieve the
   * associated communicator.
   *
   * @tparam Handle PETSc opaque handle type.
   */
  template <class Handle>
  class Object
  {
    public:
      Object() = default;

      virtual ~Object() = default;

      /**
       * @brief Retrieves the MPI communicator associated with the handle.
       * @param[out] comm Pointer to the output communicator.
       */
      void getComm(MPI_Comm* comm) const
      {
        PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject&>(this->getHandle()), comm);
        assert(ierr == PETSC_SUCCESS);
      }

      /**
       * @brief Returns a read-only reference to the underlying PETSc handle.
       */
      virtual const Handle& getHandle() const noexcept = 0;

      /**
       * @brief Returns a mutable reference to the underlying PETSc handle.
       */
      virtual Handle& getHandle() noexcept = 0;
  };
}

#endif
