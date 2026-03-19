#ifndef RODIN_PETSC_OBJECT_H
#define RODIN_PETSC_OBJECT_H

/**
 * @file
 * @brief Base helper for wrappers around PETSc opaque handles.
 */

#include <cassert>
#include <petsc.h>

namespace Rodin::PETSc
{
  /**
   * @brief Base‐class for any wrapper around a PETSc object.
   *
   * Each derived instance is automatically tracked in a per‐thread registry,
   * and any remaining PETSc handles are destroyed at PetscFinalize().
   */
  template <class Handle>
  class Object
  {
    public:
      Object() = default;

      virtual ~Object() = default;

      /**
       * @brief Retrieves the communicator associated with the PETSc handle.
       * @param[out] comm Output communicator.
       */
      void getComm(MPI_Comm* comm) const
      {
        PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject&>(this->getHandle()), comm);
        assert(ierr == PETSC_SUCCESS);
      }

      /**
       * @brief Gets read-only access to the wrapped PETSc handle.
       */
      virtual const Handle& getHandle() const noexcept = 0;

      /**
       * @brief Gets mutable access to the wrapped PETSc handle.
       */
      virtual Handle& getHandle() noexcept = 0;
  };
}

#endif
