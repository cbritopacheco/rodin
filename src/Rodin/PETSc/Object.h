#ifndef RODIN_PETSC_OBJECT_H
#define RODIN_PETSC_OBJECT_H

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

      void getComm(MPI_Comm* comm) const
      {
        PetscErrorCode ierr = PetscObjectGetComm(reinterpret_cast<PetscObject&>(this->getHandle()), comm);
        assert(ierr == PETSC_SUCCESS);
      }

      virtual const Handle& getHandle() const noexcept = 0;

      virtual Handle& getHandle() noexcept = 0;
  };
}

#endif
