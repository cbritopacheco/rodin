#include "Vector.h"

namespace Rodin::Math
{
  void duplicate(const ::Vec& src, ::Vec& dest)
  {
    PetscErrorCode ierr;
    ierr = VecDuplicate(src, &dest);
    assert(ierr == PETSC_SUCCESS);
  }

  void copy(const ::Vec& src, ::Vec& dest)
  {
    PetscErrorCode ierr;
    ierr = VecCopy(src, dest);
    assert(ierr == PETSC_SUCCESS);
  }
}
