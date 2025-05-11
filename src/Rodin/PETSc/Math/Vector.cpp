#include "Vector.h"

namespace Rodin::Math
{
  void duplicate(const ::Vec& src, ::Vec& dst)
  {
    PetscErrorCode ierr;
    ierr = VecDuplicate(src, &dst);
    assert(ierr == PETSC_SUCCESS);
  }

  void copy(const ::Vec& src, ::Vec& dst)
  {
    PetscErrorCode ierr;
    ierr = VecCopy(src, dst);
    assert(ierr == PETSC_SUCCESS);
  }
}
