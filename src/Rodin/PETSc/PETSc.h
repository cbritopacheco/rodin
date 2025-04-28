#ifndef RODIN_PETSC_PETSC_H
#define RODIN_PETSC_PETSC_H

#include <petsc.h>
#include "Rodin/Math/Vector.h"

namespace Rodin::PETSc
{
  void copyVector(const Math::Vector<PetscScalar>& source, ::Vec& destination);
}

#endif
