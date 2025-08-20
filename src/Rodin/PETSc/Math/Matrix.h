#ifndef RODIN_PETSC_MATH_MATRIX_H
#define RODIN_PETSC_MATH_MATRIX_H

#include <boost/mpi/communicator.hpp>
#include <mpi.h>
#include <petsc.h>
#include <petscmat.h>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::PETSc::Math
{
  using Matrix = ::Mat;
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<::Mat>
  {
    using ScalarType = PetscScalar;
  };
}

#endif

