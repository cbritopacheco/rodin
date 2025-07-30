#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

#include <boost/mpi/communicator.hpp>

#include <mpi.h>

#include <petsc.h>
#include <petscsys.h>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"

namespace Rodin::PETSc::Math
{
  using Vector = ::Vec;
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<::Vec>
  {
    using ScalarType = PetscScalar;
  };
}

#endif
