#ifndef RODIN_PETSC_TRAITS_H
#define RODIN_PETSC_TRAITS_H

#include <petsc.h>

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/Math/Vector.h"

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<::Mat>
  {
    using ScalarType = PetscScalar;
  };

  template <>
  struct Traits<::Vec>
  {
    using ScalarType = PetscScalar;
  };
}

namespace Rodin::PETSc
{}

#endif

