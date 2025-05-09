#ifndef RODIN_PETSC_MATH_MATRIX_H
#define RODIN_PETSC_MATH_MATRIX_H

#include <petsc.h>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Matrix.h"

namespace Rodin::Math
{
  template <class Scalar>
  void axpy(::Mat& y, Scalar alpha, const ::Mat& x, MatStructure str)
  {
    PetscErrorCode ierr;
    ierr = MatAXPY(y, alpha, x, str);
    assert(ierr == PETSC_SUCCESS);
  }

  template <class Scalar>
  void axpy(::Mat& y, Scalar alpha, const ::Mat& x)
  {
    axpy(y, alpha, x, DIFFERENT_NONZERO_PATTERN);
  }
}

#endif

