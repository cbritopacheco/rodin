#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

#include <petsc.h>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Math
{
  void duplicate(const ::Vec& src, ::Vec& dest);

  template <class Scalar, int Size>
  void duplicate(const ::Vec& src, Eigen::Vector<Scalar, Size>& dest)
  {
    PetscErrorCode ierr;
    PetscInt sz;
    ierr = VecGetLocalSize(src, &sz);
    assert(ierr == PETSC_SUCCESS);
    dest.resize(sz);
  }

  void copy(const ::Vec& src, ::Vec& dest);

  template <class Scalar, int Size>
  void copy(const ::Vec& src, Eigen::Vector<Scalar, Size>& dest)
  {
    PetscErrorCode ierr;
    const PetscScalar* vec = nullptr;
    ierr = VecGetArrayRead(src, &vec);
    assert(ierr == PETSC_SUCCESS);
    std::copy(vec, vec + dest.size(), dest.data());
    ierr = VecRestoreArrayRead(src, &vec);
    assert(ierr == PETSC_SUCCESS);
  }

  template <class Scalar, int Size>
  void axpy(Eigen::Vector<Scalar, Size>& y, Scalar alpha, const ::Vec& x)
  {
    if (alpha == Scalar(0))
      return;
    PetscErrorCode ierr;
    PetscInt sz;
    ierr = VecGetLocalSize(x, &sz);
    assert(ierr == PETSC_SUCCESS);
    assert(y.size() == sz);
    const PetscScalar* vec = nullptr;
    ierr = VecGetArrayRead(x, &vec);
    assert(ierr == PETSC_SUCCESS);
    Eigen::Map<const Eigen::Vector<PetscScalar, Size>, Eigen::Unaligned> xm(vec, sz);
    if constexpr (std::is_same_v<Scalar, PetscScalar>)
      y.noalias() += alpha * xm;
    else
      y.noalias() += alpha * xm.template cast<Scalar>();
    ierr = VecRestoreArrayRead(x, &vec);
    assert(ierr == PETSC_SUCCESS);
  }

  template <class Scalar>
  void axpy(::Vec& y, Scalar alpha, const ::Vec& x)
  {
    PetscErrorCode ierr;
    ierr = VecAXPY(y, alpha, x);
    assert(ierr == PETSC_SUCCESS);
  }
}

#endif
