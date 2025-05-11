#ifndef RODIN_PETSC_MATH_VECTOR_H
#define RODIN_PETSC_MATH_VECTOR_H

#include <petsc.h>

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"

namespace Rodin::Math
{
  void duplicate(const ::Vec& src, ::Vec& dst);

  template <class Scalar, int Size>
  void duplicate(const ::Vec& src, Eigen::Vector<Scalar, Size>& dst)
  {
    PetscErrorCode ierr;
    PetscInt sz;
    ierr = VecGetLocalSize(src, &sz);
    assert(ierr == PETSC_SUCCESS);
    dst.resize(sz);
  }

  void copy(const ::Vec& src, ::Vec& dst);

  template <class Scalar, int Size>
  void copy(const ::Vec& src, Eigen::Vector<Scalar, Size>& dst)
  {
    PetscErrorCode ierr;
    const PetscScalar* vec = nullptr;
    ierr = VecGetArrayRead(src, &vec);
    assert(ierr == PETSC_SUCCESS);
    std::copy(vec, vec + dst.size(), dst.data());
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

  template <class Scalar, int Rows, int Cols, int Options>
  void duplicate(const ::Mat& src, Eigen::Matrix<Scalar, Rows, Cols, Options>& dst)
  {
    PetscErrorCode ierr;
    PetscInt rows, cols;
    ierr = MatGetLocalSize(src, &rows, &cols);
    assert(ierr == PETSC_SUCCESS);
    dst.resize(rows, cols);
  }

  template <class Scalar>
  void duplicate(const ::Mat& src, Math::SparseMatrix<Scalar>& dst)
  {
    PetscErrorCode ierr;
    PetscInt       rows, cols;
    ierr = MatGetLocalSize(src, &rows, &cols);
    assert(ierr == PETSC_SUCCESS);
    dst.resize(static_cast<int>(rows), static_cast<int>(cols));
    for (PetscInt i = 0; i < rows; ++i)
    {
      PetscInt           nnz;
      const PetscInt    *colsIdx = nullptr;
      const PetscScalar *unused  = nullptr;
      ierr = MatGetRow(src, i, &nnz, &colsIdx, &unused);
      assert(ierr == PETSC_SUCCESS);
      for (PetscInt k = 0; k < nnz; ++k)
        dst.insert(i, colsIdx[k]) = Scalar(0);
      ierr = MatRestoreRow(src, i, &nnz, &colsIdx, &unused);
      assert(ierr == PETSC_SUCCESS);
    }
    dst.makeCompressed();
  }

  template <class Scalar>
  void copy(const ::Mat& src, Math::SparseMatrix<Scalar>& dst)
  {
    PetscErrorCode ierr;
    PetscInt       rows, cols;
    ierr = MatGetLocalSize(src, &rows, &cols);
    assert(ierr == PETSC_SUCCESS);
    for (PetscInt i = 0; i < rows; ++i)
    {
      PetscInt           nnz;
      const PetscInt    *colsIdx = nullptr;
      const PetscScalar *vals    = nullptr;
      ierr = MatGetRow(src, i, &nnz, &colsIdx, &vals);
      assert(ierr == PETSC_SUCCESS);
      for (PetscInt k = 0; k < nnz; ++k)
        dst.coeffRef(i, colsIdx[k]) = static_cast<Scalar>(vals[k]);
      ierr = MatRestoreRow(src, i, &nnz, &colsIdx, &vals);
      assert(ierr == PETSC_SUCCESS);
    }
  }
}

#endif
