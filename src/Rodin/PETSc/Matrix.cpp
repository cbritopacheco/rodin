#include "Matrix.h"

namespace Rodin::PETSc
{
  Matrix::Matrix(Matrix&& other) noexcept
  {
    m_mat = other.m_mat;
    other.m_mat = PETSC_NULLPTR;
  }

  Matrix& Matrix::operator=(Matrix&& other) noexcept
  {
    if (this != &other)
    {
      PetscErrorCode ierr = MatDestroy(&m_mat);
      if (ierr == PETSC_SUCCESS)
      {
        m_mat = other.m_mat;
        other.getMat() = PETSC_NULLPTR;
      }
    }
    return *this;
  }

  Matrix::~Matrix()
  {
    MatDestroy(&m_mat);
  }

  Mat& Matrix::getMat()
  {
    return m_mat;
  }

  const Mat& Matrix::getMat() const
  {
    return m_mat;
  }

  PetscErrorCode Matrix::setSizes(
      PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols)
  {
    return MatSetSizes(getMat(), localRows, localCols, globalRows, globalCols);
  }

  PetscErrorCode Matrix::setFromOptions()
  {
    return MatSetFromOptions(getMat());
  }

  PetscErrorCode Matrix::setPreallocationMPIAIJ(
      PetscInt d_nz, const PetscInt d_nz_indices[], PetscInt o_nz, const PetscInt o_nz_indices[])
  {
    return MatMPIAIJSetPreallocation(getMat(), d_nz, d_nz_indices, o_nz, o_nz_indices);
  }

  PetscErrorCode Matrix::setUp()
  {
    return MatSetUp(getMat());
  }

  PetscErrorCode Matrix::zeroEntries()
  {
    return MatZeroEntries(getMat());
  }

  PetscErrorCode Matrix::setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode)
  {
    return MatSetValue(getMat(), row, col, value, mode);
  }

  PetscErrorCode Matrix::setValues(
      PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
  {
    return MatSetValues(getMat(), m, idxm, n, idxn, v, addv);
  }

  PetscErrorCode Matrix::assemblyBegin(MatAssemblyType type)
  {
    return MatAssemblyBegin(getMat(), type);
  }

  PetscErrorCode Matrix::assemblyEnd(MatAssemblyType type)
  {
    return MatAssemblyEnd(getMat(), type);
  }

  PetscErrorCode Matrix::duplicate(MatDuplicateOption op, Mat* newmat)
  {
    return MatDuplicate(getMat(), op, newmat);
  }
}
