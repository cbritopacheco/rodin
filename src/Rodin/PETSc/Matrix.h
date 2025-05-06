#ifndef RODIN_PETSC_MATRIX_H
#define RODIN_PETSC_MATRIX_H

#include <petsc.h>

namespace Rodin::PETSc
{
  class Matrix
  {
    public:
      /*
       * @brief Constructor for creating a PETSc matrix.
       */
      Matrix() = default;

      Matrix(const Matrix& other) = delete;

      Matrix& operator=(const Matrix& other) = delete;

      /*
       * @brief Move constructor for transferring ownership of a matrix.
       * @param other The matrix to move from.
       *
       * Transfers ownership of the matrix from other to this object.
       */
      Matrix(Matrix&& other) noexcept;

      Matrix& operator=(Matrix&& other) noexcept;

      ~Matrix();

      PetscErrorCode setSizes(PetscInt localRows, PetscInt localCols, PetscInt globalRows, PetscInt globalCols);

      PetscErrorCode setFromOptions();

      PetscErrorCode setPreallocationMPIAIJ(
          PetscInt d_nz, const PetscInt* d_nz_indices, PetscInt o_nz, const PetscInt o_nz_indices[]);

      PetscErrorCode setUp();

      PetscErrorCode zeroEntries();

      PetscErrorCode setValue(PetscInt row, PetscInt col, PetscScalar value, InsertMode mode);

      PetscErrorCode setValues(PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv);

      PetscErrorCode assemblyBegin(MatAssemblyType type);

      PetscErrorCode assemblyEnd(MatAssemblyType type);

      PetscErrorCode duplicate(MatDuplicateOption op, Mat* newmat);

      const Mat& getMat() const;

      Mat& getMat();

    private:
      PetscErrorCode m_ierr;
      Mat m_mat;
  };
}

#endif
