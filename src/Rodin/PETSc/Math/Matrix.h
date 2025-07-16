#ifndef RODIN_PETSC_MATH_MATRIX_H
#define RODIN_PETSC_MATH_MATRIX_H

#include <boost/mpi/communicator.hpp>
#include <petsc.h>
#include <petscmat.h>
#include <petscsystypes.h>

#include "Rodin/FormLanguage/Traits.h"

#include "Rodin/PETSc/Math/Vector.h"

#include "Rodin/PETSc/Object.h"

namespace Rodin::PETSc
{
  class Matrix : public Object
  {
    public:
      class SeqAIJ
      {
        public:
          constexpr
          SeqAIJ(Matrix& mat)
            : m_mat(mat)
          {}

          Matrix& setPreallocation(PetscInt nz, const PetscInt nnz[])
          {
            PetscErrorCode ierr;
            ierr = MatSeqAIJSetPreallocation(m_mat.m_mat, nz, nnz);
            assert(ierr == PETSC_SUCCESS);
            return m_mat;
          }

        private:
          Matrix& m_mat; ///< Reference to the parent Matrix object.
      } SeqAIJ;

      class MPIAIJ
      {
        public:
          constexpr
          MPIAIJ(Matrix& mat)
            : m_mat(mat)
          {}

          Matrix& setPreallocation(
              PetscInt d_nz, const PetscInt d_nnz[], PetscInt o_nz, const PetscInt o_nnz[])
          {
            PetscErrorCode ierr;
            ierr = MatMPIAIJSetPreallocation(m_mat.m_mat, d_nz, d_nnz, o_nz, o_nnz);
            assert(ierr == PETSC_SUCCESS);
            return m_mat;
          }

        private:
          Matrix& m_mat; ///< Reference to the parent Matrix object.
      } MPIAIJ;

      using Parent = Object;

      Matrix(const boost::mpi::communicator& comm);

      Matrix(const Matrix& other);

      Matrix(Matrix&& other) noexcept;

      virtual ~Matrix() override;

      Matrix& operator=(const Matrix& other);

      Matrix& operator=(Matrix&& other) noexcept;

      Matrix& operator+=(const Matrix& other);

      Matrix& operator-=(const Matrix& other);

      Matrix& operator*=(const Matrix& rhs);

      Matrix& operator+=(const PetscScalar& other);

      Matrix& operator-=(const PetscScalar& other);

      Matrix& operator*=(const PetscScalar& rhs);

      Matrix& operator/=(const PetscScalar& rhs);

      const Matrix& duplicate(Matrix& dst) const;

      Matrix& setType(const ::MatType& type);

      Matrix& setFromOptions();

      Matrix& setSizes(
          const PetscInt& localRows, const PetscInt& localCols,
          const PetscInt& globalRows, const PetscInt& globalCols);

      Matrix& zeroEntries();

      Matrix& setUp();

      Matrix& setValue(
          const PetscInt& row, const PetscInt& col,
          const PetscScalar& value, ::InsertMode mode);

      Matrix& assemblyBegin(::MatAssemblyType type);

      Matrix& assemblyEnd(::MatAssemblyType type);

      Matrix& zeroRowsColumns(
        PetscInt n, const PetscInt rows[], PetscScalar diag, Vector& x, Vector& b);

      const Matrix& getComm(MPI_Comm& comm) const noexcept;

      ::PetscObject& getHandle() noexcept override;

      const ::PetscObject& getHandle() const noexcept;

    private:
      ::Mat m_mat;
  };
}

namespace Rodin::FormLanguage
{
  template <>
  struct Traits<PETSc::Matrix>
  {
    using ScalarType = PetscScalar;
  };
}

#endif

