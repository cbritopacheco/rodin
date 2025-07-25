#include <mpi.h>
#include <petscmat.h>
#include <petscsystypes.h>

#include "Matrix.h"

namespace Rodin::PETSc
{
  Matrix::Matrix()
    : SeqAIJ(*this),
      MPIAIJ(*this),
      m_mat(PETSC_NULLPTR)
  {}

  Matrix::Matrix(MPI_Comm comm)
    : SeqAIJ(*this),
      MPIAIJ(*this)
  {
    PetscErrorCode ierr;
    ierr = MatCreate(comm, &m_mat);
    assert(ierr == PETSC_SUCCESS);
  }

  Matrix::Matrix(const Matrix& other)
    : Object(other),
      SeqAIJ(*this),
      MPIAIJ(*this),
      m_mat(other.m_mat)
  {
    PetscErrorCode ierr;
    ierr = MatDuplicate(other.m_mat, MAT_COPY_VALUES, &m_mat);
    assert(ierr == PETSC_SUCCESS);
  }

  Matrix::Matrix(Matrix&& other) noexcept
    : Object(std::move(other)),
      SeqAIJ(*this),
      MPIAIJ(*this),
      m_mat(std::exchange(other.m_mat, PETSC_NULLPTR))
  {}

  Matrix::~Matrix()
  {
    if (m_mat != PETSC_NULLPTR)
    {
      PetscErrorCode ierr = MatDestroy(&m_mat);
      assert(ierr == PETSC_SUCCESS);
      m_mat = PETSC_NULLPTR;
    }
  }

  Matrix& Matrix::setType(const ::MatType& type)
  {
    assert(m_mat);
    PetscErrorCode ierr;
    ierr = MatSetType(m_mat, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setFromOptions()
  {
    assert(m_mat);
    PetscErrorCode ierr;
    ierr = MatSetFromOptions(m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setSizes(
      const PetscInt& localRows, const PetscInt& localCols,
      const PetscInt& globalRows, const PetscInt& globalCols)
  {
    assert(m_mat);
    PetscErrorCode ierr;
    ierr = MatSetSizes(m_mat, localRows, localCols, globalRows, globalCols);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::zeroEntries()
  {
    assert(m_mat);
    PetscErrorCode ierr = MatZeroEntries(m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setUp()
  {
    assert(m_mat);
    PetscErrorCode ierr = MatSetUp(m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setValue(
      const PetscInt& row, const PetscInt& col,
      const PetscScalar& value, ::InsertMode mode)
  {
    assert(m_mat);
    PetscErrorCode ierr = MatSetValue(m_mat, row, col, value, mode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::assemblyBegin(::MatAssemblyType type)
  {
    assert(m_mat);
    PetscErrorCode ierr = MatAssemblyBegin(m_mat, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::assemblyEnd(::MatAssemblyType type)
  {
    assert(m_mat);
    PetscErrorCode ierr = MatAssemblyEnd(m_mat, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator=(const Matrix& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      PetscErrorCode ierr;
      ierr = MatDuplicate(other.m_mat, MAT_COPY_VALUES, &m_mat);
      assert(ierr == PETSC_SUCCESS);
    }
    return *this;
  }

  Matrix& Matrix::operator=(Matrix&& other) noexcept
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      m_mat = std::exchange(other.m_mat, PETSC_NULLPTR);
    }
    return *this;
  }

  Matrix& Matrix::operator+=(const Matrix& other)
  {
    assert(m_mat && other.m_mat);
    PetscErrorCode ierr = MatAXPY(m_mat, 1.0, other.m_mat, UNKNOWN_NONZERO_PATTERN);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator-=(const Matrix& other)
  {
    assert(m_mat && other.m_mat);
    PetscErrorCode ierr = MatAXPY(m_mat, -1.0, other.m_mat, UNKNOWN_NONZERO_PATTERN);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator*=(const Matrix& rhs)
  {
    assert(m_mat && rhs.m_mat);
    PetscErrorCode ierr =
      MatMatMult(m_mat, rhs.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator*=(const PetscScalar& rhs)
  {
    assert(m_mat);
    PetscErrorCode ierr = MatScale(m_mat, rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator/=(const PetscScalar& rhs)
  {
    assert(m_mat);
    assert(rhs != PetscScalar(0));
    PetscErrorCode ierr = MatScale(m_mat, 1.0 / rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator+=(const PetscScalar& other)
  {
    assert(m_mat);
    PetscErrorCode ierr = MatShift(m_mat, other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator-=(const PetscScalar& other)
  {
    assert(m_mat);
    PetscErrorCode ierr = MatShift(m_mat, -other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::zeroRowsColumns(
      PetscInt n, const PetscInt rows[], PetscScalar diag, Vector& x, Vector& b)
  {
    assert(m_mat && x.getHandle() && b.getHandle());
    ::Vec& xh = reinterpret_cast<::Vec&>(x.getHandle());
    ::Vec& bh = reinterpret_cast<::Vec&>(b.getHandle());
    MatZeroRowsColumns(m_mat, n, rows, diag, xh, bh);
    return *this;
  }

  const Matrix& Matrix::duplicate(Matrix& dst) const
  {
    assert(m_mat);
    PetscErrorCode ierr;
    ierr = MatDuplicate(m_mat, MAT_COPY_VALUES, &dst.m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Matrix& Matrix::getComm(MPI_Comm& comm) const noexcept
  {
    assert(m_mat);
    PetscErrorCode ierr;
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(m_mat), &comm);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  ::PetscObject& Matrix::getHandle() noexcept
  {
    return reinterpret_cast<::PetscObject&>(m_mat);
  }

  const ::PetscObject& Matrix::getHandle() const noexcept
  {
    return reinterpret_cast<const ::PetscObject&>(m_mat);
  }
}

