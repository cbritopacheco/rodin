#include <petscmat.h>
#include <petscsystypes.h>

#include "Matrix.h"

namespace Rodin::PETSc
{
  Matrix::Matrix(const boost::mpi::communicator& comm)
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
    PetscErrorCode ierr;
    ierr = MatSetType(m_mat, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setFromOptions()
  {
    PetscErrorCode ierr;
    ierr = MatSetFromOptions(m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setSizes(
      const PetscInt& localRows, const PetscInt& localCols,
      const PetscInt& globalRows, const PetscInt& globalCols)
  {
    PetscErrorCode ierr;
    ierr = MatSetSizes(m_mat, localRows, localCols, globalRows, globalCols);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::zeroEntries()
  {
    PetscErrorCode ierr = MatZeroEntries(m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setUp()
  {
    PetscErrorCode ierr = MatSetUp(m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::setValue(
      const PetscInt& row, const PetscInt& col,
      const PetscScalar& value, ::InsertMode mode)
  {
    PetscErrorCode ierr = MatSetValue(m_mat, row, col, value, mode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::assemblyBegin(::MatAssemblyType type)
  {
    PetscErrorCode ierr = MatAssemblyBegin(m_mat, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::assemblyEnd(::MatAssemblyType type)
  {
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
    PetscErrorCode ierr = MatAXPY(m_mat, 1.0, other.m_mat, UNKNOWN_NONZERO_PATTERN);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator-=(const Matrix& other)
  {
    PetscErrorCode ierr = MatAXPY(m_mat, -1.0, other.m_mat, UNKNOWN_NONZERO_PATTERN);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator*=(const Matrix& rhs)
  {
    PetscErrorCode ierr =
      MatMatMult(m_mat, rhs.m_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator*=(const PetscScalar& rhs)
  {
    PetscErrorCode ierr = MatScale(m_mat, rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator/=(const PetscScalar& rhs)
  {
    assert(rhs != PetscScalar(0));
    PetscErrorCode ierr = MatScale(m_mat, 1.0 / rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator+=(const PetscScalar& other)
  {
    PetscErrorCode ierr = MatShift(m_mat, other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::operator-=(const PetscScalar& other)
  {
    PetscErrorCode ierr = MatShift(m_mat, -other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Matrix& Matrix::zeroRowsColumns(
      PetscInt n, const PetscInt rows[], PetscScalar diag, Vector& x, Vector& b)
  {
    ::Vec& xh = reinterpret_cast<::Vec&>(x.getHandle());
    ::Vec& bh = reinterpret_cast<::Vec&>(b.getHandle());
    MatZeroRowsColumns(m_mat, n, rows, diag, xh, bh);
    return *this;
  }

  const Matrix& Matrix::duplicate(Matrix& dst) const
  {
    PetscErrorCode ierr;
    ierr = MatDuplicate(m_mat, MAT_COPY_VALUES, &dst.m_mat);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Matrix& Matrix::getComm(MPI_Comm& comm) const noexcept
  {
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

