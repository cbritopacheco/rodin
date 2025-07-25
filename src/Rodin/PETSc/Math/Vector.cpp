#include "Vector.h"
#include <petscsystypes.h>
#include <petscvec.h>

namespace Rodin::PETSc
{
  Vector::Vector()
    : MPI(*this),
      m_vec(PETSC_NULLPTR)
  {}

  Vector::Vector(MPI_Comm comm)
    : MPI(*this)
  {
    PetscErrorCode ierr;
    ierr = VecCreate(comm, &m_vec);
    assert(ierr == PETSC_SUCCESS);
  }

  Vector::Vector(const Vector& other)
    : Parent(other),
      MPI(*this),
      m_vec(other.m_vec)
  {
    PetscErrorCode ierr;
    ierr = VecDuplicate(other.m_vec, &m_vec);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCopy(other.m_vec, m_vec);
    assert(ierr == PETSC_SUCCESS);
  }

  Vector::Vector(Vector&& other) noexcept
    : Parent(std::move(other)),
      MPI(*this),
      m_vec(std::exchange(other.m_vec, PETSC_NULLPTR))
  {}

  Vector::~Vector()
  {
    if (m_vec != PETSC_NULLPTR)
    {
      PetscErrorCode ierr = VecDestroy(&m_vec);
      assert(ierr == PETSC_SUCCESS);
      m_vec = PETSC_NULLPTR;
    }
  }

  Vector::operator bool() const noexcept
  {
    return m_vec != PETSC_NULLPTR;
  }

  Vector& Vector::operator=(const Vector& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      PetscErrorCode ierr;
      if (m_vec != PETSC_NULLPTR)
      {
        ierr = VecDestroy(&m_vec);
        assert(ierr == PETSC_SUCCESS);
      }
      ierr = VecDuplicate(other.m_vec, &m_vec);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecCopy(other.m_vec, m_vec);
      assert(ierr == PETSC_SUCCESS);
    }
    return *this;
  }

  Vector& Vector::operator=(Vector&& other) noexcept
  {
    if (this != &other)
    {
      Parent::operator=(std::move(other));
      if (m_vec != PETSC_NULLPTR)
      {
        PetscErrorCode ierr = VecDestroy(&m_vec);
        assert(ierr == PETSC_SUCCESS);
      }
      m_vec = std::exchange(other.m_vec, PETSC_NULLPTR);
    }
    return *this;
  }

  Vector& Vector::operator+=(const Vector& other)
  {
    assert(m_vec && other.m_vec);
    PetscErrorCode ierr = VecAXPY(m_vec, 1.0, other.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator-=(const Vector& other)
  {
    assert(m_vec && other.m_vec);
    PetscErrorCode ierr = VecAXPY(m_vec, -1.0, other.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator*=(const Vector& rhs)
  {
    assert(m_vec && rhs.m_vec);
    PetscErrorCode ierr = VecPointwiseMult(m_vec, m_vec, rhs.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator/=(const Vector& rhs)
  {
    assert(m_vec && rhs.m_vec);
    PetscErrorCode ierr = VecPointwiseDivide(m_vec, m_vec, rhs.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator*=(const PetscScalar& rhs)
  {
    assert(m_vec);
    PetscErrorCode ierr = VecScale(m_vec, rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator/=(const PetscScalar& rhs)
  {
    assert(m_vec);
    assert(rhs != PetscScalar(0));
    PetscErrorCode ierr = VecScale(m_vec, 1.0 / rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator+=(const PetscScalar& other)
  {
    assert(m_vec);
    PetscErrorCode ierr = VecShift(m_vec, other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator-=(const PetscScalar& other)
  {
    assert(m_vec);
    PetscErrorCode ierr = VecShift(m_vec, -other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::duplicate(Vector& dst) const
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecDuplicate(m_vec, &dst.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::copy(Vector& dst) const
  {
    assert(m_vec);
    assert(dst.m_vec);
    PetscErrorCode ierr;
    ierr = VecCopy(m_vec, dst.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::create(MPI_Comm comm)
  {
    PetscErrorCode ierr;
    ierr = VecCreate(comm, &m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::destroy()
  {
    PetscErrorCode ierr = VecDestroy(&m_vec);
    assert(ierr == PETSC_SUCCESS);
    m_vec = PETSC_NULLPTR;
    return *this;
  }

  Vector& Vector::setType(const ::VecType& type)
  {
    assert(m_vec);
    PetscErrorCode ierr = VecSetType(m_vec, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setFromOptions()
  {
    assert(m_vec);
    PetscErrorCode ierr = VecSetFromOptions(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setSizes(PetscInt localSize, PetscInt globalSize)
  {
    assert(m_vec);
    PetscErrorCode ierr = VecSetSizes(m_vec, localSize, globalSize);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::zeroEntries()
  {
    assert(m_vec);
    PetscErrorCode ierr = VecZeroEntries(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setValue(PetscInt idx, const PetscScalar& value, ::InsertMode mode)
  {
    assert(m_vec);
    PetscErrorCode ierr = VecSetValue(m_vec, idx, value, mode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::assemblyBegin()
  {
    assert(m_vec);
    PetscErrorCode ierr = VecAssemblyBegin(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::assemblyEnd()
  {
    assert(m_vec);
    PetscErrorCode ierr = VecAssemblyEnd(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::axpy(const PetscScalar& alpha, const Vector& x)
  {
    assert(m_vec && x.m_vec);
    PetscErrorCode ierr = VecAXPY(m_vec, alpha, x.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::aypx(const PetscScalar& beta, const Vector& x)
  {
    assert(m_vec && x.m_vec);
    PetscErrorCode ierr = VecAYPX(m_vec, beta, x.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::axpby(const PetscScalar& alpha, const PetscScalar& beta, const Vector& x)
  {
    assert(m_vec && x.m_vec);
    PetscErrorCode ierr = VecAXPBY(m_vec, alpha, beta, x.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::waxpy(const PetscScalar& alpha, const Vector& x, const Vector& y)
  {
    assert(m_vec && x.m_vec && y.m_vec);
    PetscErrorCode ierr = VecWAXPY(m_vec, alpha, x.m_vec, y.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::axpbypcz(
      const PetscScalar& alpha, const Vector& x, const PetscScalar& beta, const Vector& y, const PetscScalar& gamma)
  {
    assert(m_vec && x.m_vec && y.m_vec);
    PetscErrorCode ierr = VecAXPBYPCZ(m_vec, alpha, beta, gamma, x.m_vec, y.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::getArrayRead(const PetscScalar *a[]) const
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecGetArrayRead(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::restoreArrayRead(const PetscScalar *a[]) const
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecRestoreArrayRead(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::getArrayWrite(PetscScalar *a[])
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecGetArray(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::restoreArrayWrite(PetscScalar *a[])
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecRestoreArray(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::getLocalSize(PetscInt* size) const
  {
    assert(m_vec);
    PetscErrorCode ierr = VecGetLocalSize(m_vec, size);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::ghostUpdateBegin(::InsertMode mode, ::ScatterMode scatterMode)
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecGhostUpdateBegin(m_vec, mode, scatterMode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::ghostUpdateEnd(::InsertMode mode, ::ScatterMode scatterMode)
  {
    assert(m_vec);
    PetscErrorCode ierr;
    ierr = VecGhostUpdateEnd(m_vec, mode, scatterMode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  ::Vec& Vector::getHandle() noexcept
  {
    return m_vec;
  }

  const ::Vec& Vector::getHandle() const noexcept
  {
    return m_vec;
  }
}
