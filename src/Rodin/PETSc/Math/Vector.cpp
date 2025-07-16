#include "Vector.h"
#include <petscsystypes.h>
#include <petscvec.h>

namespace Rodin::PETSc
{
  Vector::Vector(const boost::mpi::communicator& comm)
  {
    PetscErrorCode ierr;
    ierr = VecCreate(comm, &m_vec);
    assert(ierr == PETSC_SUCCESS);
  }

  Vector::Vector(const Vector& other)
    : Object(other),
      m_vec(other.m_vec)
  {
    PetscErrorCode ierr;
    ierr = VecDuplicate(other.m_vec, &m_vec);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecCopy(other.m_vec, m_vec);
    assert(ierr == PETSC_SUCCESS);
  }

  Vector::Vector(Vector&& other) noexcept
    : Object(std::move(other)),
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

  Vector& Vector::operator=(const Vector& other)
  {
    if (this != &other)
    {
      Parent::operator=(other);
      PetscErrorCode ierr;
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
      m_vec = std::exchange(other.m_vec, PETSC_NULLPTR);
    }
    return *this;
  }

  Vector& Vector::operator+=(const Vector& other)
  {
    PetscErrorCode ierr = VecAXPY(m_vec, 1.0, other.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator-=(const Vector& other)
  {
    PetscErrorCode ierr = VecAXPY(m_vec, -1.0, other.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator*=(const Vector& rhs)
  {
    PetscErrorCode ierr = VecPointwiseMult(m_vec, m_vec, rhs.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator/=(const Vector& rhs)
  {
    PetscErrorCode ierr = VecPointwiseDivide(m_vec, m_vec, rhs.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator*=(const PetscScalar& rhs)
  {
    PetscErrorCode ierr = VecScale(m_vec, rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator/=(const PetscScalar& rhs)
  {
    assert(rhs != PetscScalar(0));
    PetscErrorCode ierr = VecScale(m_vec, 1.0 / rhs);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator+=(const PetscScalar& other)
  {
    PetscErrorCode ierr = VecShift(m_vec, other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::operator-=(const PetscScalar& other)
  {
    PetscErrorCode ierr = VecShift(m_vec, -other);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setType(const ::VecType& type)
  {
    PetscErrorCode ierr = VecSetType(m_vec, type);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setFromOptions()
  {
    PetscErrorCode ierr = VecSetFromOptions(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setSizes(PetscInt localSize, PetscInt globalSize)
  {
    PetscErrorCode ierr = VecSetSizes(m_vec, localSize, globalSize);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::zeroEntries()
  {
    PetscErrorCode ierr = VecZeroEntries(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::setValue(PetscInt idx, const PetscScalar& value, ::InsertMode mode)
  {
    PetscErrorCode ierr = VecSetValue(m_vec, idx, value, mode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::assemblyBegin()
  {
    PetscErrorCode ierr = VecAssemblyBegin(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::assemblyEnd()
  {
    PetscErrorCode ierr = VecAssemblyEnd(m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::axpy(const PetscScalar& alpha, const Vector& x)
  {
    PetscErrorCode ierr = VecAXPY(m_vec, alpha, x.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::aypx(const PetscScalar& beta, const Vector& x)
  {
    PetscErrorCode ierr = VecAYPX(m_vec, beta, x.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::axpby(const PetscScalar& alpha, const PetscScalar& beta, const Vector& x)
  {
    PetscErrorCode ierr = VecAXPBY(m_vec, alpha, beta, x.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::waxpy(const PetscScalar& alpha, const Vector& x, const Vector& y)
  {
    PetscErrorCode ierr = VecWAXPY(m_vec, alpha, x.m_vec, y.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::axpbypcz(
      const PetscScalar& alpha, const Vector& x, const PetscScalar& beta, const Vector& y, const PetscScalar& gamma)
  {
    PetscErrorCode ierr = VecAXPBYPCZ(m_vec, alpha, beta, gamma, x.m_vec, y.m_vec);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::getArrayRead(const PetscScalar *a[]) const
  {
    PetscErrorCode ierr;
    ierr = VecGetArrayRead(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::restoreArrayRead(const PetscScalar *a[]) const
  {
    PetscErrorCode ierr;
    ierr = VecRestoreArrayRead(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::getArrayWrite(PetscScalar *a[])
  {
    PetscErrorCode ierr;
    ierr = VecGetArray(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::restoreArrayWrite(PetscScalar *a[])
  {
    PetscErrorCode ierr;
    ierr = VecRestoreArray(m_vec, a);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::getLocalSize(PetscInt* size) const
  {
    PetscErrorCode ierr = VecGetLocalSize(m_vec, size);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::ghostUpdateBegin(::InsertMode mode, ::ScatterMode scatterMode)
  {
    PetscErrorCode ierr;
    ierr = VecGhostUpdateBegin(m_vec, mode, scatterMode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  Vector& Vector::ghostUpdateEnd(::InsertMode mode, ::ScatterMode scatterMode)
  {
    PetscErrorCode ierr;
    ierr = VecGhostUpdateEnd(m_vec, mode, scatterMode);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  const Vector& Vector::getComm(MPI_Comm* comm) const noexcept
  {
    PetscErrorCode ierr;
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(m_vec), comm);
    assert(ierr == PETSC_SUCCESS);
    return *this;
  }

  ::PetscObject& Vector::getHandle() noexcept
  {
    return reinterpret_cast<::PetscObject&>(m_vec);
  }

  const ::PetscObject& Vector::getHandle() const noexcept
  {
    return reinterpret_cast<const ::PetscObject&>(m_vec);
  }
}
