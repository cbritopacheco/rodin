#include "Vector.h"

namespace Rodin::PETSc
{
  Vector::Vector(Vector&& other) noexcept
  {
    m_vec = other.m_vec;
    other.m_vec = PETSC_NULLPTR;
  }

  Vector::~Vector()
  {
    VecDestroy(&m_vec);
  }

  PetscErrorCode Vector::setSizes(PetscInt localSize, PetscInt globalSize)
  {
    return VecSetSizes(m_vec, localSize, globalSize);
  }

  PetscErrorCode Vector::setFromOptions()
  {
    return VecSetFromOptions(m_vec);
  }

  PetscErrorCode Vector::zeroEntries()
  {
    return VecZeroEntries(m_vec);
  }

  PetscErrorCode Vector::setValue(PetscInt idx, PetscScalar value, InsertMode mode)
  {
    return VecSetValue(m_vec, idx, value, mode);
  }

  PetscErrorCode Vector::assemblyBegin()
  {
    return VecAssemblyBegin(m_vec);
  }

  PetscErrorCode Vector::assemblyEnd()
  {
    return VecAssemblyEnd(m_vec);
  }
}
