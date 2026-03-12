#include <cassert>
#include <petsc.h>

#include "SNES.h"

namespace Rodin::Solver
{
  SNES::SNES(MPI_Comm comm)
    : m_snes(PETSC_NULLPTR),
      m_type(PETSC_NULLPTR),
      m_abstol(PETSC_DECIDE),
      m_rtol(PETSC_DECIDE),
      m_stol(PETSC_DECIDE),
      m_maxIt(PETSC_DECIDE),
      m_maxF(PETSC_DECIDE),
      m_ksp(PETSC_NULLPTR)
  {
    PetscErrorCode ierr = SNESCreate(comm, &m_snes);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
  }

  SNES::~SNES()
  {
    if (m_snes)
    {
      PetscErrorCode ierr = SNESDestroy(&m_snes);
      assert(ierr == PETSC_SUCCESS);
      m_snes = PETSC_NULLPTR;
      (void) ierr;
    }
  }

  SNES& SNES::setType(::SNESType type) noexcept
  {
    m_type = type;
    return *this;
  }

  SNES& SNES::setTolerances(
      PetscReal abstol, PetscReal rtol, PetscReal stol, PetscInt maxIt, PetscInt maxF) noexcept
  {
    m_abstol = abstol;
    m_rtol = rtol;
    m_stol = stol;
    m_maxIt = maxIt;
    m_maxF = maxF;
    return *this;
  }

  SNES& SNES::setKSP(KSPType ksp) noexcept
  {
    m_ksp = ksp;
    return *this;
  }

  void SNES::solve(VectorType b, VectorType x)
  {
    PetscErrorCode ierr;

    ierr = SNESSetType(m_snes, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetTolerances(m_snes, m_abstol, m_rtol, m_stol, m_maxIt, m_maxF);
    assert(ierr == PETSC_SUCCESS);

    if (m_ksp)
    {
      ierr = SNESSetKSP(m_snes, m_ksp);
      assert(ierr == PETSC_SUCCESS);
    }

    ierr = SNESSetFromOptions(m_snes);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSolve(m_snes, b, x);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
  }

  ::SNES& SNES::getHandle() noexcept
  {
    return m_snes;
  }

  const ::SNES& SNES::getHandle() const noexcept
  {
    return m_snes;
  }
}
