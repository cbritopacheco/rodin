/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <petsc.h>

#include "Rodin/Variational/Problem.h"

#include "Rodin/PETSc/Math/LinearSystem.h"

#include "KSP.h"

namespace Rodin::Solver
{
  KSP::KSP(ProblemBaseType& pb)
    : Parent(pb),
      m_ksp(PETSC_NULLPTR),
      m_type(PETSC_NULLPTR),
      m_rtol(PETSC_DECIDE),
      m_abstol(PETSC_DECIDE),
      m_dtol(PETSC_DECIDE),
      m_maxIt(PETSC_DECIDE)
  {
    PetscErrorCode ierr;
    ierr = KSPCreate(pb.getLinearSystem().getCommunicator(), &m_ksp);
    assert(ierr == PETSC_SUCCESS);
  }

  KSP::~KSP()
  {
    if (m_ksp)
    {
      PetscErrorCode ierr = KSPDestroy(&m_ksp);
      assert(ierr == PETSC_SUCCESS);
      m_ksp = PETSC_NULLPTR;
    }
  }

  ::KSP& KSP::getHandle() noexcept
  {
    return m_ksp;
  }

  const ::KSP& KSP::getHandle() const noexcept
  {
    return m_ksp;
  }

  void KSP::solve(PETSc::Math::LinearSystem& axb)
  {
    auto& [a, x, b] = axb;

    PetscErrorCode ierr;

    ierr = KSPSetInitialGuessNonzero(m_ksp, PETSC_TRUE);
    assert(ierr == PETSC_SUCCESS);

    ierr = KSPSetType(m_ksp, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = KSPSetTolerances(m_ksp, m_rtol, m_abstol, m_dtol, m_maxIt);
    assert(ierr == PETSC_SUCCESS);

    if (m_preconditioner)
    {
      ierr = KSPSetOperators(m_ksp, a, *m_preconditioner);
      assert(ierr == PETSC_SUCCESS);
    }
    else
    {
      ierr = KSPSetOperators(m_ksp, a, a);
      assert(ierr == PETSC_SUCCESS);
    }

    ierr = KSPSetFromOptions(m_ksp);
    assert(ierr == PETSC_SUCCESS);

    ierr = KSPSolve(m_ksp, b, x);
    assert(ierr == PETSC_SUCCESS);
  }

  KSP& KSP::setType(::KSPType type) noexcept
  {
    m_type = type;
    return *this;
  }

  KSP& KSP::setTolerances(
      PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt  maxIt) noexcept
  {
    m_rtol = rtol;
    m_abstol = abstol;
    m_dtol = dtol;
    m_maxIt = maxIt;
    return *this;
  }

  KSP& KSP::setPreconditioner(OperatorType preconditioner) noexcept
  {
    m_preconditioner = preconditioner;
    return *this;
  }
}

