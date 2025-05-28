/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cassert>
#include <petsc.h>

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/PETSc/FormLanguage/Traits.h"

#include "Rodin/Variational/Problem.h"

#include "KSP.h"

namespace Rodin::Solver
{
  KSP::KSP(ProblemType& pb)
    : Parent(pb),
      m_ksp(PETSC_NULLPTR),
      m_type(PETSC_NULLPTR),
      m_rtol(PETSC_DECIDE),
      m_abstol(PETSC_DECIDE),
      m_dtol(PETSC_DECIDE),
      m_maxIt(PETSC_DECIDE),
      m_preconditioner(PETSC_NULLPTR)
  {
    PetscErrorCode ierr;
    MPI_Comm comm;
    ierr = PetscObjectGetComm(
        reinterpret_cast<PetscObject>(pb.getLinearSystem().getOperator()), &comm);
    ierr = KSPCreate(comm, &m_ksp);
    assert(ierr == PETSC_SUCCESS);
  }

  KSP::~KSP()
  {
    if (m_ksp)
    {
      PetscErrorCode ierr = KSPDestroy(&m_ksp);
      assert(ierr == PETSC_SUCCESS);
      m_ksp = nullptr;
    }
  }

  ::PetscObject& KSP::getHandle() noexcept
  {
    return reinterpret_cast<::PetscObject&>(m_ksp);
  }

  void KSP::solve(OperatorType& A, VectorType& x, VectorType& b)
  {
    PetscErrorCode ierr;

    if (x)
    {
      // use nonzero initial guess
      ierr = KSPSetInitialGuessNonzero(m_ksp, PETSC_TRUE);
      assert(ierr == PETSC_SUCCESS);
    }
    else
    {
      ierr = VecDuplicate(b, &x);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecZeroEntries(x);
      assert(ierr == PETSC_SUCCESS);
    }

    assert(x);

    // configure solver
    ierr = KSPSetType(m_ksp, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = KSPSetTolerances(m_ksp, m_rtol, m_abstol, m_dtol, m_maxIt);
    assert(ierr == PETSC_SUCCESS);

    // set operators (use A as both if P not set)
    if (m_preconditioner)
      ierr = KSPSetOperators(m_ksp, A, m_preconditioner);
    else
      ierr = KSPSetOperators(m_ksp, A, A);
    assert(ierr == PETSC_SUCCESS);

    // allow CLI overrides
    ierr = KSPSetFromOptions(m_ksp);
    assert(ierr == PETSC_SUCCESS);

    // solve
    ierr = KSPSolve(m_ksp, b, x);
    assert(ierr == PETSC_SUCCESS);
  }

  KSP& KSP::setType(::KSPType type) noexcept
  {
    m_type = type;
    return *this;
  }

  KSP& KSP::setTolerances(PetscReal rtol,
                          PetscReal abstol,
                          PetscReal dtol,
                          PetscInt  maxIt) noexcept
  {
    m_rtol   = rtol;
    m_abstol = abstol;
    m_dtol   = dtol;
    m_maxIt  = maxIt;
    return *this;
  }

  KSP& KSP::setPreconditioner(OperatorType preconditioner) noexcept
  {
    m_preconditioner = preconditioner;
    return *this;
  }
}

