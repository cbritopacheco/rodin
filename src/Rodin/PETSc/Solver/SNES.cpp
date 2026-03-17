#include <cassert>
#include <petsc.h>
#include <petscsnes.h>

#include "Rodin/Variational/Problem.h"
#include "Rodin/PETSc/Variational/Problem.h"

#include "SNES.h"

namespace Rodin::Solver
{
  SNES::SNES(KSP& ksp)
    : NewtonSolverParent(ksp),
      m_snes(PETSC_NULLPTR),
      m_type(SNESNEWTONLS),
      m_abstol(PETSC_DECIDE),
      m_rtol(PETSC_DECIDE),
      m_stol(PETSC_DECIDE),
      m_maxIt(PETSC_DECIDE),
      m_maxF(PETSC_DECIDE)
  {
    auto& problem = ksp.getProblem();
    auto& system = problem.getLinearSystem();
    const auto& comm = system.getCommunicator();

    PetscErrorCode ierr = SNESCreate(comm, &m_snes);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetFunction(
      m_snes,
      PETSC_NULLPTR,
      &SNES::Residual,
      this);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetJacobian(
      m_snes,
      system.getOperator(),
      system.getOperator(),
      &SNES::Jacobian,
      this);
    assert(ierr == PETSC_SUCCESS);

    // Use the user-provided KSP handle as the SNES sub-solver
    ierr = SNESSetKSP(m_snes, ksp.getHandle());
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

  PetscErrorCode SNES::Residual(::SNES, ::Vec x, ::Vec f, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);

    auto& system = self->getProblem().getLinearSystem();
    system.getSolution() = x;
    self->getProblem().assemble();

    PetscErrorCode ierr = VecCopy(system.getVector(), f);
    assert(ierr == PETSC_SUCCESS);
    if (ierr != PETSC_SUCCESS)
      return ierr;
    ierr = VecScale(f, -1.0);
    assert(ierr == PETSC_SUCCESS);
    return ierr;
  }

  PetscErrorCode SNES::Jacobian(::SNES, ::Vec x, ::Mat J, ::Mat P, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);

    auto& system = self->getProblem().getLinearSystem();
    system.getSolution() = x;
    self->getProblem().assemble();

    const auto& assembledJ = system.getOperator();
    PetscErrorCode ierr = PETSC_SUCCESS;
    // PETSc may pass callback workspace matrices distinct from the problem-managed
    // operator; copy assembled data when handles differ.
    if (J != assembledJ)
    {
      ierr = MatCopy(assembledJ, J, DIFFERENT_NONZERO_PATTERN);
      assert(ierr == PETSC_SUCCESS);
    }

    // Keep preconditioner matrix synchronized as well when PETSc uses a
    // separate handle for P.
    if (P && P != J && P != assembledJ)
    {
      ierr = MatCopy(assembledJ, P, DIFFERENT_NONZERO_PATTERN);
      assert(ierr == PETSC_SUCCESS);
    }
    return ierr;
  }

  void SNES::solve(VectorType& x)
  {
    auto& pb = this->getProblem();
    auto& system = pb.getLinearSystem();
    auto& b = system.getVector();
    auto& x0 = system.getSolution();

    PetscErrorCode ierr;

    ierr = SNESSetType(m_snes, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetTolerances(m_snes, m_abstol, m_rtol, m_stol, m_maxIt, m_maxF);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetInitialFunction(m_snes, x0);

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
