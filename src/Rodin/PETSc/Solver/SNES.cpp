#include <cassert>
#include <petscmat.h>
#include <petscvec.h>
#include <utility>
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
      m_abstol(DEFAULT_ABSTOL),
      m_rtol(DEFAULT_RTOL),
      m_stol(DEFAULT_STOL),
      m_maxIt(DEFAULT_MAXIT),
      m_maxF(DEFAULT_MAXF)
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

  SNES& SNES::setStateUpdate(StateUpdate update)
  {
    m_update = std::move(update);
    return *this;
  }

  void SNES::solve()
  {
    auto& problem = this->getProblem();
    auto& system = problem.getLinearSystem();
    auto& x = system.getSolution();
    solve(x);
  }

  PetscInt SNES::getIterationNumber() const
  {
    PetscInt its = 0;
    PetscErrorCode ierr = SNESGetIterationNumber(m_snes, &its);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return its;
  }

  ::SNESConvergedReason SNES::getConvergedReason() const
  {
    ::SNESConvergedReason reason;
    PetscErrorCode ierr = SNESGetConvergedReason(m_snes, &reason);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return reason;
  }

  bool SNES::converged() const
  {
    ::SNESConvergedReason reason;
    PetscErrorCode ierr = SNESGetConvergedReason(m_snes, &reason);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return reason > 0;
  }

  PetscErrorCode SNES::Update(::Vec x, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);
    PetscObjectState state;
    PetscErrorCode ierr = PetscObjectStateGet(reinterpret_cast<PetscObject>(x), &state);
    if (ierr)
      return ierr;
    if (self->m_updated && *self->m_updated == state)
      return PETSC_SUCCESS;
    if (self->m_update)
      self->m_update(x);
    self->m_updated = state;
    self->m_assembled.reset();
    return PETSC_SUCCESS;
  }

  PetscErrorCode SNES::Assemble(::Vec x, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);
    PetscErrorCode ierr = Update(x, ctx);
    if (ierr)
      return ierr;
    PetscObjectState state;
    ierr = PetscObjectStateGet(reinterpret_cast<PetscObject>(x), &state);
    if (ierr)
      return ierr;
    if (self->m_assembled && *self->m_assembled == state)
      return PETSC_SUCCESS;
    auto& problem = self->getProblem();
    problem.assemble();
    self->m_assembled = state;
    return PETSC_SUCCESS;
  }

  PetscErrorCode SNES::Residual(::SNES, ::Vec x, ::Vec f, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);
    PetscErrorCode ierr = Assemble(x, ctx);
    if (ierr)
      return ierr;
    auto& problem = self->getProblem();
    auto& system = problem.getLinearSystem();
    auto& b = system.getVector();
    ierr = VecCopy(b, f);
    if (ierr)
      return ierr;
    ierr = VecScale(f, -1.0);
    return ierr;
  }

  PetscErrorCode SNES::Jacobian(::SNES, ::Vec x, ::Mat J, ::Mat P, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);

    PetscErrorCode ierr = Assemble(x, ctx);
    if (ierr)
      return ierr;

    auto& problem = self->getProblem();
    auto& system = problem.getLinearSystem();
    const auto& A = system.getOperator();

    if (J != A)
    {
      ierr = MatCopy(A, J, DIFFERENT_NONZERO_PATTERN);
      if (ierr)
        return ierr;
    }

    if (P && P != J && P != A)
    {
      ierr = MatCopy(A, P, DIFFERENT_NONZERO_PATTERN);
      if (ierr)
        return ierr;
    }

    return PETSC_SUCCESS;
  }

  void SNES::solve(VectorType& x)
  {
    PetscErrorCode ierr;

    m_updated.reset();
    m_assembled.reset();

    ierr = SNESSetType(m_snes, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetTolerances(m_snes, m_abstol, m_rtol, m_stol, m_maxIt, m_maxF);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetFromOptions(m_snes);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSolve(m_snes, PETSC_NULLPTR, x);
    assert(ierr == PETSC_SUCCESS);

    ierr = Update(x, this);
    assert(ierr == PETSC_SUCCESS);

    // State is synchronized, but no assembled system should be trusted after solve.
    m_assembled.reset();
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
