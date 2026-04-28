#include <cassert>
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
      m_x(PETSC_NULLPTR),
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
    if (m_x)
    {
      PetscErrorCode ierr = VecDestroy(&m_x);
      assert(ierr == PETSC_SUCCESS);
      m_x = PETSC_NULLPTR;
      (void) ierr;
    }
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
    m_stateUpdate = std::move(update);
    return *this;
  }

  SNES& SNES::setSubVector(size_t offset, const PETSc::Math::Vector& src)
  {
    if (!m_x)
    {
      auto& system = this->getProblem().getLinearSystem();
      PetscErrorCode ierr = VecDuplicate(system.getSolution(), &m_x);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecZeroEntries(m_x);
      assert(ierr == PETSC_SUCCESS);
      (void) ierr;
    }

    PetscErrorCode ierr;

    PetscInt n = 0;
    ierr = VecGetSize(src, &n);
    assert(ierr == PETSC_SUCCESS);

    ::IS is = PETSC_NULLPTR;
    ierr = ISCreateStride(PETSC_COMM_SELF, n, static_cast<PetscInt>(offset), 1, &is);
    assert(ierr == PETSC_SUCCESS);

    ::Vec sub = PETSC_NULLPTR;
    ierr = VecGetSubVector(m_x, is, &sub);
    assert(ierr == PETSC_SUCCESS);

    ierr = VecCopy(src, sub);
    assert(ierr == PETSC_SUCCESS);

    ierr = VecRestoreSubVector(m_x, is, &sub);
    assert(ierr == PETSC_SUCCESS);

    ierr = ISDestroy(&is);
    assert(ierr == PETSC_SUCCESS);

    (void) ierr;
    return *this;
  }

  void SNES::solve()
  {
    if (!m_x)
    {
      auto& system = this->getProblem().getLinearSystem();
      PetscErrorCode ierr = VecDuplicate(system.getSolution(), &m_x);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecZeroEntries(m_x);
      assert(ierr == PETSC_SUCCESS);
      (void) ierr;
    }
    solve(m_x);
  }

  PetscInt SNES::getIterationNumber() const
  {
    PetscInt its = 0;
    PetscErrorCode ierr = SNESGetIterationNumber(m_snes, &its);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return its;
  }

  bool SNES::hasConverged() const
  {
    ::SNESConvergedReason reason;
    PetscErrorCode ierr = SNESGetConvergedReason(m_snes, &reason);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return reason > 0;
  }

  PetscErrorCode SNES::Residual(::SNES, ::Vec x, ::Vec f, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    assert(self);

    auto& system = self->getProblem().getLinearSystem();
    if (self->m_stateUpdate)
      self->m_stateUpdate(x);
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
    if (self->m_stateUpdate)
      self->m_stateUpdate(x);
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
    PetscErrorCode ierr;

    ierr = SNESSetType(m_snes, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetTolerances(m_snes, m_abstol, m_rtol, m_stol, m_maxIt, m_maxF);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetFromOptions(m_snes);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSolve(m_snes, PETSC_NULLPTR, x);
    assert(ierr == PETSC_SUCCESS);

    if (m_stateUpdate)
      m_stateUpdate(x);
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
