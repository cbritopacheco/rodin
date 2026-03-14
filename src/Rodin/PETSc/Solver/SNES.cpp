#include <cassert>
#include <petsc.h>

#include "SNES.h"

namespace Rodin::Solver
{
  SNES::SNES(MPI_Comm comm)
    : m_snes(PETSC_NULLPTR),
      m_type(SNESNEWTONLS),
      m_abstol(PETSC_DECIDE),
      m_rtol(PETSC_DECIDE),
      m_stol(PETSC_DECIDE),
      m_maxIt(PETSC_DECIDE),
      m_maxF(PETSC_DECIDE),
      m_kspHandle(PETSC_NULLPTR),
      m_functionCallback(PETSC_NULLPTR),
      m_jacobianCallback(PETSC_NULLPTR),
      m_functionContext(PETSC_NULLPTR),
      m_jacobianContext(PETSC_NULLPTR),
      m_residual(PETSC_NULLPTR),
      m_jacobianOperator(PETSC_NULLPTR),
      m_preconditionerOperator(PETSC_NULLPTR)
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
    m_kspHandle = ksp;
    return *this;
  }

  SNES& SNES::setFunction(
      FunctionCallbackType f,
      void* ctx,
      VectorType residual)
  {
    if (!f)
      return *this;
    m_functionCallback = f;
    m_functionContext = ctx;
    m_residual = residual;
    PetscErrorCode ierr = SNESSetFunction(
        m_snes,
        m_residual,
        m_functionCallback,
        m_functionContext);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return *this;
  }

  SNES& SNES::setFunction(
      FunctionCallbackType f,
      VectorType residual)
  {
    return setFunction(f, PETSC_NULLPTR, residual);
  }

  SNES& SNES::setJacobian(
      JacobianCallbackType j,
      void* ctx,
      MatrixType jacobian,
      MatrixType preconditioner)
  {
    if (!j)
      return *this;
    m_jacobianCallback = j;
    m_jacobianContext = ctx;
    m_jacobianOperator = jacobian;
    m_preconditionerOperator = preconditioner;
    PetscErrorCode ierr = SNESSetJacobian(
        m_snes,
        m_jacobianOperator,
        m_preconditionerOperator,
        m_jacobianCallback,
        m_jacobianContext);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return *this;
  }

  SNES& SNES::setJacobian(
      JacobianCallbackType j,
      MatrixType jacobian,
      MatrixType preconditioner)
  {
    return setJacobian(j, PETSC_NULLPTR, jacobian, preconditioner);
  }

  void SNES::solve(VectorType b, VectorType x)
  {
    PetscErrorCode ierr;

    ierr = SNESSetType(m_snes, m_type);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSetTolerances(m_snes, m_abstol, m_rtol, m_stol, m_maxIt, m_maxF);
    assert(ierr == PETSC_SUCCESS);

    if (m_kspHandle)
    {
      ierr = SNESSetKSP(m_snes, m_kspHandle);
      assert(ierr == PETSC_SUCCESS);
    }

    ierr = SNESSetFromOptions(m_snes);
    assert(ierr == PETSC_SUCCESS);

    ierr = SNESSolve(m_snes, b, x);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
  }

  void SNES::solve(LinearSystemType& system)
  {
    solve(system.getVector(), system.getSolution());
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
