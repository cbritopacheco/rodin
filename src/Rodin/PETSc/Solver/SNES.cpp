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
      m_functionAssembly(),
      m_jacobianAssembly(),
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
    m_functionAssembly =
      [this](VectorType& residual, const VectorType& x)
      {
        if (!m_functionCallback)
          return;
        PetscErrorCode ierr = m_functionCallback(
            m_snes,
            const_cast<VectorType>(x),
            residual,
            m_functionContext);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      };
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
    m_jacobianAssembly =
      [this](MatrixType& jacobian, const VectorType& x)
      {
        if (!m_jacobianCallback)
          return;
        const MatrixType preconditioner = m_preconditionerOperator
          ? m_preconditionerOperator
          : jacobian;
        PetscErrorCode ierr = m_jacobianCallback(
            m_snes,
            const_cast<VectorType>(x),
            jacobian,
            preconditioner,
            m_jacobianContext);
        assert(ierr == PETSC_SUCCESS);
        (void) ierr;
      };
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

  void SNES::solve(VectorType& x)
  {
    solve(PETSC_NULLPTR, x);
  }

  void SNES::solve(LinearSystemType& system)
  {
    solve(system.getVector(), system.getSolution());
  }

  const SNES::ResidualAssembly& SNES::getFunction() const
  {
    return m_functionAssembly;
  }

  const SNES::JacobianAssembly& SNES::getJacobian() const
  {
    return m_jacobianAssembly;
  }

  const SNES::KSPType& SNES::getLinearSolver() const
  {
    return m_kspHandle;
  }

  SNES::KSPType& SNES::getLinearSolver()
  {
    return m_kspHandle;
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
