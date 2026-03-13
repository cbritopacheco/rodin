#include <cassert>
#include <petsc.h>
#include <utility>

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

  SNES& SNES::setFunction(FunctionType f, VectorType residual)
  {
    assert(static_cast<bool>(f));
    if (!f)
      return *this;
    m_function = std::move(f);
    m_residual = residual;
    PetscErrorCode ierr = SNESSetFunction(m_snes, m_residual, &SNES::FunctionCallback, this);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return *this;
  }

  SNES& SNES::setJacobian(
      JacobianType j,
      MatrixType jacobian,
      MatrixType preconditioner)
  {
    assert(static_cast<bool>(j));
    if (!j)
      return *this;
    m_jacobian = std::move(j);
    m_jacobianOperator = jacobian;
    m_preconditionerOperator = preconditioner;
    PetscErrorCode ierr = SNESSetJacobian(
        m_snes,
        m_jacobianOperator,
        m_preconditionerOperator,
        &SNES::JacobianCallback,
        this);
    assert(ierr == PETSC_SUCCESS);
    (void) ierr;
    return *this;
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

  ::SNES& SNES::getHandle() noexcept
  {
    return m_snes;
  }

  const ::SNES& SNES::getHandle() const noexcept
  {
    return m_snes;
  }

  PetscErrorCode SNES::FunctionCallback(
      ::SNES snes, ::Vec x, ::Vec f, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    if (!self)
      return PETSC_ERR_ARG_NULL;
    if (!self->m_function)
      return PETSC_ERR_ARG_WRONGSTATE;
    return self->m_function(snes, x, f);
  }

  PetscErrorCode SNES::JacobianCallback(
      ::SNES snes, ::Vec x, ::Mat A, ::Mat P, void* ctx)
  {
    auto* self = static_cast<SNES*>(ctx);
    if (!self)
      return PETSC_ERR_ARG_NULL;
    if (!self->m_jacobian)
      return PETSC_ERR_ARG_WRONGSTATE;
    return self->m_jacobian(snes, x, A, P);
  }
}
