/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_SOLVER_GMRES_H
#define RODIN_PETSC_SOLVER_GMRES_H

/**
 * @file GMRES.h
 * @brief PETSc specialization of the GMRES solver.
 *
 * Provides the `GMRES<PETSc::Math::LinearSystem>` specialization that
 * sets the PETSc KSP type to `KSPGMRES`.  GMRES is applicable to
 * general (non-symmetric) linear systems
 * @f$ A\mathbf{x} = \mathbf{b} @f$.
 *
 * @see <a href="class_rodin_1_1_solver_1_1_k_s_p.html">Rodin::Solver::KSP</a>
 * @see <a href="class_rodin_1_1_solver_1_1_c_g_3_01_p_e_t_sc_1_1_math_1_1_linear_system_01_4.html">Rodin::PETSc::Solver::CG</a>
 */

#include <petscksp.h>

#include "Rodin/PETSc/Math/LinearSystem.h"
#include "Rodin/Variational/ForwardDecls.h"

#include "KSP.h"

namespace Rodin::Solver
{
  /**
   * @ingroup GMRESSpecializations
   * @brief Generalized minimal residual (GMRES) solver for PETSc linear
   *        systems.
   *
   * Specialization of @ref Rodin::Solver::GMRES that targets
   * @ref Rodin::PETSc::Math::LinearSystem and sets the underlying KSP
   * algorithm to `KSPGMRES`.  Suitable for non-symmetric and
   * indefinite systems.
   *
   * @see <a href="class_rodin_1_1_solver_1_1_k_s_p.html">Rodin::Solver::KSP</a>
   * @see <a href="class_rodin_1_1_solver_1_1_g_m_r_e_s_3_01_p_e_t_sc_1_1_math_1_1_linear_system_01_4.html">Rodin::PETSc::Solver::GMRES</a>
   */
  template <>
  class GMRES<PETSc::Math::LinearSystem> final : public KSP
  {
    public:
      /// @brief PETSc matrix type (@c Mat) for the system operator.
      using OperatorType = ::Mat;
      /// @brief PETSc vector type (@c Vec) for the RHS and solution.
      using VectorType = ::Vec;
      /// @brief Scalar type (`PetscScalar`).
      using ScalarType = PetscScalar;
      /// @brief Linear system type coupling @f$ A @f$, @f$ \mathbf{b} @f$, and @f$ \mathbf{x} @f$.
      using LinearSystemType = PETSc::Math::LinearSystem;
      /// @brief Base problem type that provides the linear system.
      using ProblemBaseType = Variational::ProblemBase<LinearSystemType>;
      /// @brief Parent class type (@ref Rodin::Solver::KSP).
      using Parent = KSP;
      using Parent::solve;

      /**
       * @brief Construct a GMRES solver and set the PETSc type to KSPGMRES.
       * @param pb The variational problem to solve.
       */
      explicit GMRES(ProblemBaseType& pb);

      /**
       * @brief Copy constructor.
       * @param other Another GMRES instance.
       */
      GMRES(const GMRES& other);

      /**
       * @brief Move constructor.
       * @param other Moved-from GMRES instance.
       */
      GMRES(GMRES&& other);

      /// @brief Creates a heap-allocated copy of this GMRES solver.
      /// @returns Pointer to the cloned GMRES instance.
      GMRES* copy() const noexcept override
      {
        return new GMRES(*this);
      }
  };
}

namespace Rodin::PETSc::Solver
{
  /**
   * @brief Convenient PETSc alias for
   * @ref Rodin::Solver::GMRES<Rodin::PETSc::Math::LinearSystem>.
   */
  using GMRES = Rodin::Solver::GMRES<Math::LinearSystem>;
}

#endif // RODIN_PETSC_SOLVER_GMRES_H
