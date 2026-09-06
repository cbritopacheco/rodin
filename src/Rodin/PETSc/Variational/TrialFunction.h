/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H

/**
 * @file TrialFunction.h
 * @brief PETSc-aware trial function wrappers.
 *
 * Provides the PETSc-aware trial-function wrapper class, a thin
 * wrapper around the core variational trial function that lives in the
 * PETSc namespace.  Its CTAD deduction guide automatically associates the
 * trial function with a PETSc-backed
 * PETSc-backed grid function as its solution type.
 *
 * @see <a href="class_rodin_1_1_p_e_t_sc_1_1_variational_1_1_test_function.html">Rodin::PETSc::Variational::TestFunction</a>
 * @see <a href="class_rodin_1_1_p_e_t_sc_1_1_variational_1_1_grid_function.html">Rodin::PETSc::Variational::GridFunction</a>
 * @see <a href="class_rodin_1_1_variational_1_1_bilinear_form_3_01_solution_00_01_trial_f_e_s_00_01_test_f_e_s_00_01_1_1_mat_01_4.html">Rodin::PETSc::Variational::BilinearForm</a>
 */

#include <petsc.h>

#include "Rodin/Variational/TrialFunction.h"

#include "GridFunction.h"

namespace Rodin::PETSc::Variational
{
  /**
   * @brief PETSc-aware trial function wrapper.
   *
   * Inherits all functionality from `Rodin::Variational::TrialFunction`
   * and triggers PETSc-specific CTAD.  The CTAD deduction guide
   * automatically sets the `Solution` template parameter to
   * `PETSc::Variational::GridFunction<FES>`, ensuring that the discrete
   * solution is stored in a PETSc @c Vec.
   *
   * @tparam Solution Solution type (typically
   *         `PETSc::Variational::GridFunction<FES>`).
   * @tparam FES      Finite element space type.
   *
   * @see <a href="_variational_2_trial_function_8h.html">Rodin::Variational::TrialFunction</a>
   */
  template <class Solution, class FES>
  class TrialFunction : public Rodin::Variational::TrialFunction<Solution, FES>
  {
    public:
      /// @brief Finite element space type (e.g. `P1<Real, Mesh<Context::Local>>`).
      using FESType =
        FES;

      /// @brief Parent class type (`Rodin::Variational::TrialFunction`).
      using Parent =
        Rodin::Variational::TrialFunction<Solution, FESType>;

      using Parent::Parent;
  };

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for PETSc::Variational::TrialFunction.
   */
  template <class FES>
  TrialFunction(const FES& fes) -> TrialFunction<GridFunction<FES>, FES>;
}

namespace Rodin::Variational
{
  /// @brief Marks PETSc trial functions as valid trial functions.
  template <class Solution, class FES>
  struct IsTrialFunction<PETSc::Variational::TrialFunction<Solution, FES>>
  {
    /// @brief True because the PETSc wrapper models a trial function.
      static constexpr Boolean Value = true;
  };
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Form-language traits specialization for PETSc trial functions.
   *
   * Provides `FESType` and `SolutionType` so that the Rodin form language
   * can determine the finite element space and solution grid function
   * associated with a PETSc trial function.
   *
   * @tparam Solution Solution grid function type.
   * @tparam FES      Finite element space type.
   */
  template <class Solution, class FES>
  struct Traits<PETSc::Variational::TrialFunction<Solution, FES>>
  {
    /// @brief Finite element space type.
    using FESType = FES;
    /// @brief Solution type.
    using SolutionType = Solution;
  };

}


#endif
