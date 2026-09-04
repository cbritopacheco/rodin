/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H

/**
 * @file TestFunction.h
 * @brief PETSc-aware test function wrappers.
 *
 * Provides the PETSc-aware test-function wrapper class, a thin
 * wrapper around the core variational test function that lives in the
 * PETSc namespace and triggers PETSc-specific template argument deduction
 * for linear forms, bilinear forms, and problems.
 *
 * @see <a href="class_rodin_1_1_p_e_t_sc_1_1_variational_1_1_trial_function.html">Rodin::PETSc::Variational::TrialFunction</a>
 * @see <a href="class_rodin_1_1_variational_1_1_linear_form_3_01_f_e_s_00_01_1_1_vec_01_4.html">Rodin::PETSc::Variational::LinearForm</a>
 * @see <a href="class_rodin_1_1_variational_1_1_bilinear_form_3_01_solution_00_01_trial_f_e_s_00_01_test_f_e_s_00_01_1_1_mat_01_4.html">Rodin::PETSc::Variational::BilinearForm</a>
 */

#include <petsc.h>

#include "Rodin/Variational/TestFunction.h"

namespace Rodin::PETSc::Variational
{
  /**
   * @brief PETSc-aware test function wrapper.
   *
   * Inherits all functionality from `Rodin::Variational::TestFunction`
   * and triggers PETSc-specific CTAD so that constructing a
   * `PETSc::Variational::TestFunction` from a finite element space
   * automatically selects the PETSc-backed linear form, bilinear form,
   * and problem specializations.
   *
   * @tparam FES Finite element space type (e.g. `P1<Real, Mesh<Context::Local>>`).
   *
   * @see <a href="_variational_2_test_function_8h.html">Rodin::Variational::TestFunction</a>
   */
  template <class FES>
  class TestFunction : public Rodin::Variational::TestFunction<FES>
  {
    public:
      /// @brief Finite element space type (e.g. `P1<Real, Mesh<Context::Local>>`).
      using FESType =
        FES;

      /// @brief Parent class type (`Rodin::Variational::TestFunction`).
      using Parent =
        Rodin::Variational::TestFunction<FESType>;

      using Parent::Parent;
  };

  /**
   * @ingroup RodinCTAD
   * @brief Deduction guide for PETSc::Variational::TestFunction.
   */
  template <class FES>
  TestFunction(const FES& fes) -> TestFunction<FES>;
}

namespace Rodin::Variational
{
  /// @brief Marks PETSc test functions as valid test functions.
  template <class FES>
  struct IsTestFunction<PETSc::Variational::TestFunction<FES>>
  {
    /// @brief True because the PETSc wrapper models a test function.
      static constexpr Boolean Value = true;
  };
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Form-language traits specialization for PETSc test functions.
   *
   * Provides the `FESType` alias so that the Rodin form language can
   * determine the finite element space associated with a PETSc test
   * function during template argument deduction.
   *
   * @tparam FES Finite element space type.
   */
  template <class FES>
  struct Traits<PETSc::Variational::TestFunction<FES>>
  {
    /// @brief Finite element space type.
    using FESType = FES;
  };
}

#endif
