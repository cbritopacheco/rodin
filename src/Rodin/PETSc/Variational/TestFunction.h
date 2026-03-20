#ifndef RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H

/**
 * @file TestFunction.h
 * @brief PETSc-aware test function wrappers.
 *
 * Provides the @ref Rodin::PETSc::Variational::TestFunction class, a thin
 * wrapper around @ref Rodin::Variational::TestFunction that lives in the
 * PETSc namespace and triggers PETSc-specific template argument deduction
 * for linear forms, bilinear forms, and problems.
 *
 * @see Rodin::PETSc::Variational::TrialFunction,
 *      Rodin::PETSc::Variational::LinearForm,
 *      Rodin::PETSc::Variational::BilinearForm
 */

#include <petsc.h>

#include "Rodin/Variational/TestFunction.h"

namespace Rodin::PETSc::Variational
{
  /**
   * @brief PETSc-aware test function wrapper.
   *
   * Inherits all functionality from @ref Rodin::Variational::TestFunction
   * and triggers PETSc-specific CTAD so that constructing a
   * `PETSc::Variational::TestFunction` from a finite element space
   * automatically selects the PETSc-backed linear form, bilinear form,
   * and problem specializations.
   *
   * @tparam FES Finite element space type (e.g. `P1<Real, Mesh<Context::Local>>`).
   *
   * @see Rodin::Variational::TestFunction
   */
  template <class FES>
  class TestFunction : public Rodin::Variational::TestFunction<FES>
  {
    public:
      /// @brief Finite element space type (e.g. `P1<Real, Mesh<Context::Local>>`).
      using FESType =
        FES;

      /// @brief Parent class type (@ref Rodin::Variational::TestFunction).
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
