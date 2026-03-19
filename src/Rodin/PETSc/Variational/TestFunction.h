#ifndef RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TESTFUNCTION_H

/**
 * @file
 * @brief PETSc-aware test function wrappers.
 */

#include <petsc.h>

#include "Rodin/Variational/TestFunction.h"

namespace Rodin::PETSc::Variational
{
  /**
   * @brief PETSc-aware test function wrapper.
   *
   * Inherits all functionality from @ref Rodin::Variational::TestFunction
   * and triggers PETSc-specific template argument deduction.
   *
   * @tparam FES Finite element space type.
   */
  template <class FES>
  class TestFunction : public Rodin::Variational::TestFunction<FES>
  {
    public:
      /// @brief Finite element space type.
      using FESType =
        FES;

      /// @brief Parent class type.
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
   * @brief Form-language traits for PETSc test functions.
   */
  template <class FES>
  struct Traits<PETSc::Variational::TestFunction<FES>>
  {
    /// @brief Finite element space type.
    using FESType = FES;
  };
}

#endif
