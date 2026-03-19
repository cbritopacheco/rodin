#ifndef RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H

/**
 * @file
 * @brief PETSc-aware trial function wrappers.
 */

#include <petsc.h>

#include "Rodin/Variational/TrialFunction.h"

#include "GridFunction.h"

namespace Rodin::PETSc::Variational
{
  /**
   * @brief PETSc-aware trial function wrapper.
   *
   * Inherits all functionality from @ref Rodin::Variational::TrialFunction
   * and triggers PETSc-specific template argument deduction.
   *
   * @tparam Solution Solution type (typically a PETSc GridFunction).
   * @tparam FES      Finite element space type.
   */
  template <class Solution, class FES>
  class TrialFunction : public Rodin::Variational::TrialFunction<Solution, FES>
  {
    public:
      /// @brief Finite element space type.
      using FESType =
        FES;

      /// @brief Parent class type.
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
    static constexpr Boolean Value = true;
  };
}

namespace Rodin::FormLanguage
{
  /**
   * @brief Form-language traits for PETSc trial functions.
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
