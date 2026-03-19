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
   * @brief PETSc namespace wrapper around @ref Rodin::Variational::TrialFunction.
   */
  template <class Solution, class FES>
  class TrialFunction : public Rodin::Variational::TrialFunction<Solution, FES>
  {
    public:
      using FESType =
        FES;

      using Parent =
        Rodin::Variational::TrialFunction<Solution, FESType>;

      using Parent::Parent;
  };

  template <class FES>
  TrialFunction(const FES& fes) -> TrialFunction<GridFunction<FES>, FES>;
}

namespace Rodin::Variational
{
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
    using FESType = FES;
    using SolutionType = Solution;
  };

}


#endif
