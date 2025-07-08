#ifndef RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H

#include <petsc.h>

#include "Rodin/Variational/TrialFunction.h"

#include "GridFunction.h"

namespace Rodin::PETSc
{
  template <class FES>
  class TrialFunction : public Variational::TrialFunction<GridFunction<FES>, FES>
  {
    public:
      using Parent = Variational::TrialFunction<GridFunction<FES>, FES>;
      using FESType = FES;
      using SolutionType = GridFunction<FES>;

      TrialFunction(const FESType& fes)
        : Parent(fes)
      {}

      TrialFunction(SolutionType& gf, const FESType& fes)
        : Parent(gf, fes)
      {}

      TrialFunction(SolutionType&& gf, const FESType& fes)
        : Parent(std::move(gf), fes)
      {}
  };

  template <class FES>
  TrialFunction(const FES& fes) -> TrialFunction<FES>;

  template <class FES>
  TrialFunction(GridFunction<FES>& gf, const FES& fes) -> TrialFunction<FES>;

  template <class FES>
  TrialFunction(GridFunction<FES>&& gf, const FES& fes) -> TrialFunction<FES>;
}

namespace Rodin::FormLanguage
{
  template <class FES>
  struct Traits<PETSc::TrialFunction<FES>>
  {
    using FESType = FES;
    using SolutionType = typename PETSc::GridFunction<FES>;
  };
}

namespace Rodin::Variational
{
  template <class FES>
  struct IsTrialFunction<PETSc::TrialFunction<FES>>
  {
    static constexpr Boolean Value = true;
  };
}


#endif
