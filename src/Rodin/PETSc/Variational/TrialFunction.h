#ifndef RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_TRIALFUNCTION_H

#include <petsc.h>

#include "Rodin/Variational/TrialFunction.h"

#include "GridFunction.h"

namespace Rodin::PETSc::Variational
{
  template <class FES>
  using TrialFunction = Rodin::Variational::TrialFunction<PETSc::Variational::GridFunction<FES>, FES>;
}

namespace Rodin::Variational
{
  template <class FES>
  struct IsTrialFunction<PETSc::Variational::TrialFunction<FES>>
  {
    static constexpr Boolean Value = true;
  };
}


#endif
