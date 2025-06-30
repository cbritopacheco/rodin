#ifndef RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H
#define RODIN_PETSC_VARIATIONAL_GRIDFUNCTION_H

#include <petsc.h>
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Variational
{
  template <class FES>
  class GridFunction<FES, ::Vec> final
    : public GridFunctionBase<GridFunction<FES, ::Vec>>
  {
  };
}

#endif
