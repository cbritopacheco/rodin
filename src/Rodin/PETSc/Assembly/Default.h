#ifndef RODIN_PETSC_ASSEMBLY_DEFAULT_H
#define RODIN_PETSC_ASSEMBLY_DEFAULT_H

#include "Rodin/Assembly/Default.h"

namespace Rodin::PETSc::Assembly
{
  template <class ... Ts>
  using Default = Rodin::Assembly::Default<Ts...>;
}

#endif
