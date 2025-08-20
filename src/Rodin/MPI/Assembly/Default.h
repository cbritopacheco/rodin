#ifndef RODIN_MPI_ASSEMBLY_DEFAULT_H
#define RODIN_MPI_ASSEMBLY_DEFAULT_H

#include "Rodin/MPI/Context.h"
#include "Rodin/Assembly/Default.h"

#include "MPI.h"

namespace Rodin::Assembly
{
  template <>
  class Default<Context::MPI>
  {
    public:
      template <class LinearAlgebraType, class Object>
      using Type = MPI<LinearAlgebraType, Object>;
  };

  template <>
  class Default<Context::MPI, Context::MPI>
  {
    public:
      template <class LinearAlgebraType, class Object>
      using Type = MPI<LinearAlgebraType, Object>;
  };
}

#endif
