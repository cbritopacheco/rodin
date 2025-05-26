#ifndef RODIN_ASSEMBLY_DEFAULT_H
#define RODIN_ASSEMBLY_DEFAULT_H

#include "Rodin/Context/Local.h"

#include "ForwardDecls.h"

#ifdef RODIN_USE_OPENMP

#include "OpenMP.h"

namespace Rodin::Assembly
{
  template <>
  class Default<Context::Local>
  {
    public:
      template <class LinearAlgebraType, class Object>
      using Type = OpenMP<LinearAlgebraType, Object>;
  };

  template <>
  class Default<Context::Local, Context::Local>
  {
    public:
      template <class LinearAlgebraType, class Object>
      using Type = OpenMP<LinearAlgebraType, Object>;
 };
}

#else

#include "Sequential.h"

namespace Rodin::Assembly
{
  template <>
  class Default<Context::Local>
  {
    public:
      template <class LinearAlgebraType, class Object>
      using Type = Sequential<LinearAlgebraType, Object>;
  };

  template <>
  class Default<Context::Local, Context::Local>
  {
    public:
      template <class LinearAlgebraType, class Object>
      using Type = Sequential<LinearAlgebraType, Object>;
 };
}

#endif

#endif
