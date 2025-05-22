#ifndef RODIN_ASSEMBLY_DEFAULT_H
#define RODIN_ASSEMBLY_DEFAULT_H

#include "Rodin/Context/Sequential.h"

#include "ForwardDecls.h"

#include "Sequential.h"
#include "Multithreaded.h"

namespace Rodin::Assembly
{
  template <>
  class Default<Context::Local>
  {
    public:
#ifdef RODIN_MULTITHREADED
      template <class LinearAlgebraType, class Object>
      using Type = Multithreaded<LinearAlgebraType, Object>;
#else
      template <class LinearAlgebraType, class Object>
      using Type = Sequential<LinearAlgebraType, Object>;
#endif
  };

  template <>
  class Default<Context::Local, Context::Local>
  {
    public:
#ifdef RODIN_MULTITHREADED
      template <class LinearAlgebraType, class Object>
      using Type = Multithreaded<LinearAlgebraType, Object>;
#else
      template <class LinearAlgebraType, class Object>
      using Type = Sequential<LinearAlgebraType, Object>;
#endif
 };

}

#endif
