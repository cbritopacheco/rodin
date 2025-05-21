#ifndef RODIN_ASSEMBLY_DEFAULT_H
#define RODIN_ASSEMBLY_DEFAULT_H

#include "ForwardDecls.h"
#include "Rodin/Context/Sequential.h"

namespace Rodin::Assembly
{
  template <class ... Ts>
  class Default;

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
