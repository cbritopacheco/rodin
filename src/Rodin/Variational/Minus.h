#ifndef RODIN_VARIATIONAL_MINUS_H
#define RODIN_VARIATIONAL_MINUS_H

#include <type_traits>

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Sum<ScalarFunctionBase, ScalarFunctionBase>>
   operator-(const ScalarFunctionBase& lhs, T v)
   {
      return lhs - ScalarFunction(v);
   }

   Sum<ScalarFunctionBase, ScalarFunctionBase>
   operator-(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs);

   Sum<VectorFunctionBase, VectorFunctionBase>
   operator-(const VectorFunctionBase& lhs, const VectorFunctionBase& rhs);

   template <ShapeFunctionSpaceType Space>
   Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>
   operator-(const ShapeFunctionBase<Space>& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Sum(lhs, UnaryMinus(rhs));
   }
}

#endif
