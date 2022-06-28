#include "Minus.h"

#include "Sum.h"
#include "UnaryMinus.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   Sum<ScalarFunctionBase, ScalarFunctionBase>
   operator-(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Sum<ScalarFunctionBase, ScalarFunctionBase>(lhs, UnaryMinus(rhs));
   }

   Sum<VectorFunctionBase, VectorFunctionBase>
   operator-(const VectorFunctionBase& lhs, const VectorFunctionBase& rhs)
   {
      return Sum<VectorFunctionBase, VectorFunctionBase>(lhs, UnaryMinus(rhs));
   }
}
