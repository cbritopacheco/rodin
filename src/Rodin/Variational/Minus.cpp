#include "Minus.h"

#include "Sum.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
   Sum<FunctionBase, FunctionBase>
   operator-(const FunctionBase& lhs, const FunctionBase& rhs)
   {
      return Sum(lhs, UnaryMinus(rhs));
   }
}
