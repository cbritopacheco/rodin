#include "Product.h"

#include "../TrialFunction.h"
#include "../TestFunction.h"

namespace Rodin::Variational::FormLanguage
{
   Product<TrialFunctionBase, TestFunctionBase>
   operator*(const TrialFunctionBase& lhs, const TestFunctionBase& rhs)
   {
      return Product(lhs, rhs);
   }

   Product<ScalarCoefficientBase, TestFunctionBase>
   operator*(const ScalarCoefficientBase& lhs, const TestFunctionBase& rhs)
   {
      return Product(lhs, rhs);
   }
}
