#ifndef RODIN_VARIATIONAL_MINUS_H
#define RODIN_VARIATIONAL_MINUS_H

#include <type_traits>

#include "ForwardDecls.h"
#include "Sum.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Sum<FunctionBase, FunctionBase>>
   operator-(const FunctionBase& lhs, T v)
   {
      return Sum(lhs, UnaryMinus(ScalarFunction(v)));
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Sum<FunctionBase, FunctionBase>>
   operator-(T v, const FunctionBase& rhs)
   {
      return Sum(ScalarFunction(v), UnaryMinus(rhs));
   }

   Sum<FunctionBase, FunctionBase>
   operator-(const FunctionBase& lhs, const FunctionBase& rhs);

   template <ShapeFunctionSpaceType Space>
   Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>
   operator-(const ShapeFunctionBase<Space>& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Sum(lhs, UnaryMinus(rhs));
   }

   FormLanguage::List<BilinearFormIntegratorBase>
   operator-(
         const BilinearFormIntegratorBase& lhs,
         const BilinearFormIntegratorBase& rhs);

   FormLanguage::List<BilinearFormIntegratorBase>
   operator-(
         const BilinearFormIntegratorBase& lhs,
         const FormLanguage::List<BilinearFormIntegratorBase>& rhs);

   FormLanguage::List<BilinearFormIntegratorBase>
   operator-(
         const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
         const BilinearFormIntegratorBase& rhs);

   FormLanguage::List<BilinearFormIntegratorBase>
   operator-(
         const FormLanguage::List<BilinearFormIntegratorBase>& lhs,
         const FormLanguage::List<BilinearFormIntegratorBase>& rhs);
}

#endif
