#ifndef RODIN_VARIATIONAL_MINUS_H
#define RODIN_VARIATIONAL_MINUS_H

#include <type_traits>

#include "ForwardDecls.h"
#include "Sum.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
  template <class LHSDerived, class Number, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator-(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Sum(lhs, UnaryMinus(ScalarFunction(rhs)));
  }

  template <class Number, class RHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator-(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(UnaryMinus(ScalarFunction(lhs)), rhs);
  }

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator-(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Sum(lhs, UnaryMinus(rhs));
  }

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator-(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs, const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
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
