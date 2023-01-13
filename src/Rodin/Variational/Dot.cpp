#include "Exceptions.h"
#include "Dot.h"

namespace Rodin::Variational
{
   // ---- Dot<FunctionBase, FunctionBase> -----------------------------------
   Dot<FunctionBase, FunctionBase>::Dot(const FunctionBase& a, const FunctionBase& b)
      : m_a(a.copy()), m_b(b.copy())
   {
      if (a.getRangeShape() != b.getRangeShape())
         RangeShapeMismatchException(a.getRangeShape(), b.getRangeShape()).raise();
   }

   Dot<FunctionBase, FunctionBase>::Dot(const Dot& other)
      :  ScalarFunctionBase(other),
         m_a(other.m_a->copy()), m_b(other.m_b->copy())
   {}

   Dot<FunctionBase, FunctionBase>::Dot(Dot&& other)
      :  ScalarFunctionBase(std::move(other)),
         m_a(std::move(other.m_a)), m_b(std::move(other.m_b))
   {}

   Dot<FunctionBase, FunctionBase>&
   Dot<FunctionBase, FunctionBase>::traceOf(Geometry::Attribute attrs)
   {
      ScalarFunctionBase::traceOf(attrs);
      m_a->traceOf(attrs);
      m_b->traceOf(attrs);
      return *this;
   }
}
