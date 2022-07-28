#include "MatrixFunction.h"
#include "Exceptions.h"

#include "Trace.h"

namespace Rodin::Variational
{
   Trace::Trace(const FunctionBase& m)
      : m_matrix(m.copy())
   {
      if (m.getRangeShape().height() != m.getRangeShape().width())
         NotSquareRangeShapeException(m.getRangeShape()).raise();
   }

   Trace::Trace(const Trace& other)
      :  ScalarFunctionBase(other),
         m_matrix(other.m_matrix->copy())
   {}

   Trace::Trace(Trace&& other)
      :  ScalarFunctionBase(other),
         m_matrix(std::move(other.m_matrix))
   {}

   double Trace::getValue(
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip) const
   {
      mfem::DenseMatrix mat;
      m_matrix->getValue(mat, trans, ip);
      return mat.Trace();
   }
}
