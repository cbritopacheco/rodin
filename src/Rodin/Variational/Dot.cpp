#include "MatrixCoefficient.h"

#include "Dot.h"

namespace Rodin::Variational
{
   Dot<MatrixCoefficientBase>::Dot(const MatrixCoefficientBase& a, const MatrixCoefficientBase& b)
      : m_a(a.copy()), m_b(b.copy())
   {}

   Dot<MatrixCoefficientBase>::Dot(const Dot& other)
      : m_a(other.m_a->copy()), m_b(other.m_b->copy())
   {}

   double Dot<MatrixCoefficientBase>::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      mfem::DenseMatrix ma, mb;
      m_a->getValue(ma, trans, ip);
      m_b->getValue(mb, trans, ip);
      return ma * mb;
   }
}
