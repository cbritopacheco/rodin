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
   Dot<FunctionBase, FunctionBase>::traceOf(const std::set<int>& attrs)
   {
      ScalarFunctionBase::traceOf(attrs);
      m_a->traceOf(attrs);
      m_b->traceOf(attrs);
      return *this;
   }

   double Dot<FunctionBase, FunctionBase>::getValue(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const
   {
      assert(m_a->getRangeShape() == m_b->getRangeShape());
      mfem::DenseMatrix a, b;
      m_a->getValue(a, trans, ip);
      m_b->getValue(b, trans, ip);
      return a * b;
   }

   // ---- Dot<BasisOperator, BasisOperator> ---------------------------------
   Dot<BasisOperator, BasisOperator>::Dot(
         std::unique_ptr<BasisOperator> lhs, std::unique_ptr<BasisOperator> rhs)
      : m_result(lhs->getDOFs(), rhs->getDOFs())
   {
      if (lhs->isDense() && rhs->isDense())
      {
         const auto& denseLhs = static_cast<const DenseBasisOperator&>(*lhs);
         const auto& denseRhs = static_cast<const DenseBasisOperator&>(*rhs);

         for (int i = 0; i < denseLhs.getDOFs(); i++)
            for (int j = 0; j < denseRhs.getDOFs(); j++)
               m_result(i, j) = denseLhs(i) * denseRhs(j);
      }
      else if (lhs->isDense() && rhs->isSparse())
      {
         auto& denseLhs = static_cast<DenseBasisOperator&>(*lhs);
         const auto& sparseRhs = static_cast<const SparseBasisOperator&>(*rhs);

         for (int i = 0; i < denseLhs.getDOFs(); i++)
         {
            for (int j = 0; j < sparseRhs.getDOFs(); j++)
            {
               lhs->Transpose();
               mfem::DenseMatrix* mult = mfem::Mult(sparseRhs(j), denseLhs(i));
               m_result(i, j) = mult->Trace();
               delete mult;
            }
         }
      }
      else if (lhs->isSparse() && rhs->isDense())
      {
         auto& denseRhs = static_cast<DenseBasisOperator&>(*rhs);
         const auto& sparseLhs = static_cast<const SparseBasisOperator&>(*lhs);

         for (int i = 0; i < sparseLhs.getDOFs(); i++)
         {
            for (int j = 0; j < denseRhs.getDOFs(); j++)
            {
               lhs->Transpose();
               mfem::DenseMatrix* mult = mfem::Mult(sparseLhs(i), denseRhs(j));
               m_result(i, j) = mult->Trace();
               delete mult;
            }
         }
      }
      else if (lhs->isSparse() && rhs->isSparse())
      {
         const auto& sparseLhs = static_cast<const SparseBasisOperator&>(*lhs);
         const auto& sparseRhs = static_cast<const SparseBasisOperator&>(*rhs);

         for (int i = 0; i < sparseLhs.getDOFs(); i++)
         {
            for (int j = 0; j < sparseRhs.getDOFs(); j++)
            {
               mfem::SparseMatrix* mult = mfem::TransposeMult(sparseLhs(i), sparseLhs(j));
               mfem::Vector diag;
               mult->GetDiag(diag);
               m_result(i, j) = diag.Sum();
               delete mult;
            }
         }
      }
      else
      {
         assert(false);
      }
   }

   Dot<BasisOperator, BasisOperator>::Dot(const Dot& other)
      : m_result(other.m_result)
   {}

   Dot<BasisOperator, BasisOperator>::Dot(Dot&& other)
      : m_result(std::move(other.m_result))
   {}
}
