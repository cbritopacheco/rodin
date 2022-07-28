#include "Rodin/Alert.h"

#include "Function.h"
#include "Transpose.h"

namespace Rodin::Variational
{
   Transpose<FunctionBase> FunctionBase::T() const
   {
      return Transpose<FunctionBase>(*this);
   }

   mfem::ElementTransformation& FunctionBase::getTraceElementTrans(
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const
   {
      const auto& traceDomain = getTraceDomain();
      switch (trans.ElementType)
      {
         case mfem::ElementTransformation::BDR_ELEMENT:
         {
            int fn = trans.mesh->GetBdrFace(trans.ElementNo);
            if (trans.mesh->FaceIsInterior(fn))
            {
               if (traceDomain.empty())
               {
                  Alert::Warning()
                     << __func__
                     << " called with an empty trace domain. May result in segfault."
                     << Alert::Raise;
                  return trans;
               }

               mfem::FaceElementTransformations* ft =
                  trans.mesh->GetFaceElementTransformations(fn);
               ft->SetAllIntPoints(&ip);
               auto& trans1 = ft->GetElement1Transformation();
               auto& trans2 = ft->GetElement2Transformation();
               if (traceDomain.count(trans1.Attribute))
                  return trans1;
               else if (traceDomain.count(trans2.Attribute))
                  return trans2;
               else
               {
                  /* The boundary over which we are evaluating must be
                   * the interface between the trace domain and some
                   * other domain, i.e. it is not the boundary that was
                   * specified!
                   */
                  assert(false);
                  return trans;
               }
            }
            else
            {
               return trans;
            }
            break;
         }
         default:
         {
            return trans;
         }
      }
   }

   RangeType FunctionBase::getRangeType() const
   {
      auto shape = getRangeShape();
      if (shape.height() == 1 && shape.width() == 1)
      {
         return RangeType::Scalar;
      }
      else if (shape.height() > 1 && shape.width() == 1)
      {
         return RangeType::Vector;
      }
      else
      {
         return RangeType::Matrix;
      }
   }

   std::variant<
      Internal::ScalarProxyFunction,
      Internal::VectorProxyFunction,
      Internal::MatrixProxyFunction> FunctionBase::build() const
   {
      switch (getRangeType())
      {
         case RangeType::Scalar:
            return Internal::ScalarProxyFunction(*this);
         case RangeType::Vector:
            return Internal::VectorProxyFunction(*this);
         case RangeType::Matrix:
            return Internal::MatrixProxyFunction(*this);
      }

      // The following return is needed to prevent compiler warnings/errors
      return Internal::ScalarProxyFunction(*this);
   }
}

namespace Rodin::Variational::Internal
{
   ProxyFunction<RangeType::Scalar>::ProxyFunction(const FunctionBase& s)
      : m_s(s)
   {
      assert(s.getRangeType() == RangeType::Scalar);
   }

   ProxyFunction<RangeType::Scalar>::ProxyFunction(const ProxyFunction& other)
      : mfem::Coefficient(other),
        m_s(other.m_s)
   {}

   ProxyFunction<RangeType::Scalar>::ProxyFunction(ProxyFunction&& other)
      : mfem::Coefficient(std::move(other)),
        m_s(other.m_s)
   {}

   double ProxyFunction<RangeType::Scalar>::Eval(
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip)
   {
      mfem::DenseMatrix v;
      m_s.getValue(v, trans, ip);
      return v(0, 0);
   }

   ProxyFunction<RangeType::Vector>::ProxyFunction(const FunctionBase& s)
      :  mfem::VectorCoefficient(
            s.getRangeShape().height() == 1 ?
               s.getRangeShape().width() : s.getRangeShape().height()),
         m_s(s)
   {
      assert(s.getRangeType() == RangeType::Vector);
   }

   ProxyFunction<RangeType::Vector>::ProxyFunction(const ProxyFunction& other)
      : mfem::VectorCoefficient(other),
        m_s(other.m_s)
   {}

   ProxyFunction<RangeType::Vector>::ProxyFunction(ProxyFunction&& other)
      : mfem::VectorCoefficient(std::move(other)),
        m_s(other.m_s)
   {}

   void ProxyFunction<RangeType::Vector>::Eval(
         mfem::Vector& value,
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip)
   {
      mfem::DenseMatrix v;
      m_s.getValue(v, trans, ip);
      value.SetDataAndSize(v.GetMemory(), v.NumCols() * v.NumRows());
      value.MakeDataOwner();
      v.GetMemory().Reset();
   }

   ProxyFunction<RangeType::Matrix>::ProxyFunction(const FunctionBase& s)
      :  mfem::MatrixCoefficient(s.getRangeShape().height(), s.getRangeShape().width()),
         m_s(s)
   {
      assert(s.getRangeType() == RangeType::Matrix);
   }

   ProxyFunction<RangeType::Matrix>::ProxyFunction(const ProxyFunction& other)
      : mfem::MatrixCoefficient(other),
        m_s(other.m_s)
   {}

   ProxyFunction<RangeType::Matrix>::ProxyFunction(ProxyFunction&& other)
      : mfem::MatrixCoefficient(std::move(other)),
        m_s(other.m_s)
   {}

   void ProxyFunction<RangeType::Matrix>::Eval(
         mfem::DenseMatrix& value, mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip)
   {
      m_s.getValue(value, trans, ip);
   }
}

