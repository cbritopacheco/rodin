#include "Rodin/Alert.h"
#include "Rodin/Geometry/SubMesh.h"

#include "Function.h"
#include "Transpose.h"

namespace Rodin::Variational
{
   Transpose<FunctionBase> FunctionBase::T() const
   {
      return Transpose<FunctionBase>(*this);
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

   Internal::MFEMFunction FunctionBase::build(const Geometry::MeshBase& mesh) const
   {
      switch (getRangeType())
      {
         case RangeType::Scalar:
            return Internal::MFEMFunction(new Internal::ScalarProxyFunction(mesh, *this));
         case RangeType::Vector:
            return Internal::MFEMFunction(new Internal::VectorProxyFunction(mesh, *this));
         case RangeType::Matrix:
            return Internal::MFEMFunction(new Internal::MatrixProxyFunction(mesh, *this));
      }
      // The following return is needed to prevent compiler warnings/errors
      return Internal::MFEMFunction(nullptr);
   }
}

namespace Rodin::Variational::Internal
{
   ProxyFunction<RangeType::Scalar>::ProxyFunction(
         const Geometry::MeshBase& mesh, const FunctionBase& s)
      :  m_mesh(mesh),
         m_s(s)
   {
      assert(s.getRangeType() == RangeType::Scalar);
   }

   ProxyFunction<RangeType::Scalar>::ProxyFunction(const ProxyFunction& other)
      : mfem::Coefficient(other),
        m_mesh(other.m_mesh),
        m_s(other.m_s)
   {}

   ProxyFunction<RangeType::Scalar>::ProxyFunction(ProxyFunction&& other)
      : mfem::Coefficient(std::move(other)),
        m_mesh(other.m_mesh),
        m_s(other.m_s)
   {}

   double ProxyFunction<RangeType::Scalar>::Eval(
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip)
   {
      double res = 0;
      switch (trans.ElementType)
      {
         case mfem::ElementTransformation::ELEMENT:
         {
            res = m_s.getValue(Geometry::Point(*m_mesh.getElement(trans.ElementNo), ip)).scalar();
            break;
         }
         case mfem::ElementTransformation::BDR_ELEMENT:
         {
            int faceIdx = m_mesh.getHandle().GetBdrFace(trans.ElementNo);
            res = m_s.getValue(Geometry::Point(*m_mesh.getFace(faceIdx), ip)).scalar();
            break;
         }
         case mfem::ElementTransformation::FACE:
         {
            res = m_s.getValue(Geometry::Point(*m_mesh.getFace(trans.ElementNo), ip)).scalar();
            break;
         }
      }
      return res;
   }

   ProxyFunction<RangeType::Vector>::ProxyFunction(
         const Geometry::MeshBase& mesh, const FunctionBase& s)
      :  mfem::VectorCoefficient(
            s.getRangeShape().height() == 1 ?
               s.getRangeShape().width() : s.getRangeShape().height()),
         m_mesh(mesh),
         m_s(s)
   {
      assert(s.getRangeType() == RangeType::Vector);
   }

   ProxyFunction<RangeType::Vector>::ProxyFunction(const ProxyFunction& other)
      : mfem::VectorCoefficient(other),
        m_mesh(other.m_mesh),
        m_s(other.m_s)
   {}

   ProxyFunction<RangeType::Vector>::ProxyFunction(ProxyFunction&& other)
      : mfem::VectorCoefficient(std::move(other)),
        m_mesh(other.m_mesh),
        m_s(other.m_s)
   {}

   void ProxyFunction<RangeType::Vector>::Eval(
         mfem::Vector& value,
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip)
   {
      switch (trans.ElementType)
      {
         case mfem::ElementTransformation::ELEMENT:
         {
            value = m_s.getValue(Geometry::Point(*m_mesh.getElement(trans.ElementNo), ip)).vector();
            break;
         }
         case mfem::ElementTransformation::BDR_ELEMENT:
         {
            int faceIdx = m_mesh.getHandle().GetBdrFace(trans.ElementNo);
            value = m_s.getValue(Geometry::Point(*m_mesh.getFace(faceIdx), ip)).vector();
            break;
         }
         case mfem::ElementTransformation::FACE:
         {
            value = m_s.getValue(Geometry::Point(*m_mesh.getFace(trans.ElementNo), ip)).vector();
            break;
         }
      }
   }

   ProxyFunction<RangeType::Matrix>::ProxyFunction(
         const Geometry::MeshBase& mesh, const FunctionBase& s)
      :  mfem::MatrixCoefficient(s.getRangeShape().height(), s.getRangeShape().width()),
         m_mesh(mesh),
         m_s(s)
   {
      assert(s.getRangeType() == RangeType::Matrix);
   }

   ProxyFunction<RangeType::Matrix>::ProxyFunction(const ProxyFunction& other)
      : mfem::MatrixCoefficient(other),
        m_mesh(other.m_mesh),
        m_s(other.m_s)
   {}

   ProxyFunction<RangeType::Matrix>::ProxyFunction(ProxyFunction&& other)
      : mfem::MatrixCoefficient(std::move(other)),
        m_mesh(other.m_mesh),
        m_s(other.m_s)
   {}

   void ProxyFunction<RangeType::Matrix>::Eval(
         mfem::DenseMatrix& value, mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip)
   {
      switch (trans.ElementType)
      {
         case mfem::ElementTransformation::ELEMENT:
         {
            value = m_s.getValue(Geometry::Point(*m_mesh.getElement(trans.ElementNo), ip)).matrix();
            break;
         }
         case mfem::ElementTransformation::BDR_ELEMENT:
         {
            int faceIdx = m_mesh.getHandle().GetBdrFace(trans.ElementNo);
            value = m_s.getValue(Geometry::Point(*m_mesh.getFace(faceIdx), ip)).matrix();
            break;
         }
         case mfem::ElementTransformation::FACE:
         {
            value = m_s.getValue(Geometry::Point(*m_mesh.getFace(trans.ElementNo), ip)).matrix();
            break;
         }
      }
   }
}

