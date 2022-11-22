#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   void
   Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   ::getElementMatrix(const Bilinear::Assembly::Native& as) const
   {
      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      as.matrix.SetSize(test.getDOFs(as.element), trial.getDOFs(as.element));
      as.matrix = 0.0;

      mfem::ElementTransformation& trans = as.element.getTransformation();

      const int order = getIntegrationOrder(*this, as);
      const mfem::IntegrationRule* ir = &mfem::IntRules.Get(trans.GetGeometryType(), order);

      mfem::DenseMatrix tmp;
      ShapeComputator shapeCompute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         m_prod.getElementMatrix(tmp, shapeCompute, as.element);
         mfem::Add(as.matrix, tmp, trans.Weight() * ip.weight, as.matrix);
      }
   }

   void
   InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   ::getElementMatrix(const Bilinear::Assembly::Native& as) const
   {
      assert(false);
      // auto& trial = m_prod.getLHS();
      // auto& test = m_prod.getRHS();

      // const auto& data = std::get<Bilinear::Assembly::Native::Interface>(as.data);

      // as.mat.SetSize(
      //       test.getDOFs(data.test1, data.test2, data.trans)
      //       test.getDOFs(data.test1, data.trans) + test.getDOFs(data.test2, data.trans),
      //       trial.getDOFs(data.trial1, data.trans) + trial.getDOFs(data.trial2, data.trans));
      // as.mat = 0.0;

      // const int order = getIntegrationOrder(as);
      // const mfem::IntegrationRule* ir =
      //    &mfem::IntRules.Get(data.trans.GetGeometryType(), order);

      // mfem::DenseMatrix mat;
      // ShapeComputator shapeCompute;
      // for (int i = 0; i < ir->GetNPoints(); i++)
      // {
      //    const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      //    data.trans.SetIntPoint(&ip);
      //    m_prod.getElementMatrix(mat, shapeCompute, data);
      //    mfem::Add(as.mat, mat, data.trans.Weight() * ip.weight, as.mat);
      // }
   }

   void Integral<ShapeFunctionBase<TestSpace>>::getElementVector(
         const Linear::Assembly::Native& as) const
   {
      mfem::Vector& vec = as.vec;

      auto& test = *m_integrand;

      assert(test.getRangeType() == RangeType::Scalar);

      vec.SetSize(test.getDOFs(as.element));
      vec = 0.0;

      mfem::ElementTransformation& trans = as.element.getTransformation();

      const int order = getIntegrationOrder(*this, as);
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);

      ShapeComputator compute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         DenseBasisOperator testOp;
         test.getOperator(testOp, compute, as.element);
         testOp *= trans.Weight() * ip.weight;
         testOp.addToVector(vec);
      }
   }
}
