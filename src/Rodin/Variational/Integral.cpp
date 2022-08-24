#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   void
   Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   ::getElementMatrix(const Bilinear::Assembly::Common& as) const
   {
      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      as.mat.SetSize(test.getDOFs(as.trial, as.trans), trial.getDOFs(as.test, as.trans));
      as.mat = 0.0;

      const int order = getIntegrationOrder(as);
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(as.trans.GetGeometryType(), order);
      ShapeComputator shapeCompute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         as.trans.SetIntPoint(&ip);
         mfem::Add(as.mat,
               m_prod.getElementMatrix(as.trial, as.test, as.trans, ip, shapeCompute),
               as.trans.Weight() * ip.weight,
               as.mat);
      }
   }

   void Integral<ShapeFunctionBase<TestSpace>>::getElementVector(
         const Linear::Assembly::Common& as) const
   {
      const mfem::FiniteElement& fe = as.fe;
      mfem::ElementTransformation& trans = as.trans;
      mfem::Vector& vec = as.vec;

      auto& test = *m_integrand;

      assert(test.getRangeType() == RangeType::Scalar);

      vec.SetSize(test.getDOFs(fe, trans));
      vec = 0.0;

      const int order = getIntegrationOrder(as);
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);

      ShapeComputator compute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         DenseBasisOperator testOp;
         test.getOperator(testOp, fe, trans, ip, compute);
         testOp *= trans.Weight() * ip.weight;
         testOp.addToVector(vec);
      }
   }
}
