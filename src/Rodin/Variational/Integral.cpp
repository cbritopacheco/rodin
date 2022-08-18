#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   void
   Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   ::getElementMatrix(
         const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) const
   {
      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      mat.SetSize(test.getDOFs(trialElement, trans), trial.getDOFs(testElement, trans));
      mat = 0.0;

      int order = getIntegrationOrder(trialElement, testElement, trans);
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);
      ShapeComputator shapeCompute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         mfem::Add(mat,
               m_prod.getElementMatrix(trialElement, testElement, trans, ip, shapeCompute),
               trans.Weight() * ip.weight,
               mat);
      }
   }

   void Integral<ShapeFunctionBase<TestSpace>>::getElementVector(const Assembly::Common& as) const
   {
      const mfem::FiniteElement& fe = as.fe;
      mfem::ElementTransformation& trans = as.trans;
      mfem::Vector& vec = as.vec;

      auto& test = *m_integrand;

      assert(test.getRangeType() == RangeType::Scalar);

      vec.SetSize(test.getDOFs(fe, trans));
      vec = 0.0;

      int order = getIntegrationOrder(as);
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
