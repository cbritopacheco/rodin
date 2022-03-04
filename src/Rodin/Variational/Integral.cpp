#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   void
   Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
   ::getElementMatrix(
         const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
   {
      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      mat.SetSize(test.getDOFs(trialElement, trans), trial.getDOFs(testElement, trans));
      mat = 0.0;

      int order = trialElement.GetOrder() + testElement.GetOrder() + trans.OrderW();
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         mfem::Add(mat,
               m_prod.getElementMatrix(trialElement, testElement, trans),
               trans.Weight() * ip.weight,
               mat);
      }
   }

   void Integral<ShapeFunctionBase<Test>>::getElementVector(
         const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans, mfem::Vector& vec)
   {
      auto& test = *m_test;

      assert(test.getRows(fe, trans) == 1);
      assert(test.getColumns(fe, trans) == 1);

      vec.SetSize(test.getDOFs(fe, trans));
      vec = 0.0;

      int order = fe.GetOrder() + trans.OrderW();
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         auto testOp = test.getOperator(fe, trans);
         (*testOp) *= trans.Weight() * ip.weight;
         testOp->AddToVector(vec);
      }
   }
}
