#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   void
   Integral<FormLanguage::Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
   ::getElementMatrix(
         const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
   {
      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      assert(trial.getColumns(trialElement, trans) == test.getColumns(testElement, trans));

      mat.SetSize(test.getRows(trialElement, trans), trial.getRows(testElement, trans));
      mat = 0.0;

      mfem::DenseMatrix trialOp;
      mfem::DenseMatrix testOp;
      trialOp.SetSize(trial.getRows(trialElement, trans), trial.getColumns(trialElement, trans));
      testOp.SetSize(test.getRows(testElement, trans), test.getColumns(testElement, trans));

      int order = trialElement.GetOrder() + testElement.GetOrder() + trans.OrderW();
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         trial.getOperator(trialElement, trans, trialOp);
         test.getOperator(testElement, trans, testOp);
         mfem::AddMult_a_ABt(trans.Weight() * ip.weight, testOp, trialOp, mat);
      }
   }

   void Integral<ShapeFunctionBase<Test>>::getElementVector(
         const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans, mfem::Vector& vec)
   {
      auto& test = *m_test;

      assert(test.getColumns(fe, trans) == 1);

      vec.SetSize(test.getRows(fe, trans));
      vec = 0.0;

      mfem::DenseMatrix testOp;
      testOp.SetSize(test.getRows(fe, trans), 1);

      int order = fe.GetOrder() + trans.OrderW();
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         test.getOperator(fe, trans, testOp);
         testOp *= trans.Weight() * ip.weight;
         testOp.AddToVector(0, vec);
      }
   }
}
