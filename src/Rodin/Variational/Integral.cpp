#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarCoefficient.h"

namespace Rodin::Variational
{
   // ---- FormLanguage::Product<TrialFunctionBase, TestFunctionBase> --------
   void
   Integral<FormLanguage::Product<TrialFunctionBase, TestFunctionBase>>
   ::getElementMatrix(
         const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
   {
      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      int trialDof = trialElement.GetDof(),
          testDof = testElement.GetDof();
      int spaceDim = trans.GetSpaceDim();

      mat.SetSize(testDof, trialDof);
      mat = 0.0;

      assert(trial.getValueType() == test.getValueType());

      switch (test.getValueType())
      {
         case ShapeFunction::Scalar:
         {
            break;
         }
         case ShapeFunction::Vector:
         {
            mfem::DenseMatrix trialValue(trialDof, spaceDim);
            mfem::DenseMatrix testValue(testDof, spaceDim);
            trialValue.SetSize(trialDof, spaceDim);
            testValue.SetSize(testDof, spaceDim);

            int order = trialElement.GetOrder() + testElement.GetOrder() + trans.OrderW() - 1;
            const mfem::IntegrationRule* ir =
               &mfem::IntRules.Get(trans.GetGeometryType(), order);
            for (int i = 0; i < ir->GetNPoints(); i++)
            {
               const mfem::IntegrationPoint &ip = ir->IntPoint(i);
               trans.SetIntPoint(&ip);
               trial.getValue(trialElement, trans, trialValue);
               test.getValue(testElement, trans, testValue);
               double w = trans.Weight() * ip.weight;
               mfem::AddMult_a_ABt(w, testValue, trialValue, mat);
            }
            break;
         }
         case ShapeFunction::Matrix:
         {
            break;
         }
      }
   }

   // ---- FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase> ----
   void
   Integral<FormLanguage::Product<ScalarCoefficientBase, TestFunctionBase>>
   ::getElementVector(
         const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans, mfem::Vector& vec)
   {
      auto& scalar = m_prod.getLHS();
      auto& test = m_prod.getRHS();
      int dof = fe.GetDof();

      switch (test.getValueType())
      {
         case ShapeFunction::Scalar:
         {
            vec.SetSize(dof);
            vec = 0.0;

            mfem::Vector value;
            value.SetSize(dof);

            int order = fe.GetOrder() + trans.OrderW() - 1;
            const mfem::IntegrationRule* ir =
               &mfem::IntRules.Get(trans.GetGeometryType(), order);

            for (int i = 0; i < ir->GetNPoints(); i++)
            {
               const mfem::IntegrationPoint& ip = ir->IntPoint(i);
               trans.SetIntPoint (&ip);
               double val = ip.weight * trans.Weight() * scalar.getValue(trans, ip);

               test.getValue(fe, trans, vec);
               add(vec, ip.weight * val, value, vec);
            }
            break;
         }
         case ShapeFunction::Vector:
         {
            Alert::Exception()
               << "Mismatching dimensions. Here goeth my AST log."
               << Alert::Raise;
            break;
         }
         case ShapeFunction::Matrix:
         {
            Alert::Exception()
               << "Mismatching dimensions. Here goeth my AST log."
               << Alert::Raise;
            break;
         }
      }
   }
}
