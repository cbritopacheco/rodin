#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"

namespace Rodin::Variational
{
   Integral<FormLanguage::Product<TrialFunctionBase, TestFunctionBase>>
   ::Integral(const FormLanguage::Product<TrialFunctionBase, TestFunctionBase>& prod)
      : m_prod(prod)
   {}

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

            int order = trialElement.GetOrder() + testElement.GetOrder() + trans.OrderW();
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
}
