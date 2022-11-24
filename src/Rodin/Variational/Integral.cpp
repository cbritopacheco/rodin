#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   mfem::DenseMatrix
   Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   ::getElementMatrix(const Geometry::SimplexBase& element) const
   {
      mfem::DenseMatrix mat;

      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      mat.SetSize(test.getDOFs(element), trial.getDOFs(element));
      mat = 0.0;

      mfem::ElementTransformation& trans = element.getTransformation();

      const int order = getIntegrationOrder(trial.getFiniteElementSpace(), test.getFiniteElementSpace(), element);
      const mfem::IntegrationRule* ir = &mfem::IntRules.Get(trans.GetGeometryType(), order);

      mfem::DenseMatrix tmp;
      ShapeComputator shapeCompute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         m_prod.getElementMatrix(tmp, shapeCompute, element);
         mfem::Add(mat, tmp, trans.Weight() * ip.weight, mat);
      }

      return mat;
   }

   mfem::Vector Integral<ShapeFunctionBase<TestSpace>>::getElementVector(
         const Geometry::SimplexBase& element) const
   {
      mfem::Vector vec;

      auto& test = *m_integrand;

      assert(test.getRangeType() == RangeType::Scalar);

      vec.SetSize(test.getDOFs(element));
      vec = 0.0;

      mfem::ElementTransformation& trans = element.getTransformation();

      const int order = getIntegrationOrder(test.getFiniteElementSpace(), element);
      const mfem::IntegrationRule* ir =
         &mfem::IntRules.Get(trans.GetGeometryType(), order);

      ShapeComputator compute;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(i);
         trans.SetIntPoint(&ip);
         DenseBasisOperator testOp;
         test.getOperator(testOp, compute, element);
         testOp *= trans.Weight() * ip.weight;
         testOp.addToVector(vec);
      }

      return vec;
   }
}
