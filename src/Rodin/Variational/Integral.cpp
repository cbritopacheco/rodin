#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   mfem::DenseMatrix
   Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   ::getElementMatrix(const Geometry::Simplex& element) const
   {
      mfem::DenseMatrix mat;

      auto& trial = m_prod.getLHS();
      auto& test = m_prod.getRHS();

      auto& trans = element.getTransformation();

      mat.SetSize(test.getDOFs(element), trial.getDOFs(element));
      mat = 0.0;

      const int order = getIntegrationOrder(
            trial.getFiniteElementSpace(), test.getFiniteElementSpace(), element);

      mfem::DenseMatrix tmp;
      ShapeComputator shapeCompute;
      for (const auto& p : element.getIntegrationRule(order))
      {
         m_prod.getElementMatrix(tmp, shapeCompute, p);
         mfem::Add(mat, tmp, trans.Weight() * trans.GetIntPoint().weight, mat);
      }

      return mat;
   }

   mfem::Vector Integral<ShapeFunctionBase<TestSpace>>::getElementVector(
         const Geometry::Simplex& element) const
   {
      mfem::Vector vec;

      auto& test = *m_integrand;

      auto& trans = element.getTransformation();

      assert(test.getRangeType() == RangeType::Scalar);

      vec.SetSize(test.getDOFs(element));
      vec = 0.0;

      const int order = getIntegrationOrder(test.getFiniteElementSpace(), element);

      ShapeComputator compute;
      for (const auto& p : element.getIntegrationRule(order))
      {
         DenseBasisOperator testOp;
         test.getOperator(testOp, compute, p);
         testOp *= trans.Weight() * trans.GetIntPoint().weight;
         testOp.addToVector(vec);
      }

      return vec;
   }
}
