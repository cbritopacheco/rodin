#include "Integral.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  Math::Matrix
  Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  ::getMatrix(const Geometry::Simplex& element) const
  {
    const auto& trial = m_prod.getLHS();
    const auto& test = m_prod.getRHS();
    auto& trans = element.getTransformation();
    const size_t order = getIntegrationOrder(
        trial.getFiniteElementSpace(), test.getFiniteElementSpace(), element);
    ShapeComputator shapeCompute;

    Math::Matrix res = Math::Matrix::Zero(test.getDOFs(element), trial.getDOFs(element));
    for (const auto& p : element.getIntegrationRule(order))
    {
      auto tmp = m_prod.getMatrix(shapeCompute, p);
      res += trans.Weight() * trans.GetIntPoint().weight * tmp;
    }
    return res;
  }

  Math::Vector
  Integral<ShapeFunctionBase<TestSpace>>
  ::getVector(const Geometry::Simplex& simplex) const
  {
    const auto& test = *m_integrand;
    assert(test.getRangeType() == RangeType::Scalar);
    auto& trans = simplex.getTransformation();
    const size_t order = getIntegrationOrder(test.getFiniteElementSpace(), simplex);
    ShapeComputator compute;

    Math::Vector res = Math::Vector::Zero(test.getDOFs(simplex));
    for (const auto& p : simplex.getIntegrationRule(order))
    {
      const TensorBasis basis =
        trans.Weight() * trans.GetIntPoint().weight * test.getOperator(compute, p);
      res += Eigen::Map<const Math::Vector>(basis.data(), basis.size());
    }
    return res;
  }
}
