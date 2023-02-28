#ifndef RODIN_VARIATIONAL_GAUSSIANQUADRATURE_H
#define RODIN_VARIATIONAL_GAUSSIANQUADRATURE_H

#include "Dot.h"
#include "ShapeFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{
  template <class Integrand>
  class GaussianQuadrature;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class GaussianQuadrature<
    Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;
      using Integrand = Dot<LHS, RHS>;
      using Parent = BilinearFormIntegratorBase;

      constexpr
      GaussianQuadrature(const LHS& lhs, const RHS& rhs)
        : GaussianQuadrature(Dot(lhs, rhs))
      {}

      constexpr
      GaussianQuadrature(const Integrand& prod)
        : BilinearFormIntegratorBase(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
          m_prod(prod.copy())
      {}

      constexpr
      GaussianQuadrature(const GaussianQuadrature& other)
        : BilinearFormIntegratorBase(other),
          m_prod(other.m_prod->copy())
      {}

      constexpr
      GaussianQuadrature(GaussianQuadrature&& other)
        : BilinearFormIntegratorBase(std::move(other)),
          m_prod(std::move(other.m_prod))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_prod);
        return *m_prod;
      }

      Math::Matrix getMatrix(const Geometry::Simplex& simplex) const final override
      {
        const auto& integrand = getIntegrand();
        const auto& trial = integrand.getLHS();
        const auto& test = integrand.getRHS();
        const auto& trans = simplex.getTransformation();
        const size_t order =
          trial.getFiniteElementSpace().getOrder(simplex) +
          test.getFiniteElementSpace().getOrder(simplex) +
          simplex.getTransformation().getHandle().OrderW();
        const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(static_cast<mfem::Geometry::Type>(simplex.getGeometry()), order);
        Math::Matrix res = Math::Matrix::Zero(test.getDOFs(simplex), trial.getDOFs(simplex));
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
          const auto& ip = ir->IntPoint(i);
          Geometry::Point p(simplex, trans, Internal::ip2vec(ir->IntPoint(i), simplex.getDimension()));
          res += ip.weight * p.getDistortion() * integrand.getMatrix(p);
        }
        return res;
      }

      virtual Region getRegion() const override = 0;

      virtual GaussianQuadrature* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_prod;
  };

  template <class NestedDerived, class FES>
  class GaussianQuadrature<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = LinearFormIntegratorBase;

      template <class LHSDerived, class RHSDerived>
      constexpr
      GaussianQuadrature(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : GaussianQuadrature(Dot(lhs, rhs))
      {}

      constexpr
      GaussianQuadrature(const Integrand& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      GaussianQuadrature(const GaussianQuadrature& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      GaussianQuadrature(GaussianQuadrature&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      Math::Vector getVector(const Geometry::Simplex& simplex) const final override
      {
        const auto& integrand = getIntegrand();
        assert(integrand.getRangeType() == RangeType::Scalar);
        const auto& trans = simplex.getTransformation();
        const size_t order =
          integrand.getFiniteElementSpace().getOrder(simplex) + trans.getHandle().OrderW();
        const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(static_cast<mfem::Geometry::Type>(simplex.getGeometry()), order);
        Math::Vector res = Math::Vector::Zero(integrand.getDOFs(simplex));
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
          const auto& ip = ir->IntPoint(i);
          Geometry::Point p(simplex, trans, Internal::ip2vec(ip, simplex.getDimension()));
          res += ip.weight * p.getDistortion() * integrand.getTensorBasis(p).getVector();
        }
        return res;
      }

      virtual Region getRegion() const override = 0;

      virtual GaussianQuadrature* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };
}

#endif
