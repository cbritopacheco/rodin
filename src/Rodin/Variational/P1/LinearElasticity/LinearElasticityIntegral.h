#ifndef RODIN_VARIATIONAL_P1_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_P1_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

#include "Rodin/Variational/LinearElasticity/LinearElasticityIntegral.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/QF1P1.h"

namespace Rodin::Variational
{
  template <class Context, class MuDerived, class LambdaDerived>
  class LinearElasticityIntegrator<VectorP1<Context>, MuDerived, LambdaDerived> final
    : public LocalBilinearFormIntegratorBase
  {
    public:
      using FESType = VectorP1<Context>;

      using MuType = FunctionBase<MuDerived>;

      using LambdaType = FunctionBase<LambdaDerived>;

      using Parent = LocalBilinearFormIntegratorBase;

    private:
        using MuRangeType = typename FormLanguage::Traits<MuType>::RangeType;

        using LambdaRangeType = typename FormLanguage::Traits<LambdaType>::RangeType;

        static_assert(std::is_same_v<MuRangeType, Scalar>);
        static_assert(std::is_same_v<LambdaRangeType, Scalar>);

    public:
      LinearElasticityIntegrator(
          const TrialFunction<FESType>& u, const TestFunction<FESType>& v,
          const LambdaType& lambda, const MuType& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()), m_mu(mu.copy()),
          m_fes(u.getFiniteElementSpace())
      {
        assert(u.getFiniteElementSpace() == v.getFiniteElementSpace());
      }

      LinearElasticityIntegrator(const LinearElasticityIntegrator& other)
        : Parent(other),
          m_lambda(other.m_lambda->copy()), m_mu(other.m_mu->copy()),
          m_fes(other.m_fes)
      {}

      LinearElasticityIntegrator(LinearElasticityIntegrator&& other)
        : Parent(std::move(other)),
          m_lambda(std::move(other.m_lambda)), m_mu(std::move(other.mu)),
          m_fes(std::move(other.m_fes))
      {}

      void assemble(const Geometry::Polytope& polytope) override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& fe = m_fes.get().getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        const QF::QF1P1 qf(polytope.getGeometry());
        assert(qf.getSize() == 1);
        const Scalar w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        const Geometry::Point p(polytope, trans, std::cref(rc));
        const Scalar distortion = p.getDistortion();
        const Scalar mu = getMu().getValue(p);
        const Scalar lambda = getLambda().getValue(p);
        auto& res = getMatrix();
        res = Math::Matrix<Scalar>::Zero(dofs, dofs);
        m_rjac.resize(dofs);
        m_pjac.resize(dofs);
        m_psym.resize(dofs);
        m_pdiv.resize(dofs);
        for (size_t local = 0; local < dofs; local++)
        {
          fe.getJacobian(local)(m_rjac[local], rc);
          m_pjac[local] = m_rjac[local] * p.getJacobianInverse();
          m_psym[local] = m_pjac[local] + m_pjac[local].transpose();
          m_pdiv[local] = m_pjac[local].trace();
        }
        for (size_t i = 0; i < dofs; i++)
        {
          const Scalar lhs = lambda * m_pdiv[i];
          const Scalar rhs = m_pdiv[i];
          res(i, i) += w * distortion * lhs * rhs;
        }
        for (size_t i = 0; i < dofs; i++)
        {
          const auto lhs = mu * m_psym[i];
          const auto rhs = 0.5 * m_psym[i];
          res(i, i) += w * distortion * (lhs.array() * rhs.array()).rowwise().sum().colwise().sum().value();
        }
        for (size_t i = 0; i < dofs; i++)
        {
          const Scalar lhs = lambda * m_pdiv[i];
          for (size_t j = 0; j < i; j++)
          {
            const Scalar rhs = m_pdiv[j];
            res(i, j) += w * distortion * lhs * rhs;
          }
        }
        for (size_t i = 0; i < dofs; i++)
        {
          const auto lhs = mu * m_psym[i];
          for (size_t j = 0; j < i; j++)
          {
            const auto rhs = 0.5 * m_psym[j];
            res(i, j) += w * distortion * (lhs.array() * rhs.array()).rowwise().sum().colwise().sum().value();
          }
        }
        res.template triangularView<Eigen::Upper>() = res.transpose();
      }

      inline
      constexpr
      const MuType& getMu() const
      {
        assert(m_mu);
        return *m_mu;
      }

      inline
      constexpr
      const LambdaType& getLambda() const
      {
        assert(m_lambda);
        return *m_lambda;
      }

      inline
      constexpr
      const FESType& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Cells;
      }

      inline LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }

    private:
      std::unique_ptr<LambdaType> m_lambda;
      std::unique_ptr<MuType> m_mu;
      std::reference_wrapper<const FESType> m_fes;

      std::vector<Math::SpatialMatrix<Scalar>> m_rjac;
      std::vector<Math::SpatialMatrix<Scalar>> m_pjac;
      std::vector<Math::SpatialMatrix<Scalar>> m_psym;
      Math::Vector<Scalar> m_pdiv;
  };
}
#endif

