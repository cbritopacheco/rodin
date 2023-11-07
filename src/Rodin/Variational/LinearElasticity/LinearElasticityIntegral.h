#ifndef RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  template <class FESType, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator final
    : public BilinearFormIntegratorBase
  {
    public:
      using FES = FESType;
      using Parent = BilinearFormIntegratorBase;
      using Mu = FunctionBase<MuDerived>;
      using Lambda = FunctionBase<LambdaDerived>;

    private:
        using MuRange = typename FormLanguage::Traits<Mu>::RangeType;
        using LambdaRange = typename FormLanguage::Traits<Lambda>::RangeType;
        static_assert(std::is_same_v<MuRange, Scalar>);
        static_assert(std::is_same_v<LambdaRange, Scalar>);

    public:
      LinearElasticityIntegrator(
          const TrialFunction<FES>& u, const TestFunction<FES>& v,
          const Lambda& lambda, const Mu& mu)
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
        const QF::GenericPolytopeQuadrature qf(polytope.getGeometry());
        auto& res = getMatrix();
        res = Math::Matrix::Zero(dofs, dofs);
        m_rjac.resize(dofs);
        m_pjac.resize(dofs);
        m_psym.resize(dofs);
        m_pdiv.resize(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          const auto& rc = qf.getPoint(k);
          const Scalar w = qf.getWeight(k);
          const Geometry::Point p(polytope, trans, std::cref(rc));
          const Scalar distortion = p.getDistortion();
          const Scalar mu = getMu().getValue(p);
          const Scalar lambda = getLambda().getValue(p);
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
        }
        res.template triangularView<Eigen::Upper>() = res.transpose();
      }

      inline
      constexpr
      const Mu& getMu() const
      {
        assert(m_mu);
        return *m_mu;
      }

      inline
      constexpr
      const Lambda& getLambda() const
      {
        assert(m_lambda);
        return *m_lambda;
      }

      inline
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_fes.get();
      }

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Domain;
      }

      inline LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }

    private:
      std::unique_ptr<Lambda> m_lambda;
      std::unique_ptr<Mu> m_mu;
      std::reference_wrapper<const FES> m_fes;

      std::vector<Math::SpatialMatrix> m_rjac;
      std::vector<Math::SpatialMatrix> m_pjac;
      std::vector<Math::SpatialMatrix> m_psym;
      Math::Vector m_pdiv;
  };

  template <class FES, class LambdaDerived, class MuDerived>
  LinearElasticityIntegrator(
      const TrialFunction<FES>&, const TestFunction<FES>&,
      const FunctionBase<LambdaDerived>&, const FunctionBase<MuDerived>&)
    -> LinearElasticityIntegrator<FES, LambdaDerived, MuDerived>;

  template <class FES>
  class LinearElasticityIntegral final
  {
    public:
      LinearElasticityIntegral(const TrialFunction<FES>& u, const TestFunction<FES>& v)
        : m_u(u), m_v(v)
      {}

      template <class L, class M>
      inline
      constexpr
      auto
      operator()(const L& lambda, const M& mu) const
      {
        return LinearElasticityIntegrator(m_u.get(), m_v.get(),
            ScalarFunction<L>(lambda), ScalarFunction<M>(mu));
      }

      template <class LambdaDerived, class MuDerived>
      inline
      constexpr
      auto
      operator()(const FunctionBase<LambdaDerived>& lambda, const FunctionBase<MuDerived>& mu) const
      {
        return LinearElasticityIntegrator(m_u.get(), m_v.get(), lambda, mu);
      }

    private:
      std::reference_wrapper<const TrialFunction<FES>> m_u;
      std::reference_wrapper<const TestFunction<FES>>  m_v;
  };

  template <class FES>
  LinearElasticityIntegral(const TrialFunction<FES>&, const TestFunction<FES>&)
    -> LinearElasticityIntegral<FES>;
}

#endif

