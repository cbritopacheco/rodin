#ifndef RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  template <class FES, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator final
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;

      using Mu = FunctionBase<MuDerived>;

      using Lambda = FunctionBase<LambdaDerived>;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

    private:
        using MuRangeType = typename FormLanguage::Traits<Mu>::RangeType;

        using LambdaRangeType = typename FormLanguage::Traits<Lambda>::RangeType;

        static_assert(std::is_same_v<MuRangeType, ScalarType>);

        static_assert(std::is_same_v<LambdaRangeType, ScalarType>);

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

      virtual const Geometry::Polytope& getPolytope() const override
      {
        return m_polytope.value().get();
      }

      virtual LinearElasticityIntegrator& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_polytope = polytope;
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& fes = getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t order = fe.getOrder();
        m_qf.reset(new QF::GenericPolytopeQuadrature(order, polytope.getGeometry()));
        m_ps.clear();
        m_ps.reserve(m_qf->getSize());
        for (size_t i = 0; i < m_qf->getSize(); i++)
          m_ps.emplace_back(polytope, trans, std::cref(m_qf->getPoint(i)));
        return *this;
      }

      virtual ScalarType integrate(size_t tr, size_t te) override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& fes = getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const auto& qf = *m_qf;
        ScalarType res = 0;
        for (size_t i = 0; i < m_qf->getSize(); i++)
        {
          const auto& p = m_ps[i];
          const auto& w = qf.getWeight(i);
          const auto& rc = qf.getPoint(i);
          const ScalarType mu = getMu().getValue(p);
          const ScalarType lambda = getLambda().getValue(p);
          if (tr == te)
          {
            fe.getJacobian(tr)(m_jac1, rc);
            const auto jac = m_jac1 * p.getJacobianInverse();
            const auto sym = jac + jac.transpose();
            const ScalarType div = jac.trace();
            res += w * p.getDistortion() * (
                getLambda().getValue(p) * div * div + 0.5 * mu * sym.squaredNorm());
          }
          else
          {
            fe.getJacobian(tr)(m_jac1, rc);
            const auto jac1 = m_jac1 * p.getJacobianInverse();
            const auto sym1 = jac1 + jac1.transpose();
            const ScalarType div1 = jac1.trace();

            fe.getJacobian(te)(m_jac2, rc);
            const auto jac2 = m_jac2 * p.getJacobianInverse();
            const auto sym2 = jac2 + jac2.transpose();
            const ScalarType div2 = jac2.trace();
            res += w * p.getDistortion() * (
                getLambda().getValue(p) * div1 * div2 + (
                  (mu * sym2).array() * (0.5 * sym1).array()).rowwise().sum().colwise().sum().value());
          }
        }
        return res;
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
        return Integrator::Region::Cells;
      }

      virtual LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }

    private:
      std::unique_ptr<Lambda> m_lambda;
      std::unique_ptr<Mu> m_mu;
      std::reference_wrapper<const FES> m_fes;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;
      std::vector<Geometry::Point> m_ps;

      Math::SpatialMatrix<Real> m_jac1, m_jac2;
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
            RealFunction<L>(lambda), RealFunction<M>(mu));
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

