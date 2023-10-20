#ifndef RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

namespace Rodin::Variational
{
  template <class TrialFES, class TestFES, class MuDerived, class LambdaDerived>
  class LinearElasticityIntegrator final
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;
      using Mu = FunctionBase<MuDerived>;
      using Lambda = FunctionBase<LambdaDerived>;

      LinearElasticityIntegrator(
          const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v,
          const Lambda& lambda, const Mu& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()), m_mu(mu.copy()),
          m_fes(u.getFiniteElementSpace())
      {
        using MuRange = typename FormLanguage::Traits<Mu>::RangeType;
        using LambdaRange = typename FormLanguage::Traits<Lambda>::RangeType;
        static_assert(std::is_same_v<TrialFES, TestFES>);
        static_assert(std::is_same_v<MuRange, Scalar>);
        static_assert(std::is_same_v<LambdaRange, Scalar>);
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

      void assemble(const Geometry::Polytope& polytope) override;

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
      const TrialFES& getFiniteElementSpace() const
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
      std::reference_wrapper<const TrialFES> m_fes;
  };

  template <class TrialFES, class TestFES, class LambdaDerived, class MuDerived>
  LinearElasticityIntegrator(
      const TrialFunction<TrialFES>&, const TestFunction<TestFES>&,
      const FunctionBase<LambdaDerived>&, const FunctionBase<MuDerived>&)
    -> LinearElasticityIntegrator<TrialFES, TestFES, LambdaDerived, MuDerived>;

  template <class TrialFES, class TestFES>
  class LinearElasticityIntegral final
  {
    public:
      LinearElasticityIntegral(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
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

      template <class MuDerived, class LambdaDerived>
      inline
      constexpr
      auto
      operator()(const FunctionBase<LambdaDerived>& lambda, const FunctionBase<MuDerived>& mu) const
      {
        return LinearElasticityIntegrator(m_u.get(), m_v.get(), lambda, mu);
      }

    private:
      std::reference_wrapper<const TrialFunction<TrialFES>> m_u;
      std::reference_wrapper<const TestFunction<TestFES>>   m_v;
  };

  template <class TrialFES, class TestFES>
  LinearElasticityIntegral(const TrialFunction<TrialFES>&, const TestFunction<TestFES>&)
    -> LinearElasticityIntegral<TrialFES, TestFES>;
}

#endif

