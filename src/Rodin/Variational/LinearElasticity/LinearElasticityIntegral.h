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
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<FES>::NumberType>
  {
    public:
      using FESType = FES;

      using NumberType = typename FormLanguage::Traits<FES>::NumberType;

      using Mu = FunctionBase<MuDerived>;

      using Lambda = FunctionBase<LambdaDerived>;

      using Parent = LocalBilinearFormIntegratorBase<NumberType>;

    private:
        using MuRangeType = typename FormLanguage::Traits<Mu>::RangeType;

        using LambdaRangeType = typename FormLanguage::Traits<Lambda>::RangeType;

        static_assert(std::is_same_v<MuRangeType, NumberType>);

        static_assert(std::is_same_v<LambdaRangeType, NumberType>);

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
        return *this;
      }

      virtual NumberType integrate(size_t tr, size_t te) override
      {
        assert(false);
        return 0;
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

