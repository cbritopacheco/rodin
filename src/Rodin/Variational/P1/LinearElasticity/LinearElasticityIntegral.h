#ifndef RODIN_VARIATIONAL_P1_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_P1_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

#include "Rodin/Variational/LinearElasticity/LinearElasticityIntegral.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"
#include "Rodin/Variational/P1/P1Element.h"

namespace Rodin::Variational
{
  template <class Context, class MuDerived, class LambdaDerived>
  class LinearElasticityIntegrator<ScalarP1<Context>, ScalarP1<Context>, MuDerived, LambdaDerived> final
    : public BilinearFormIntegratorBase
  {
    public:
      using FES = ScalarP1<Context>;
      using Parent = BilinearFormIntegratorBase;
      using Mu = FunctionBase<MuDerived>;
      using Lambda = FunctionBase<LambdaDerived>;

      LinearElasticityIntegrator(
          const TrialFunction<FES>& u, const TestFunction<FES>& v,
          const Lambda& lambda, const Mu& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()), m_mu(mu.copy()),
          m_fes(u.getFiniteElementSpace())
      {
        using MuRange = typename FormLanguage::Traits<Mu>::RangeType;
        using LambdaRange = typename FormLanguage::Traits<Lambda>::RangeType;
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

      void assemble(const Geometry::Polytope& polytope) override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& fe = m_fes.get().getFiniteElement(d, idx);
        const size_t n = fe.getCount();
        const size_t vc = Geometry::Polytope::getVertexCount(polytope.getGeometry());
        const size_t order = 2 * n + vc;
        const QF::QFGG qf(order, polytope.getGeometry());
        auto& res = getMatrix();
        res = Math::Matrix::Zero(n, n);
        for (size_t i = 0; i < qf.getSize(); i++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(i)));
          for (size_t local = 0; local < fe.getCount(); local++)
          {
            const auto basis = fe.getJacobian(local);
            const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
            const auto jacobian = p.getJacobianInverse().transpose() * basis(rc);
            const auto divergence = jacobian.diagonal();
            const auto symmetrized = jacobian + jacobian.transpose();
            res.noalias() += p.getDistortion() * getLambda().getValue(p) * divergence * divergence.transpose();
            res.noalias() += 0.5 * getMu().getValue(p) * symmetrized * symmetrized.transpose();
          }
        }
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
  };
}
#endif

