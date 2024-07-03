#ifndef RODIN_VARIATIONAL_P1_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_P1_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

#include "Rodin/Variational/LinearElasticity/LinearElasticityIntegral.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"
#include "Rodin/Variational/P1/P1Element.h"
#include "Rodin/QF/QF1P1.h"

namespace Rodin::Variational
{
  template <class MuDerived, class LambdaDerived, class Range, class Mesh>
  class LinearElasticityIntegrator<P1<Range, Mesh>, MuDerived, LambdaDerived> final
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<P1<Range, Mesh>>::ScalarType>
  {
    public:
      using FESType = P1<Range, Mesh>;

      using MuType = FunctionBase<MuDerived>;

      using LambdaType = FunctionBase<LambdaDerived>;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

    private:
        using MuRangeType = typename FormLanguage::Traits<MuType>::RangeType;

        using LambdaRangeType = typename FormLanguage::Traits<LambdaType>::RangeType;

        static_assert(std::is_same_v<Range, Math::Vector<ScalarType>>);

        static_assert(std::is_same_v<MuRangeType, ScalarType>);

        static_assert(std::is_same_v<LambdaRangeType, ScalarType>);

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

      inline
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      LinearElasticityIntegrator& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
        const auto& trans = polytope.getTransformation();
        m_qf.emplace(polytope.getGeometry());
        assert(m_qf->getSize() == 1);
        m_p.emplace(polytope, trans, std::cref(m_qf->getPoint(0)));
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        const auto& polytope = getPolytope();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& fes = getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const auto& qf = *m_qf;
        assert(qf.getSize() == 1);
        const auto& p = m_p.value();
        const auto& w = qf.getWeight(0);
        const auto& rc = qf.getPoint(0);
        const ScalarType mu = getMu().getValue(p);
        const ScalarType lambda = getLambda().getValue(p);
        if (tr == te)
        {
          fe.getJacobian(tr)(m_jac1, rc);
          const auto jac = m_jac1 * p.getJacobianInverse();
          const auto sym = jac + jac.transpose();
          const ScalarType div = jac.trace();
          return w * p.getDistortion() * (
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
          return w * p.getDistortion() * (
              getLambda().getValue(p) * div1 * div2 + (
                (mu * sym2).array() * (0.5 * sym1).array()).rowwise().sum().colwise().sum().value());
        }
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

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Math::SpatialMatrix<Real> m_jac1, m_jac2;
  };
}
#endif

