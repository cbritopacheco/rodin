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
      using ScalarType = Real;

      using TrialFESType = P1<Range, Mesh>;

      using TestFESType = P1<Range, Mesh>;

      using MuType = FunctionBase<MuDerived>;

      using LambdaType = FunctionBase<LambdaDerived>;

      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

    private:
        using MuRangeType = typename FormLanguage::Traits<MuType>::RangeType;

        using LambdaRangeType = typename FormLanguage::Traits<LambdaType>::RangeType;

        static_assert(std::is_same_v<Range, Math::Vector<ScalarType>>);

        static_assert(std::is_same_v<MuRangeType, ScalarType>);

        static_assert(std::is_same_v<LambdaRangeType, ScalarType>);

    public:
      LinearElasticityIntegrator(
          const TrialFunction<TrialFESType>& u, const TestFunction<TestFESType>& v,
          const LambdaType& lambda, const MuType& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()), m_mu(mu.copy()),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace())
      {}

      LinearElasticityIntegrator(const LinearElasticityIntegrator& other)
        : Parent(other),
          m_lambda(other.m_lambda->copy()), m_mu(other.m_mu->copy()),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes)
      {}

      LinearElasticityIntegrator(LinearElasticityIntegrator&& other)
        : Parent(std::move(other)),
          m_lambda(std::move(other.m_lambda)), m_mu(std::move(other.mu)),
          m_trialfes(std::move(other.m_trialfes)),
          m_testfes(std::move(other.m_testfes))
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
        m_weight = m_qf->getWeight(0);
        m_distortion = m_p->getDistortion();
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trialfes = m_trialfes.get();
        const auto& testfes = m_testfes.get();
        const auto& rc = m_qf->getPoint(0);
        const auto& p = *m_p;
        const ScalarType mu = getMu().getValue(p);
        const ScalarType lambda = getLambda().getValue(p);
        if (trialfes == testfes)
        {
          const auto& fes = trialfes;
          const auto& fe = fes.getFiniteElement(d, idx);
          m_matrix.resize(fe.getCount(), fe.getCount());

          m_jac1.resize(fe.getCount());
          for (size_t i = 0; i < fe.getCount(); i++)
            fe.getJacobian(i)(m_jac1[i], rc);

          for (size_t i = 0; i < fe.getCount(); i++)
          {
            const auto jac = m_jac1[i] * p.getJacobianInverse();
            const auto sym = jac + jac.transpose();
            const ScalarType div = jac.trace();
            m_matrix(i, i) = lambda * div * div + 0.5 * mu * sym.squaredNorm();
          }

          for (size_t i = 0; i < fe.getCount(); i++)
          {
            const auto jac2 = m_jac1[i] * p.getJacobianInverse();
            const auto sym2 = jac2 + jac2.transpose();
            const ScalarType div2 = jac2.trace();
            for (size_t j = 0; j < i; j++)
            {
              const auto jac1 = m_jac1[j] * p.getJacobianInverse();
              const auto sym1 = jac1 + jac1.transpose();
              const ScalarType div1 = jac1.trace();
              m_matrix(i, j) =
                lambda * div1 * div2 + ((mu * sym2).array() * (0.5 * sym1).array()).rowwise().sum().colwise().sum().value();
            }
          }
          m_matrix.template triangularView<Eigen::Upper>() = m_matrix.transpose();
        }
        else
        {
          assert(false);
        }
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        return m_weight * m_distortion * m_matrix(te, tr);
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
      std::reference_wrapper<const TrialFESType> m_trialfes;
      std::reference_wrapper<const TestFESType> m_testfes;

      std::optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      std::optional<QF::QF1P1> m_qf;
      std::optional<Geometry::Point> m_p;

      Real m_distortion;
      Real m_weight;

      std::vector<Math::SpatialMatrix<ScalarType>> m_jac1, m_jac2;

      Math::Matrix<ScalarType> m_matrix;
  };
}
#endif

