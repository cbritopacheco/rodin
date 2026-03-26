#ifndef RODIN_SOLID_LINEAR_P1_LINEARELASTICITYINTEGRAL_H
#define RODIN_SOLID_LINEAR_P1_LINEARELASTICITYINTEGRAL_H

/**
 * @file LinearElasticityIntegral.h
 * @brief P1 specialization of the linear elasticity bilinear form integrator.
 *
 * This file provides an optimized implementation of the linear elasticity
 * integrator for P1 vector finite element spaces.
 *
 * ## Bilinear Form
 * The integrator assembles the elasticity bilinear form:
 * @f[
 *   a(\mathbf{u}, \mathbf{v}) = \int_\Omega 2\mu \, \boldsymbol{\varepsilon}(\mathbf{u}) : \boldsymbol{\varepsilon}(\mathbf{v}) + \lambda (\nabla \cdot \mathbf{u})(\nabla \cdot \mathbf{v}) \, dx
 * @f]
 *
 * For P1 elements, this uses centroid quadrature (exact for the bilinear form)
 * and exploits the constant gradient property of P1 basis functions.
 *
 * ## Usage
 * @code{.cpp}
 * P1 Vh(mesh, spaceDim);  // Vector P1 space
 * TrialFunction u(Vh);
 * TestFunction v(Vh);
 *
 * Real mu = 1.0, lambda = 1.0;
 * Problem problem(u, v);
 * problem = Solid::Linear::LinearElasticityIntegral(u, v)(lambda, mu);
 * @endcode
 *
 * @see LinearElasticityIntegral, P1
 */

#include "Rodin/Solid/Linear/LinearElasticityIntegral.h"
#include "Rodin/Variational/P1/ForwardDecls.h"
#include "Rodin/Variational/P1/P1.h"
#include "Rodin/QF/Centroid.h"

namespace Rodin::Variational
{
  /**
   * @ingroup LinearElasticitySpecializations
   * @brief P1 linear elasticity bilinear form integrator.
   *
   * Specialized integrator for the linear elasticity problem with P1 vector
   * elements. Uses centroid quadrature which is exact for the constant
   * Jacobians of P1 basis functions.
   *
   * ## Local Stiffness Matrix
   * For element @f$ K @f$, computes the local matrix entries:
   * @f[
   *   A^K_{ij} = |K| \left[ 2\mu \, \boldsymbol{\varepsilon}(\phi_j) : \boldsymbol{\varepsilon}(\phi_i) + \lambda (\nabla \cdot \phi_j)(\nabla \cdot \phi_i) \right]
   * @f]
   *
   * @tparam Solution Solution type tag
   * @tparam MuDerived Derived type of shear modulus function
   * @tparam LambdaDerived Derived type of Lamé parameter function
   * @tparam Range Value range type
   * @tparam Mesh Mesh type
   */
  template <class Solution, class MuDerived, class LambdaDerived, class Range, class Mesh>
  class LinearElasticityIntegrator<Solution, P1<Range, Mesh>, MuDerived, LambdaDerived> final
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<P1<Range, Mesh>>::ScalarType>
  {
    public:
      using ScalarType = typename FormLanguage::Traits<P1<Range, Mesh>>::ScalarType;

      using TrialFESType = P1<Range, Mesh>;
      using TestFESType  = P1<Range, Mesh>;

      using MuType     = FunctionBase<MuDerived>;
      using LambdaType = FunctionBase<LambdaDerived>;

      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

    private:
      using MuRangeType =
        typename FormLanguage::Traits<MuType>::RangeType;

      using LambdaRangeType =
        typename FormLanguage::Traits<LambdaType>::RangeType;

      static_assert(std::is_same_v<Range, Math::Vector<ScalarType>>);

    public:
      LinearElasticityIntegrator(
          const TrialFunction<Solution, TrialFESType>& u,
          const TestFunction<TestFESType>& v,
          const LambdaType& lambda,
          const MuType& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()),
          m_mu(mu.copy()),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace()),
          m_qf(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      LinearElasticityIntegrator(const LinearElasticityIntegrator& other)
        : Parent(other),
          m_lambda(other.m_lambda ? other.m_lambda->copy() : nullptr),
          m_mu(other.m_mu ? other.m_mu->copy() : nullptr),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_polytope(other.m_polytope),
          m_qf(nullptr),
          m_ps(),
          m_matrix(),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {}

      LinearElasticityIntegrator(LinearElasticityIntegrator&& other)
        : Parent(std::move(other)),
          m_lambda(std::move(other.m_lambda)),
          m_mu(std::move(other.m_mu)),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_ps(std::move(other.m_ps)),
          m_matrix(std::move(other.m_matrix)),
          m_set(std::exchange(other.m_set, false)),
          m_order(std::exchange(other.m_order, 0)),
          m_geometry(std::exchange(other.m_geometry, Geometry::Polytope::Type::Point))
      {}

      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      LinearElasticityIntegrator& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d      = polytope.getDimension();
        const Index idx     = polytope.getIndex();
        const auto geometry = polytope.getGeometry();

        const auto& trialfes = m_trialfes.get();
        const auto& testfes  = m_testfes.get();

        const auto& trialfe = trialfes.getFiniteElement(d, idx);
        const auto& testfe  = testfes.getFiniteElement(d, idx);

        const size_t ntr = trialfe.getCount();
        const size_t nte = testfe.getCount();

        const size_t lambdaOrder =
          getLameFirstParameter().getOrder(polytope).value_or(size_t(0));

        const size_t muOrder =
          getShearModulus().getOrder(polytope).value_or(size_t(0));

        // Upper bound based on the actual finite elements on the current
        // polytope. This intentionally queries the FE order directly instead
        // of assuming anything from the abstract "P1" name.
        const size_t order =
          std::max(lambdaOrder, muOrder) + trialfe.getOrder() + testfe.getOrder();

        const bool recompute =
          !m_set || m_order != order || m_geometry != geometry;

        if (recompute)
        {
          m_set      = true;
          m_order    = order;
          m_geometry = geometry;

          m_qf = &QF::GenericPolytopeQuadrature::get(order, geometry);

          m_ps.clear();
          m_ps.reserve(m_qf->getSize());
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps.emplace_back(polytope, m_qf->getPoint(qp));
        }
        else
        {
          assert(m_qf);
          for (size_t qp = 0; qp < m_qf->getSize(); ++qp)
            m_ps[qp].setPolytope(polytope);
        }

        m_matrix.resize(
          static_cast<Eigen::Index>(nte),
          static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        const bool symmetric = (trialfes == testfes);

        if (symmetric)
        {
          for (size_t qp = 0; qp < m_ps.size(); ++qp)
          {
            const auto& p  = m_ps[qp];
            const auto& rc = m_qf->getPoint(qp);

            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const auto lambda = getLameFirstParameter().getValue(p);
            const auto mu     = getShearModulus().getValue(p);

            for (size_t i = 0; i < nte; ++i)
            {
              const auto& basis_i = testfe.getBasis(i);

              Math::SpatialMatrix<ScalarType> jac_i;
              jac_i.resize(d, d);
              for (size_t r = 0; r < d; ++r)
              {
                for (size_t c = 0; c < d; ++c)
                  jac_i(r, c) = basis_i.template getDerivative<1>(r, c)(rc);
              }
              jac_i *= p.getJacobianInverse();

              const auto sym_i = jac_i + jac_i.adjoint();
              const ScalarType div_i = jac_i.trace();

              m_matrix(i, i) +=
                  wdet
                * (
                    Math::dot(lambda * div_i, div_i)
                    + static_cast<ScalarType>(0.5) * Math::dot(mu * sym_i, sym_i)
                  );

              for (size_t j = 0; j < i; ++j)
              {
                const auto& basis_j = trialfe.getBasis(j);

                Math::SpatialMatrix<ScalarType> jac_j;
                jac_j.resize(d, d);
                for (size_t r = 0; r < d; ++r)
                {
                  for (size_t c = 0; c < d; ++c)
                    jac_j(r, c) = basis_j.template getDerivative<1>(r, c)(rc);
                }
                jac_j *= p.getJacobianInverse();

                const auto sym_j = jac_j + jac_j.adjoint();
                const ScalarType div_j = jac_j.trace();

                m_matrix(i, j) +=
                    wdet
                  * (
                      Math::dot(lambda * div_j, div_i)
                      + static_cast<ScalarType>(0.5) * Math::dot(mu * sym_j, sym_i)
                    );
              }
            }
          }

          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.adjoint();
        }
        else
        {
          for (size_t qp = 0; qp < m_ps.size(); ++qp)
          {
            const auto& p  = m_ps[qp];
            const auto& rc = m_qf->getPoint(qp);

            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const auto lambda = getLameFirstParameter().getValue(p);
            const auto mu     = getShearModulus().getValue(p);

            for (size_t i = 0; i < nte; ++i)
            {
              const auto& basis_i = testfe.getBasis(i);

              Math::SpatialMatrix<ScalarType> jac_i;
              jac_i.resize(d, d);
              for (size_t r = 0; r < d; ++r)
              {
                for (size_t c = 0; c < d; ++c)
                  jac_i(r, c) = basis_i.template getDerivative<1>(r, c)(rc);
              }
              jac_i *= p.getJacobianInverse();

              const auto sym_i = jac_i + jac_i.adjoint();
              const ScalarType div_i = jac_i.trace();

              for (size_t j = 0; j < ntr; ++j)
              {
                const auto& basis_j = trialfe.getBasis(j);

                Math::SpatialMatrix<ScalarType> jac_j;
                jac_j.resize(d, d);
                for (size_t r = 0; r < d; ++r)
                {
                  for (size_t c = 0; c < d; ++c)
                    jac_j(r, c) = basis_j.template getDerivative<1>(r, c)(rc);
                }
                jac_j *= p.getJacobianInverse();

                const auto sym_j = jac_j + jac_j.adjoint();
                const ScalarType div_j = jac_j.trace();

                m_matrix(i, j) +=
                    wdet
                  * (
                      Math::dot(lambda * div_j, div_i)
                      + static_cast<ScalarType>(0.5) * Math::dot(mu * sym_j, sym_i)
                    );
              }
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
      }

      /**
       * @brief Gets the shear modulus function.
       * @returns Reference to shear modulus (second Lamé parameter) function
       */
      constexpr
      const MuType& getShearModulus() const
      {
        assert(m_mu);
        return *m_mu;
      }

      /**
       * @brief Gets the first Lamé parameter function.
       * @returns Reference to first Lamé parameter function
       */
      constexpr
      const LambdaType& getLameFirstParameter() const
      {
        assert(m_lambda);
        return *m_lambda;
      }

      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }

    private:
      std::unique_ptr<LambdaType> m_lambda;
      std::unique_ptr<MuType> m_mu;

      std::reference_wrapper<const TrialFESType> m_trialfes;
      std::reference_wrapper<const TestFESType>  m_testfes;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;

      const QF::QuadratureFormulaBase* m_qf;
      std::vector<Geometry::Point> m_ps;

      Math::Matrix<ScalarType> m_matrix;

      bool m_set;
      size_t m_order;
      Geometry::Polytope::Type m_geometry;
  };
}
#endif
