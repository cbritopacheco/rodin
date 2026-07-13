/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
#include "Rodin/Math/Traits.h"

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
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<P1<Range, Mesh>>::ScalarType;

      /// @brief Trial finite element space type.
      using TrialFESType = P1<Range, Mesh>;
      /// @brief Test finite element space type.
      using TestFESType  = P1<Range, Mesh>;

      /// @brief Shear modulus function type.
      using MuType     = FunctionBase<MuDerived>;
      /// @brief First Lamé parameter function type.
      using LambdaType = FunctionBase<LambdaDerived>;

      /// @brief Parent class type.
      using Parent =
        LocalBilinearFormIntegratorBase<ScalarType>;

    private:
      using MuRangeType =
        typename FormLanguage::Traits<MuType>::RangeType;

      using LambdaRangeType =
        typename FormLanguage::Traits<LambdaType>::RangeType;

      static_assert(FormLanguage::IsVectorRange<Range>::Value);

    public:
      /// @brief Constructs the P1 linear elasticity integrator.
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
          m_quadrature(nullptr)
      {}

      /// @brief Copy constructor.
      LinearElasticityIntegrator(const LinearElasticityIntegrator& other)
        : Parent(other),
          m_lambda(other.m_lambda ? other.m_lambda->copy() : nullptr),
          m_mu(other.m_mu ? other.m_mu->copy() : nullptr),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_polytope(other.m_polytope),
          m_qf(nullptr),
          m_quadrature(nullptr),
          m_matrix()
      {}

      /// @brief Move constructor.
      LinearElasticityIntegrator(LinearElasticityIntegrator&& other)
        : Parent(std::move(other)),
          m_lambda(std::move(other.m_lambda)),
          m_mu(std::move(other.m_mu)),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_polytope(std::move(other.m_polytope)),
          m_qf(std::exchange(other.m_qf, nullptr)),
          m_quadrature(std::exchange(other.m_quadrature, nullptr)),
          m_matrix(std::move(other.m_matrix))
      {}

      /// @brief Returns the current polytope.
      const Geometry::Polytope& getPolytope() const final override
      {
        return m_polytope.value().get();
      }

      /// @brief Sets the current polytope and assembles the local matrix.
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

        const size_t order =
          std::max(lambdaOrder, muOrder) + trialfe.getOrder() + testfe.getOrder();

        m_qf = &QF::PolytopeQuadratureFormula::get(order, geometry);
        m_quadrature = &polytope.getQuadrature(*m_qf);

        m_matrix.resize(
          static_cast<Eigen::Index>(nte),
          static_cast<Eigen::Index>(ntr));
        m_matrix.setZero();

        const bool symmetric = (trialfes == testfes);

        if (symmetric)
        {
          for (size_t qp = 0; qp < m_quadrature->getSize(); ++qp)
          {
            const auto& p  = m_quadrature->getPoint(qp);
            const Variational::IntegrationPoint ip(p, m_qf, qp);
            const auto& rc = m_qf->getPoint(qp);

            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const auto lambda = getLameFirstParameter().getValue(ip);
            const auto mu     = getShearModulus().getValue(ip);

            for (size_t i = 0; i < nte; ++i)
            {
              const auto& basisI = testfe.getBasis(i);

              Math::SpatialMatrix<ScalarType> jacI;
              jacI.resize(d, d);
              for (size_t r = 0; r < d; ++r)
              {
                for (size_t c = 0; c < d; ++c)
                  jacI(r, c) = basisI.template getDerivative<1>(r, c)(rc);
              }
              jacI *= p.getJacobianInverse();

              const auto symI = jacI + jacI.adjoint();
              const ScalarType divI = jacI.trace();

              m_matrix(i, i) += wdet *
                (Math::dot(lambda * divI, divI) +
                  static_cast<ScalarType>(0.5) * Math::dot(mu * symI, symI));

              for (size_t j = 0; j < i; ++j)
              {
                const auto& basisJ = trialfe.getBasis(j);

                Math::SpatialMatrix<ScalarType> jacJ;
                jacJ.resize(d, d);
                for (size_t r = 0; r < d; ++r)
                {
                  for (size_t c = 0; c < d; ++c)
                    jacJ(r, c) = basisJ.template getDerivative<1>(r, c)(rc);
                }
                jacJ *= p.getJacobianInverse();

                const auto symJ = jacJ + jacJ.adjoint();
                const ScalarType divJ = jacJ.trace();

                m_matrix(i, j) += wdet *
                  (Math::dot(lambda * divJ, divI) +
                    static_cast<ScalarType>(0.5) * Math::dot(mu * symJ, symI));
              }
            }
          }

          m_matrix.template triangularView<Eigen::Upper>() =
            m_matrix.adjoint();
        }
        else
        {
          for (size_t qp = 0; qp < m_quadrature->getSize(); ++qp)
          {
            const auto& p  = m_quadrature->getPoint(qp);
            const Variational::IntegrationPoint ip(p, m_qf, qp);
            const auto& rc = m_qf->getPoint(qp);

            const ScalarType wdet =
              static_cast<ScalarType>(m_qf->getWeight(qp) * p.getDistortion());

            const auto lambda = getLameFirstParameter().getValue(ip);
            const auto mu     = getShearModulus().getValue(ip);

            for (size_t i = 0; i < nte; ++i)
            {
              const auto& basisI = testfe.getBasis(i);

              Math::SpatialMatrix<ScalarType> jacI;
              jacI.resize(d, d);
              for (size_t r = 0; r < d; ++r)
              {
                for (size_t c = 0; c < d; ++c)
                  jacI(r, c) = basisI.template getDerivative<1>(r, c)(rc);
              }
              jacI *= p.getJacobianInverse();

              const auto symI = jacI + jacI.adjoint();
              const ScalarType divI = jacI.trace();

              for (size_t j = 0; j < ntr; ++j)
              {
                const auto& basisJ = trialfe.getBasis(j);

                Math::SpatialMatrix<ScalarType> jacJ;
                jacJ.resize(d, d);
                for (size_t r = 0; r < d; ++r)
                {
                  for (size_t c = 0; c < d; ++c)
                    jacJ(r, c) = basisJ.template getDerivative<1>(r, c)(rc);
                }
                jacJ *= p.getJacobianInverse();

                const auto symJ = jacJ + jacJ.adjoint();
                const ScalarType divJ = jacJ.trace();

                m_matrix(i, j) += wdet *
                  (Math::dot(lambda * divJ, divI) +
                    static_cast<ScalarType>(0.5) * Math::dot(mu * symJ, symI));
              }
            }
          }
        }

        return *this;
      }

      /// @brief Returns an entry of the current element stiffness matrix.
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

      /// @brief Returns the integration region.
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
      const Geometry::PolytopeQuadrature* m_quadrature;

      Math::Matrix<ScalarType> m_matrix;
  };
}
#endif
