/**
 * @file
 * @brief Linear elasticity bilinear form integrators.
 *
 * Implements the linear elasticity bilinear form with Lamé parameters
 * for solid mechanics problems. This generic implementation works with
 * any vector-valued finite element space.
 */
#ifndef RODIN_SOLID_LINEAR_LINEARELASTICITYINTEGRAL_H
#define RODIN_SOLID_LINEAR_LINEARELASTICITYINTEGRAL_H

#include "Rodin/QF/GenericPolytopeQuadrature.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Variational/Function.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
  /**
   * @brief Bilinear form integrator for linear elasticity.
   *
   * Assembles the bilinear form:
   * @f[
   *   a(u, v) = \int_\Omega [\lambda (\nabla \cdot u)(\nabla \cdot v) + 2\mu \varepsilon(u) : \varepsilon(v)] \, dx
   * @f]
   * where @f$ \varepsilon(u) = \frac{1}{2}(\nabla u + (\nabla u)^T) @f$ is the symmetric strain tensor.
   *
   * This generic implementation works for any vector-valued finite element
   * space by querying basis function derivatives through the standard FE
   * element interface. Quadrature order is derived automatically from the
   * polynomial order of the FE basis and the Lamé parameter functions.
   *
   * @tparam Solution Solution type
   * @tparam FES Finite element space type
   * @tparam LambdaDerived Type of first Lamé parameter function
   * @tparam MuDerived Type of second Lamé parameter (shear modulus) function
   */
  template <class Solution, class FES, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator final
    : public LocalBilinearFormIntegratorBase<typename FormLanguage::Traits<FES>::ScalarType>
  {
    public:
      /// Scalar type for computations
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;

      /// Trial function space type
      using TrialFESType = FES;

      /// Test function space type
      using TestFESType = FES;

      /// Second Lamé parameter (shear modulus) type
      using MuType = FunctionBase<MuDerived>;

      /// First Lamé parameter type
      using LambdaType = FunctionBase<LambdaDerived>;

      /// Parent class type
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

    private:
      /// Range type for shear modulus function
      using MuRangeType = typename FormLanguage::Traits<MuType>::RangeType;

      /// Range type for first Lamé parameter function
      using LambdaRangeType = typename FormLanguage::Traits<LambdaType>::RangeType;

      static_assert(std::is_same_v<MuRangeType, ScalarType>);

      static_assert(std::is_same_v<LambdaRangeType, ScalarType>);

    public:
      /**
       * @brief Constructs the linear elasticity integrator.
       * @param u Trial function (displacement)
       * @param v Test function
       * @param lambda First Lamé parameter function
       * @param mu Second Lamé parameter (shear modulus) function
       */
      LinearElasticityIntegrator(
          const TrialFunction<Solution, FES>& u, const TestFunction<FES>& v,
          const LambdaType& lambda, const MuType& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()), m_mu(mu.copy()),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace()),
          m_qf(nullptr),
          m_set(false),
          m_order(0),
          m_geometry(Geometry::Polytope::Type::Point)
      {
        assert(u.getFiniteElementSpace() == v.getFiniteElementSpace());
      }

      /**
       * @brief Copy constructor.
       * @param other Integrator to copy
       */
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

      /**
       * @brief Move constructor.
       * @param other Integrator to move
       */
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

      const Geometry::Polytope& getPolytope() const override
      {
        return m_polytope.value().get();
      }

      /**
       * @brief Sets the current polytope and computes the local stiffness matrix.
       *
       * Computes the local elasticity matrix over the given element using
       * quadrature derived from the FE basis order and Lamé parameter orders.
       * At each quadrature point, builds the physical-space Jacobian of each
       * basis function via the reference-to-physical map, then accumulates:
       * @f[
       *   A^K_{ij} = \sum_q w_q |J_q|
       *     \left[ \lambda \, (\nabla \cdot \phi_j)(\nabla \cdot \phi_i)
       *            + \mu \, (\nabla \phi_j + (\nabla \phi_j)^T) : (\nabla \phi_i + (\nabla \phi_i)^T) / 2 \right]
       * @f]
       *
       * When trial and test spaces coincide, exploits symmetry to halve work.
       *
       * @param polytope Current mesh element
       * @returns Reference to this integrator
       */
      LinearElasticityIntegrator& setPolytope(const Geometry::Polytope& polytope) override
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

      ScalarType integrate(size_t tr, size_t te) override
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

      /**
       * @brief Gets the integration region.
       * @returns Region::Cells for volume integration
       */
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /**
       * @brief Copies the integrator.
       * @returns Pointer to copied integrator
       */
      LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }

    private:
      std::unique_ptr<LambdaType> m_lambda;  ///< First Lamé parameter
      std::unique_ptr<MuType> m_mu;          ///< Shear modulus (second Lamé parameter)

      std::reference_wrapper<const TrialFESType> m_trialfes;  ///< Trial FE space
      std::reference_wrapper<const TestFESType> m_testfes;    ///< Test FE space

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;

      const QF::QuadratureFormulaBase* m_qf;  ///< Quadrature formula
      std::vector<Geometry::Point> m_ps;       ///< Quadrature points

      Math::Matrix<ScalarType> m_matrix;  ///< Local stiffness matrix

      bool m_set;                         ///< Whether QF is initialized
      size_t m_order;                     ///< Cached quadrature order
      Geometry::Polytope::Type m_geometry;  ///< Cached geometry type
  };

  /**
   * @brief Deduction guide for LinearElasticityIntegrator.
   */
  template <class Solution, class FES, class LambdaDerived, class MuDerived>
  LinearElasticityIntegrator(
      const TrialFunction<Solution, FES>&, const TestFunction<FES>&,
      const FunctionBase<LambdaDerived>&, const FunctionBase<MuDerived>&)
    -> LinearElasticityIntegrator<Solution, FES, LambdaDerived, MuDerived>;

  /**
   * @brief Helper class to construct linear elasticity integrators.
   *
   * This class provides a convenient interface for constructing linear elasticity
   * integrators with material parameters.
   *
   * @tparam Solution Solution type
   * @tparam FES Finite element space type
   */
  template <class Solution, class FES>
  class LinearElasticityIntegral final
  {
    public:
      /**
       * @brief Constructs the linear elasticity integral.
       * @param u Trial function (displacement)
       * @param v Test function
       */
      LinearElasticityIntegral(const TrialFunction<Solution, FES>& u, const TestFunction<FES>& v)
        : m_u(u), m_v(v)
      {}

      /**
       * @brief Creates integrator with constant Lamé parameters.
       * @param lambda First Lamé parameter (scalar constant)
       * @param mu Second Lamé parameter/shear modulus (scalar constant)
       * @returns LinearElasticityIntegrator with constant parameters
       */
      template <class L, class M>
      constexpr
      auto
      operator()(const L& lambda, const M& mu) const
      {
        return LinearElasticityIntegrator(m_u.get(), m_v.get(),
            RealFunction<L>(lambda), RealFunction<M>(mu));
      }

      /**
       * @brief Creates integrator with function-valued Lamé parameters.
       * @param lambda First Lamé parameter function
       * @param mu Second Lamé parameter/shear modulus function
       * @returns LinearElasticityIntegrator with spatially-varying parameters
       */
      template <class LambdaDerived, class MuDerived>
      constexpr
      auto
      operator()(const FunctionBase<LambdaDerived>& lambda, const FunctionBase<MuDerived>& mu) const
      {
        return LinearElasticityIntegrator(m_u.get(), m_v.get(), lambda, mu);
      }

    private:
      std::reference_wrapper<const TrialFunction<Solution, FES>> m_u;  ///< Trial function
      std::reference_wrapper<const TestFunction<FES>>  m_v;             ///< Test function
  };

  /**
   * @brief Deduction guide for LinearElasticityIntegral.
   */
  template <class Solution, class FES>
  LinearElasticityIntegral(const TrialFunction<Solution, FES>&, const TestFunction<FES>&)
    -> LinearElasticityIntegral<Solution, FES>;
}

#endif
