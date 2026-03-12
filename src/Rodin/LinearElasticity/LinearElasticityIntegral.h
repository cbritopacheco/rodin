/**
 * @file
 * @brief Linear elasticity bilinear form integrators.
 *
 * Implements the linear elasticity bilinear form with Lamé parameters
 * for solid mechanics problems.
 */
#ifndef RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H
#define RODIN_VARIATIONAL_LINEARELASTICITY_LINEARELASTICITYINTEGRAL_H

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
      using ScalarType = Real;

      /// Trial function space type
      using TrialFESType = FES;

      /// Test function space type
      using TestFESType = FES;

      /// Second Lamé parameter (shear modulus) type
      using Mu = FunctionBase<MuDerived>;

      /// First Lamé parameter type  
      using Lambda = FunctionBase<LambdaDerived>;

      /// Parent class type
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

    private:
        /// Range type for shear modulus function
        using MuRangeType = typename FormLanguage::Traits<Mu>::RangeType;

        /// Range type for first Lamé parameter function
        using LambdaRangeType = typename FormLanguage::Traits<Lambda>::RangeType;

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
          const Lambda& lambda, const Mu& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()), m_mu(mu.copy()),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace())
      {
        assert(u.getFiniteElementSpace() == v.getFiniteElementSpace());
      }

      /**
       * @brief Copy constructor.
       * @param other Integrator to copy
       */
      LinearElasticityIntegrator(const LinearElasticityIntegrator& other)
        : Parent(other),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes)
      {}

      /**
       * @brief Move constructor.
       * @param other Integrator to move
       */
      LinearElasticityIntegrator(LinearElasticityIntegrator&& other)
        : Parent(std::move(other)),
          m_lambda(std::move(other.m_lambda)), m_mu(std::move(other.mu)),
          m_trialfes(std::move(other.m_trialfes)),
          m_testfes(std::move(other.m_testfes))
      {}

      const Geometry::Polytope& getPolytope() const override
      {
        return m_polytope.value().get();
      }

      LinearElasticityIntegrator& setPolytope(const Geometry::Polytope& polytope) override
      {
        assert(false);
        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) override
      {
        assert(false);
        return NAN;
      }

      /**
       * @brief Gets the shear modulus function.
       * @returns Reference to shear modulus (second Lamé parameter) function
       */
      constexpr
      const Mu& getMu() const
      {
        assert(m_mu);
        return *m_mu;
      }

      /**
       * @brief Gets the first Lamé parameter function.
       * @returns Reference to first Lamé parameter function
       */
      constexpr
      const Lambda& getLambda() const
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
      std::unique_ptr<Lambda> m_lambda;  ///< First Lamé parameter
      std::unique_ptr<Mu> m_mu;          ///< Shear modulus (second Lamé parameter)

      std::reference_wrapper<const TrialFESType> m_trialfes;  ///< Trial FE space
      std::reference_wrapper<const TestFESType> m_testfes;    ///< Test FE space

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;  ///< Current element
      std::unique_ptr<QF::QuadratureFormulaBase> m_qf;  ///< Quadrature formula
      std::vector<Geometry::Point> m_ps;  ///< Quadrature points

      Math::SpatialMatrix<Real> m_jac1, m_jac2;  ///< Jacobian matrices
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

