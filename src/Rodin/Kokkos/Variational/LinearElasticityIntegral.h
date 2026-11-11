/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_LINEARELASTICITYINTEGRAL_H
#define RODIN_KOKKOS_VARIATIONAL_LINEARELASTICITYINTEGRAL_H

/**
 * @file
 * @brief Device linear elasticity integrator on a resident Kokkos mesh.
 */

#include <algorithm>
#include <memory>

#include "Rodin/Solid/Linear/LinearElasticityIntegral.h"
#include "Rodin/Variational/Div.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/Transpose.h"

#include "H1.h"
#include "Integral.h"

namespace Rodin::Kokkos::Variational
{
  /**
   * @brief Linear elasticity integrator for vector spaces on a resident
   * Kokkos mesh.
   *
   * Assembles
   * @f[
   *   a(u, v) = \int_\Omega \lambda (\nabla \cdot u)(\nabla \cdot v)
   *     + \mu (\nabla u + \nabla u^T) : \nabla v \ dx
   * @f]
   * on the device as the sum of two device integrals, which is the bilinear
   * form of the host integrators. Both integrals use the quadrature order of
   * the host,
   * @f$ \max(k_\lambda, k_\mu) + k_u + k_v @f$, with @f$ k_\lambda @f$ and
   * @f$ k_\mu @f$ the orders of the Lamé parameters, zero when they report
   * none, and @f$ k_u @f$, @f$ k_v @f$ the degrees of the trial and test
   * elements.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos. The
   * specializations of Rodin::Variational::LinearElasticityIntegrator for the
   * vector P1 and H1 spaces of a resident mesh derive from this class.
   *
   * @tparam Solution Solution type of the trial function.
   * @tparam FES Vector device space of the trial and test functions.
   * @tparam LambdaDerived Type of the first Lamé parameter.
   * @tparam MuDerived Type of the second Lamé parameter.
   */
  template <class Solution, VectorDeviceSpace FES, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegratorBase
    : public ::Rodin::Variational::LocalBilinearFormIntegratorBase<Real>,
      public DeviceBilinearFormIntegrator<Internal::BackendOf<FES>>
  {
    public:
      /// @brief Scalar type of the entries.
      using ScalarType = Real;

      /// @brief Finite element space of the trial and test functions.
      using FESType = FES;

      /// @brief Trial function space type.
      using TrialFESType = FESType;

      /// @brief Test function space type.
      using TestFESType = FESType;

      /// @brief First Lamé parameter type.
      using LambdaType = ::Rodin::Variational::FunctionBase<LambdaDerived>;

      /// @brief Second Lamé parameter (shear modulus) type.
      using MuType = ::Rodin::Variational::FunctionBase<MuDerived>;

      /// @brief Trial function type.
      using TrialFunctionType = ::Rodin::Variational::TrialFunction<Solution, FESType>;

      /// @brief Test function type.
      using TestFunctionType = ::Rodin::Variational::TestFunction<FESType>;

      /// @brief Host integrator interface.
      using Parent = ::Rodin::Variational::LocalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent = DeviceBilinearFormIntegrator<Internal::BackendOf<FES>>;

      /// @brief Kokkos memory space.
      using MemorySpace = typename DeviceParent::MemorySpace;

      /// @brief Device integral of the dilatation term.
      using DilatationType =
        decltype(::Rodin::Variational::Integral(std::declval<const LambdaType&>() *
            ::Rodin::Variational::Div(std::declval<const TrialFunctionType&>()),
          ::Rodin::Variational::Div(std::declval<const TestFunctionType&>())));

      /// @brief Device integral of the shear term.
      using ShearType =
        decltype(::Rodin::Variational::Integral(std::declval<const MuType&>() *
            (::Rodin::Variational::Jacobian(std::declval<const TrialFunctionType&>()) +
              ::Rodin::Variational::Transpose(::Rodin::Variational::Jacobian(
                std::declval<const TrialFunctionType&>()))),
          ::Rodin::Variational::Jacobian(std::declval<const TestFunctionType&>())));

      /**
       * @brief Constructs the linear elasticity integrator.
       * @param u Trial function (displacement).
       * @param v Test function.
       * @param lambda First Lamé parameter.
       * @param mu Second Lamé parameter (shear modulus).
       */
      LinearElasticityIntegratorBase(const TrialFunctionType& u,
        const TestFunctionType& v, const LambdaType& lambda, const MuType& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()),
          m_mu(mu.copy()),
          m_dilatation(::Rodin::Variational::Integral(
            lambda * ::Rodin::Variational::Div(u), ::Rodin::Variational::Div(v))),
          m_shear(::Rodin::Variational::Integral(mu *
              (::Rodin::Variational::Jacobian(u) +
                ::Rodin::Variational::Transpose(::Rodin::Variational::Jacobian(u))),
            ::Rodin::Variational::Jacobian(v)))
      {
        const auto order =
          [lambda = m_lambda, mu = m_mu, &trialfes = u.getFiniteElementSpace(),
            &testfes = v.getFiniteElementSpace()](
            const ::Rodin::Geometry::Polytope& polytope) {
            const size_t d = polytope.getDimension();
            const Index i = polytope.getIndex();
            const auto trialDegree = trialfes.getFiniteElement(d, i).getDegree();
            const auto testDegree = testfes.getFiniteElement(d, i).getDegree();
            if (!trialDegree || !testDegree)
            {
              Alert::Exception()
                << "LinearElasticityIntegrator requires finite elements of known degree."
                << Alert::Raise;
            }
            return std::max(lambda->getOrder(polytope).value_or(size_t(0)),
                     mu->getOrder(polytope).value_or(size_t(0))) +
              *trialDegree + *testDegree;
          };
        m_dilatation.setOrder(order);
        m_shear.setOrder(order);
      }

      /**
       * @brief Copy constructor.
       * @param other Integrator to copy.
       */
      LinearElasticityIntegratorBase(
        const LinearElasticityIntegratorBase& other) = default;

      /**
       * @brief Move constructor.
       * @param other Integrator to move from.
       */
      LinearElasticityIntegratorBase(LinearElasticityIntegratorBase&& other) = default;

      /// @brief Gets the first Lamé parameter. @return First Lamé parameter.
      const LambdaType& getLameFirstParameter() const
      {
        return *m_lambda;
      }

      /// @brief Gets the shear modulus. @return Second Lamé parameter.
      const MuType& getShearModulus() const
      {
        return *m_mu;
      }

      /**
       * @brief Unavailable; the integrator has no host evaluation.
       * @return Never returns; the call aborts after raising.
       * @throws Alert::Exception Always.
       */
      const ::Rodin::Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Unavailable; the integrator has no host evaluation.
       * @param polytope Polytope that would be integrated, ignored.
       * @return Never returns; the call aborts after raising.
       * @throws Alert::Exception Always.
       */
      LinearElasticityIntegratorBase& setPolytope(
        [[maybe_unused]] const ::Rodin::Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /**
       * @brief Unavailable; the integrator has no host evaluation.
       * @param tr Local trial index, ignored.
       * @param te Local test index, ignored.
       * @return Never returns; the call aborts after raising.
       * @throws Alert::Exception Always.
       */
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      DeviceLayout getLayout() const override
      {
        auto layout = m_dilatation.getLayout();
        const auto& attributes = this->getAttributes();
        layout.attributes.assign(attributes.begin(), attributes.end());
        return layout;
      }

      /**
       * @brief Locates the entries of every selected cell.
       * @param diagonal Whether to list only the diagonals of the component
       * blocks, ignored since both terms couple every pair of components.
       * @return Number of entries, shared by both terms.
       */
      size_t prepare([[maybe_unused]] bool diagonal) override
      {
        restrict();
        const size_t count = m_dilatation.prepare(false);
        [[maybe_unused]] const size_t shear = m_shear.prepare(false);
        assert(shear == count);
        return count;
      }

      /// @brief Gets the spacing of the entries, shared by both terms. @return Spacing.
      size_t getSpacing() const override
      {
        assert(m_shear.getSpacing() == m_dilatation.getSpacing());
        return m_dilatation.getSpacing();
      }

      /**
       * @brief Integrates both terms on the device.
       * @param values Entry values, the entry of rank r being written at
       * offset + getEntryPosition(r, count, getSpacing()).
       * @param offset Position of the first entry.
       * @param sign Factor of the contributions.
       * @param accumulate Whether the contributions of the dilatation term
       * are added to the values instead of replacing them.
       * @param rows Row indices, test degrees of freedom, written as the
       * values unless empty.
       * @param cols Column indices, trial degrees of freedom, written as the
       * values unless empty.
       *
       * The shear term accumulates onto the entries of the dilatation term,
       * which share their arrangement.
       */
      void fill(const ::Kokkos::View<Real*, MemorySpace>& values, size_t offset,
        Real sign, bool accumulate, const ::Kokkos::View<Index*, MemorySpace>& rows,
        const ::Kokkos::View<Index*, MemorySpace>& cols) const override
      {
        m_dilatation.fill(values, offset, sign, accumulate, rows, cols);
        m_shear.fill(values, offset, sign, true, {}, {});
      }

      /// @brief Gets the integration region. @return Cells.
      ::Rodin::Geometry::Region getRegion() const override
      {
        return ::Rodin::Geometry::Region::Cells;
      }

    private:
      /// @brief Restricts both terms to the attributes of the integrator.
      void restrict()
      {
        const auto& attributes = this->getAttributes();
        if (!attributes.empty())
        {
          m_dilatation.over(attributes);
          m_shear.over(attributes);
        }
      }

      std::shared_ptr<const LambdaType> m_lambda;
      std::shared_ptr<const MuType> m_mu;
      DilatationType m_dilatation;
      ShearType m_shear;
  };
}

namespace Rodin::Variational
{
  /**
   * @ingroup LinearElasticitySpecializations
   * @brief Linear elasticity integrator for vector P1 spaces on a resident
   * Kokkos mesh.
   * @see Rodin::Kokkos::Variational::LinearElasticityIntegratorBase
   */
  template <class Solution, class Backend, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator<Solution,
    Rodin::Kokkos::Variational::DeviceVectorP1<Backend>, LambdaDerived, MuDerived>
    final
    : public Rodin::Kokkos::Variational::LinearElasticityIntegratorBase<Solution,
        Rodin::Kokkos::Variational::DeviceVectorP1<Backend>, LambdaDerived, MuDerived>
  {
    public:
      /// @brief Device integrator.
      using Parent = Rodin::Kokkos::Variational::LinearElasticityIntegratorBase<Solution,
        Rodin::Kokkos::Variational::DeviceVectorP1<Backend>, LambdaDerived, MuDerived>;

      using Parent::Parent;

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      LinearElasticityIntegrator& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }
  };

  /**
   * @ingroup LinearElasticitySpecializations
   * @brief Linear elasticity integrator for vector H1 spaces on a resident
   * Kokkos mesh.
   * @see Rodin::Kokkos::Variational::LinearElasticityIntegratorBase
   */
  template <class Solution, size_t K, class Backend, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator<Solution,
    H1<K, Math::SpatialVector<Real>,
      Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>,
    LambdaDerived, MuDerived>
    final : public Rodin::Kokkos::Variational::LinearElasticityIntegratorBase<Solution,
              H1<K, Math::SpatialVector<Real>,
                Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>,
              LambdaDerived, MuDerived>
  {
    public:
      /// @brief Device integrator.
      using Parent = Rodin::Kokkos::Variational::LinearElasticityIntegratorBase<Solution,
        H1<K, Math::SpatialVector<Real>,
          Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>,
        LambdaDerived, MuDerived>;

      using Parent::Parent;

      /**
       * @brief Restricts the integration to mesh attributes, as the parent
       * over().
       * @param[in] attrs Attribute, attributes or set of attributes of the
       * integrated polytopes.
       * @return Reference to this integrator, which keeps its device type.
       */
      template <class... Attributes>
      LinearElasticityIntegrator& over(const Attributes&... attrs)
      {
        Parent::over(attrs...);
        return *this;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }
  };
}

#endif
