/*
 *          Copyright Carlos BRITO PACHECO 2026.
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

#include "Integral.h"

namespace Rodin::Variational
{
  /**
   * @ingroup LinearElasticitySpecializations
   * @brief Linear elasticity integrator for vector P1 spaces on a resident
   * Kokkos mesh.
   *
   * Assembles
   * @f[
   *   a(u, v) = \int_\Omega \lambda (\nabla \cdot u)(\nabla \cdot v)
   *     + \mu (\nabla u + \nabla u^T) : \nabla v \ dx
   * @f]
   * on the device as the sum of two device integrals, which is the bilinear
   * form of the host P1 specialization. Both integrals use the quadrature
   * order of the host,
   * @f$ \max(k_\lambda, k_\mu) + k_u + k_v @f$, with @f$ k_\lambda @f$ and
   * @f$ k_\mu @f$ the orders of the Lamé parameters, zero when they report
   * none, and @f$ k_u @f$, @f$ k_v @f$ the degrees of the trial and test
   * elements.
   *
   * The host evaluation interface (setPolytope(), integrate()) is not
   * available and raises; the form is assembled by Assembly::Kokkos.
   *
   * @tparam Solution Solution type of the trial function.
   * @tparam Backend Kokkos backend of the mesh context.
   * @tparam LambdaDerived Type of the first Lamé parameter.
   * @tparam MuDerived Type of the second Lamé parameter.
   */
  template <class Solution, class Backend, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator<Solution,
    Rodin::Kokkos::Variational::DeviceVectorP1<Backend>, LambdaDerived, MuDerived>
    final : public LocalBilinearFormIntegratorBase<Real>,
            public Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>
  {
    public:
      /// @brief Scalar type of the entries.
      using ScalarType = Real;

      /// @brief Finite element space of the trial and test functions.
      using FESType = Rodin::Kokkos::Variational::DeviceVectorP1<Backend>;

      /// @brief Trial function space type.
      using TrialFESType = FESType;

      /// @brief Test function space type.
      using TestFESType = FESType;

      /// @brief First Lamé parameter type.
      using LambdaType = FunctionBase<LambdaDerived>;

      /// @brief Second Lamé parameter (shear modulus) type.
      using MuType = FunctionBase<MuDerived>;

      /// @brief Trial function type.
      using TrialFunctionType = TrialFunction<Solution, FESType>;

      /// @brief Test function type.
      using TestFunctionType = TestFunction<FESType>;

      /// @brief Host integrator interface.
      using Parent = LocalBilinearFormIntegratorBase<Real>;

      /// @brief Device integrator interface.
      using DeviceParent =
        Rodin::Kokkos::Variational::DeviceBilinearFormIntegrator<Backend>;

      /// @brief Kokkos memory space.
      using MemorySpace = typename DeviceParent::MemorySpace;

      /// @brief Device integral of the dilatation term.
      using DilatationType = decltype(Integral(
        std::declval<const LambdaType&>() * Div(std::declval<const TrialFunctionType&>()),
        Div(std::declval<const TestFunctionType&>())));

      /// @brief Device integral of the shear term.
      using ShearType = decltype(Integral(std::declval<const MuType&>() *
          (Jacobian(std::declval<const TrialFunctionType&>()) +
            Transpose(Jacobian(std::declval<const TrialFunctionType&>()))),
        Jacobian(std::declval<const TestFunctionType&>())));

      /**
       * @brief Constructs the linear elasticity integrator.
       * @param u Trial function (displacement).
       * @param v Test function.
       * @param lambda First Lamé parameter.
       * @param mu Second Lamé parameter (shear modulus).
       */
      LinearElasticityIntegrator(const TrialFunctionType& u, const TestFunctionType& v,
        const LambdaType& lambda, const MuType& mu)
        : Parent(u, v),
          m_lambda(lambda.copy()),
          m_mu(mu.copy()),
          m_dilatation(Integral(lambda * Div(u), Div(v))),
          m_shear(Integral(mu * (Jacobian(u) + Transpose(Jacobian(u))), Jacobian(v)))
      {
        const auto order =
          [lambda = m_lambda, mu = m_mu, &trialfes = u.getFiniteElementSpace(),
            &testfes = v.getFiniteElementSpace()](const Geometry::Polytope& polytope) {
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
      LinearElasticityIntegrator(const LinearElasticityIntegrator& other) = default;

      /**
       * @brief Move constructor.
       * @param other Integrator to move from.
       */
      LinearElasticityIntegrator(LinearElasticityIntegrator&& other) = default;

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

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      const Geometry::Polytope& getPolytope() const override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      LinearElasticityIntegrator& setPolytope(
        [[maybe_unused]] const Geometry::Polytope& polytope) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Unavailable; the integrator has no host evaluation. @throws Alert::Exception Always.
      Real integrate([[maybe_unused]] size_t tr, [[maybe_unused]] size_t te) override
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Kokkos integrals are evaluated on the device by Assembly::Kokkos."
          << Alert::Raise;
        std::abort();
      }

      /// @brief Gets the arrangement of the entries. @return Entry layout.
      Rodin::Kokkos::Variational::DeviceLayout getLayout() const override
      {
        restrict();
        return m_dilatation.getLayout();
      }

      /**
       * @brief Locates the entries of every selected cell.
       * @return Number of entries, shared by both terms.
       */
      size_t prepare() const override
      {
        restrict();
        const size_t count = m_dilatation.prepare();
        [[maybe_unused]] const size_t shear = m_shear.prepare();
        assert(shear == count);
        return count;
      }

      /**
       * @brief Integrates both terms on the device.
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
      Geometry::Region getRegion() const override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Copies the integrator. @return Heap allocated copy.
      LinearElasticityIntegrator* copy() const noexcept override
      {
        return new LinearElasticityIntegrator(*this);
      }

    private:
      /// @brief Restricts both terms to the attributes of the integrator.
      void restrict() const
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
      mutable DilatationType m_dilatation;
      mutable ShearType m_shear;
  };
}

#endif
