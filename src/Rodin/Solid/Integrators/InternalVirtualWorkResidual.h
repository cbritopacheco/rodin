/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalVirtualWorkResidual.h
 * @brief Residual integrator for the internal virtual work in hyperelastic formulations.
 *
 * Assembles the nonlinear residual contribution:
 * @f[
 *   \delta W^{\text{int}}(\mathbf{v}) = \int_{\Omega_0} \mathbf{P}(\mathbf{u}) : \nabla_0 \mathbf{v} \, dX
 * @f]
 * where @f$ \mathbf{P} @f$ is the first Piola-Kirchhoff stress and
 * @f$ \mathbf{v} @f$ is the test function (virtual displacement).
 *
 * The integrator is generic: it obtains the finite element basis from the
 * FE space (not hardcoded to P1), supports arbitrary quadrature rules, and
 * builds a ConstitutivePoint (composed over Geometry::Point) at each
 * quadrature point for constitutive evaluation.
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORKRESIDUAL_H
#define RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORKRESIDUAL_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <type_traits>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/LinearFormIntegrator.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Geometry/Point.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"
#include "Rodin/Solid/Local/Input.h"
#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Linear form integrator for the internal virtual work residual in
   * hyperelastic problems.
   *
   * Computes the element-level contribution:
   * @f[
   *   (\delta W^{\text{int}})^e_i = \int_{K} \mathbf{P} : \nabla_0 \phi_i \, dX
   * @f]
   *
   * The test function and current displacement state are kept as their
   * concrete Rodin types. The displacement gradient and test-basis gradients
   * are evaluated through @c Variational::Jacobian(...), so PETSc or other
   * storage backends work as long as they provide the usual Rodin traits and
   * Jacobian specialization. An optional Input can inject auxiliary data
   * (fiber directions, activation, etc.) into the ConstitutivePoint at each
   * quadrature point.
   *
   * Typically constructed via @c InternalVirtualWork::Residual(v) rather
   * than directly.
   *
   * @tparam LawDerived The hyperelastic constitutive law type
   * @tparam TestFunctionType The test function type (backend-generic)
   * @tparam DisplacementType The current displacement grid-function type
   */
  template <class LawDerived, class TestFunctionType, class DisplacementType>
  class InternalVirtualWorkResidual final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Parent class type.
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using TestType = TestFunctionType;
      using StateType = DisplacementType;

      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::InternalVirtualWorkResidual expects a Rodin test function.");

      /**
       * @brief Constructs the internal virtual work residual integrator.
       * @param law The constitutive law (stored by value)
       * @param v The test function
       * @param displacement Current displacement grid function
       */
      InternalVirtualWorkResidual(
          const LawDerived& law,
          const TestType& v,
          const StateType& displacement)
        : Parent(v),
          m_law(law),
          m_test(v),
          m_displacement(displacement),
          m_testfes(v.getFiniteElementSpace()),
          m_statefes(displacement.getFiniteElementSpace()),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      InternalVirtualWorkResidual(const InternalVirtualWorkResidual& other)
        : Parent(other),
          m_law(other.m_law),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_testfes(other.m_testfes),
          m_statefes(other.m_statefes),
          m_quadOrder(other.m_quadOrder),
          m_input(other.m_input)      {}

      /**
       * @brief Rebinds the current displacement state.
       * @param gf Displacement grid function with the same concrete state type
       * @returns Reference to this object for chaining
       */
      InternalVirtualWorkResidual& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        m_statefes = std::cref(gf.getFiniteElementSpace());
        return *this;
      }

      /**
       * @brief Sets the quadrature order.
       *
       * When set to 0 (the default), the quadrature order is determined
       * automatically from the finite element approximation order as
       * @c 2 * fe.getOrder().
       *
       * @param order Polynomial order for exact integration (0 = auto)
       * @returns Reference to this object for chaining
       */
      InternalVirtualWorkResidual& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      /**
       * @brief Sets an input for auxiliary constitutive data.
       *
       * The input is called at each quadrature point after the
       * ConstitutivePoint has been constructed with geometric context and
       * kinematics, allowing injection of fiber directions, activation
       * parameters, region-wise material properties, etc.
       *
       * @param input A callable with signature void(ConstitutivePoint&)
       * @returns Reference to this object for chaining
       */
      InternalVirtualWorkResidual& setInput(InputFunction input)
      {
        m_input = std::move(input);
        return *this;
      }

      InternalVirtualWorkResidual& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& testFES = m_testfes.get();
        const auto& stateFES = m_statefes.get();
        const size_t vdim = testFES.getVectorDimension();

        const auto& testFE = testFES.getFiniteElement(d, idx);
        const auto& stateFE = stateFES.getFiniteElement(d, idx);
        const size_t testDofs = testFE.getCount();

        // Determine effective quadrature order
        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * std::max(testFE.getOrder(), stateFE.getOrder());
        const auto& qf = QF::PolytopeQuadratureFormula::get(effectiveOrder, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        // Zero element vector
        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        auto stateGradient = Variational::Jacobian(m_displacement.get());
        auto testGradient = Variational::Jacobian(m_test.get());

        // Loop over quadrature points
        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const ScalarType wq = qf.getWeight(q);

          const ScalarType distortion = pt.getDistortion();

          testGradient.setIntegrationPoint(ip);
          const auto H = stateGradient.getValue(ip);

          // Build ConstitutivePoint composed over Geometry::Point
          KinematicState state(d);
          state.setDisplacementGradient(H);

          ConstitutivePoint cp(pt, state);
          cp.set<Tags::CellIndex>(idx);
          cp.set<Tags::QuadraturePointIndex>(q);

          // Invoke input for auxiliary data injection
          if (m_input)
            m_input(cp);

          typename LawType::Cache cache;
          m_law.setCache(cache, cp);

          Math::SpatialMatrix<ScalarType> P;
          m_law.getFirstPiolaKirchhoffStress(P, cache, cp);

          // Accumulate into element vector: R_te = P : grad psi_te
          for (size_t te = 0; te < testDofs; ++te)
          {
            const auto gradTest = testGradient.getBasis(te);
            ScalarType val = 0;
            for (size_t c = 0; c < vdim; ++c)
              for (size_t k = 0; k < d; ++k)
                val += P(c, k) * gradTest(c, k);
            m_elemVec(te) += wq * distortion * val;
          }
        }

        return *this;
      }

      ScalarType integrate(size_t te) final override
      {
        return m_elemVec(te);
      }

      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return m_polytope->get();
      }

      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      InternalVirtualWorkResidual* copy() const noexcept final override
      {
        return new InternalVirtualWorkResidual(*this);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      void checkCompatibility(const StateType& displacement) const
      {
        const auto& testFES = m_testfes.get();
        const auto& stateFES = displacement.getFiniteElementSpace();

        assert(&stateFES.getMesh() == &testFES.getMesh());
        assert(stateFES.getVectorDimension() == testFES.getVectorDimension());
        (void) testFES;
        (void) stateFES;
      }

      LawType m_law;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      std::reference_wrapper<const TestFESType> m_testfes;
      std::reference_wrapper<const StateFESType> m_statefes;
      size_t m_quadOrder;
      InputFunction m_input;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /// CTAD deduction guide for InternalVirtualWorkResidual
  template <class LawDerived, class TestFunctionType, class DisplacementType>
  InternalVirtualWorkResidual(const LawDerived&,
                               const TestFunctionType&,
                               const DisplacementType&)
    -> InternalVirtualWorkResidual<
         LawDerived,
         std::decay_t<TestFunctionType>,
         std::decay_t<DisplacementType>>;
}

#endif
