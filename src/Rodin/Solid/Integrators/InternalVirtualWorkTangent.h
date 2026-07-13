/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalVirtualWorkTangent.h
 * @brief Tangent integrator for the internal virtual work in hyperelastic formulations.
 *
 * Evaluates the bilinear form arising from the linearization of the
 * internal virtual work:
 * @f[
 *   D(\delta W^{\text{int}})[\Delta\mathbf{u}, \mathbf{v}]
 *     = \int_{\Omega_0} D\mathbf{P}[\nabla_0 \Delta\mathbf{u}]
 *       : \nabla_0 \mathbf{v} \, dX
 * @f]
 * where @f$ D\mathbf{P}[\cdot] = \mathbf{A} : (\cdot) @f$ denotes the action of
 * the material tangent @f$ \mathbf{A} = \partial\mathbf{P}/\partial\mathbf{F} @f$.
 *
 * The integrator is generic: it obtains the finite element basis from the
 * FE space (not hardcoded to P1), supports arbitrary quadrature rules, and
 * builds a ConstitutivePoint (composed over Geometry::Point) at each
 * quadrature point for constitutive evaluation.
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORKTANGENT_H
#define RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORKTANGENT_H

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
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/TrialFunction.h"
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
   * @brief Local bilinear form integrator for the internal virtual work tangent
   * in hyperelastic problems.
   *
   * Computes the element-level tangent stiffness matrix:
   * @f[
   *   K^e_{ij} = \int_{K} D\mathbf{P}[\nabla_0 \phi_j] : \nabla_0 \phi_i \, dX
   * @f]
   *
   * The material tangent action @f$ \mathbf{A} : \delta\mathbf{F} @f$ encodes
   * both material and geometric stiffness (no separate geometric-stiffness
   * integrator is needed when using the first Piola-Kirchhoff formulation).
   *
   * Typically constructed via @c InternalVirtualWork::Tangent(u, v) rather
   * than directly.
   *
   * @tparam LawDerived The hyperelastic constitutive law type
   * @tparam TrialFunctionType The trial function type (backend-generic)
   * @tparam TestFunctionType The test function type (backend-generic)
   * @tparam DisplacementType The current displacement grid-function type
   */
  template <class LawDerived, class TrialFunctionType,
            class TestFunctionType, class DisplacementType>
  class InternalVirtualWorkTangent final
    : public Variational::LocalBilinearFormIntegratorBase<Real>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Parent class type.
      using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using TrialType = TrialFunctionType;
      using TestType = TestFunctionType;
      using StateType = DisplacementType;

      using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

      static_assert(Variational::IsTrialFunction<TrialType>::Value,
        "Solid::InternalVirtualWorkTangent expects a Rodin trial function.");
      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::InternalVirtualWorkTangent expects a Rodin test function.");

      /**
       * @brief Constructs the internal virtual work tangent integrator.
       * @param law The constitutive law (stored by value)
       * @param u The trial function
       * @param v The test function
       * @param displacement Current displacement grid function
       */
      InternalVirtualWorkTangent(
          const LawDerived& law,
          const TrialType& u,
          const TestType& v,
          const StateType& displacement)
        : Parent(u, v),
          m_law(law),
          m_trial(u),
          m_test(v),
          m_displacement(displacement),
          m_trialfes(u.getFiniteElementSpace()),
          m_testfes(v.getFiniteElementSpace()),
          m_statefes(displacement.getFiniteElementSpace()),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      InternalVirtualWorkTangent(const InternalVirtualWorkTangent& other)
        : Parent(other),
          m_law(other.m_law),
          m_trial(other.m_trial),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_trialfes(other.m_trialfes),
          m_testfes(other.m_testfes),
          m_statefes(other.m_statefes),
          m_quadOrder(other.m_quadOrder),
          m_input(other.m_input)
      {}

      /**
       * @brief Rebinds the linearization point.
       * @param gf Displacement grid function with the same concrete state type
       * @returns Reference to this object for chaining
       */
      InternalVirtualWorkTangent& setDisplacement(const StateType& gf)
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
      InternalVirtualWorkTangent& setQuadratureOrder(size_t order)
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
      InternalVirtualWorkTangent& setInput(InputFunction input)
      {
        m_input = std::move(input);
        return *this;
      }

      InternalVirtualWorkTangent& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trialFES = m_trialfes.get();
        const auto& testFES = m_testfes.get();
        const auto& stateFES = m_statefes.get();
        const size_t vdim = testFES.getVectorDimension();

        const auto& trialFE = trialFES.getFiniteElement(d, idx);
        const auto& testFE = testFES.getFiniteElement(d, idx);
        const auto& stateFE = stateFES.getFiniteElement(d, idx);
        const size_t trialDofs = trialFE.getCount();
        const size_t testDofs = testFE.getCount();

        // Determine effective quadrature order
        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * std::max({ trialFE.getOrder(),
                           testFE.getOrder(),
                           stateFE.getOrder() });
        const auto& qf = QF::PolytopeQuadratureFormula::get(effectiveOrder, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        // Zero element stiffness matrix
        m_matrix.resize(testDofs, trialDofs);
        m_matrix.setZero();

        auto stateGradient = Variational::Jacobian(m_displacement.get());
        auto trialGradient = Variational::Jacobian(m_trial.get());
        auto testGradient = Variational::Jacobian(m_test.get());

        // Loop over quadrature points
        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const ScalarType wq = qf.getWeight(q);

          const ScalarType distortion = pt.getDistortion();

          trialGradient.setIntegrationPoint(ip);
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

          // Build element stiffness at this quadrature point.
          // For trial DOF tr, the deformation gradient perturbation dF equals
          // the physical Jacobian of the trial basis function.
          for (size_t tr = 0; tr < trialDofs; ++tr)
          {
            const auto dF = trialGradient.getBasis(tr);

            // Compute material tangent action dP = A : dF
            Math::SpatialMatrix<ScalarType> dP;
            m_law.getMaterialTangent(dP, cache, cp, dF);

            // K_{te,tr} += wq * distortion * (dP : grad psi_te)
            for (size_t te = 0; te < testDofs; ++te)
            {
              const auto gradTest = testGradient.getBasis(te);
              ScalarType val = 0;
              for (size_t c = 0; c < vdim; ++c)
                for (size_t k = 0; k < d; ++k)
                  val += dP(c, k) * gradTest(c, k);
              m_matrix(te, tr) += wq * distortion * val;
            }
          }
        }

        return *this;
      }

      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
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

      InternalVirtualWorkTangent* copy() const noexcept final override
      {
        return new InternalVirtualWorkTangent(*this);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      void checkCompatibility(const StateType& displacement) const
      {
        const auto& trialFES = m_trialfes.get();
        const auto& testFES = m_testfes.get();
        const auto& stateFES = displacement.getFiniteElementSpace();

        assert(&trialFES.getMesh() == &testFES.getMesh());
        assert(&stateFES.getMesh() == &testFES.getMesh());
        assert(trialFES.getVectorDimension() == testFES.getVectorDimension());
        assert(stateFES.getVectorDimension() == testFES.getVectorDimension());
        (void) trialFES;
        (void) testFES;
        (void) stateFES;
      }

      LawType m_law;
      std::reference_wrapper<const TrialType> m_trial;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      std::reference_wrapper<const TrialFESType> m_trialfes;
      std::reference_wrapper<const TestFESType> m_testfes;
      std::reference_wrapper<const StateFESType> m_statefes;
      size_t m_quadOrder;
      InputFunction m_input;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Matrix<ScalarType> m_matrix;
  };

  /// CTAD deduction guide for InternalVirtualWorkTangent
  template <class LawDerived, class TrialFunctionType,
            class TestFunctionType, class DisplacementType>
  InternalVirtualWorkTangent(const LawDerived&,
                              const TrialFunctionType&,
                              const TestFunctionType&,
                              const DisplacementType&)
    -> InternalVirtualWorkTangent<
         LawDerived,
         std::decay_t<TrialFunctionType>,
         std::decay_t<TestFunctionType>,
         std::decay_t<DisplacementType>>;
}

#endif
