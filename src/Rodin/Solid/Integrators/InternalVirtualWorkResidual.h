/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalVirtualWorkResidual.h
 * @brief Residual integrators for the internal virtual work in hyperelastic
 * formulations (pure displacement and mixed displacement-pressure), written
 * in terms of the first Piola-Kirchhoff stress.
 *
 * Pure displacement:
 * @f[
 *   \delta W^{\text{int}}(\mathbf{v}) = \int_{\Omega_0} \mathbf{P}(\mathbf{u}) : \nabla_0 \mathbf{v} \, dX
 * @f]
 *
 * Mixed u-p formulation with @f$ \mathbf{P} = \mathbf{P}_{\text{iso}} +
 * p\,J\,\mathbf{F}^{-T} @f$:
 * - Momentum residual @f$ R_u @f$:
 *   InternalVirtualWorkResidual(law, v, displacement, pressure)
 * - Incompressibility residual @f$ R_p(q) = \int_{\Omega_0} (J - 1)\, q \, dX @f$:
 *   InternalVirtualWorkResidualP(q, displacement)
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
  template <class... Args>
  class InternalVirtualWorkResidual;

  /**
   * @brief Momentum residual (pure displacement formulation).
   *
   * @f[
   *   (\delta W^{\text{int}})^e_i = \int_{K} \mathbf{P} : \nabla_0 \phi_i \, dX
   * @f]
   */
  template <class LawDerived, class TestFunctionType, class DisplacementType>
  class InternalVirtualWorkResidual<LawDerived, TestFunctionType,
                                    DisplacementType> final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using TestType = TestFunctionType;
      using StateType = DisplacementType;

      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::InternalVirtualWorkResidual expects a Rodin test function.");

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

      InternalVirtualWorkResidual& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        m_statefes = std::cref(gf.getFiniteElementSpace());
        return *this;
      }

      InternalVirtualWorkResidual& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

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

        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * std::max(testFE.getOrder(), stateFE.getOrder());
        const auto& qf = QF::PolytopeQuadratureFormula::get(effectiveOrder, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        auto stateGradient = Variational::Jacobian(m_displacement.get());
        auto testGradient = Variational::Jacobian(m_test.get());

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const ScalarType wq = qf.getWeight(q);

          const ScalarType distortion = pt.getDistortion();

          testGradient.setIntegrationPoint(ip);
          const auto H = stateGradient.getValue(ip);

          KinematicState state(d);
          state.setDisplacementGradient(H);

          ConstitutivePoint cp(pt, state);
          cp.set<Tags::CellIndex>(idx);
          cp.set<Tags::QuadraturePointIndex>(q);

          if (m_input)
            m_input(cp);

          typename LawType::Cache cache;
          m_law.setCache(cache, cp);

          Math::SpatialMatrix<ScalarType> P;
          m_law.getFirstPiolaKirchhoffStress(P, cache, cp);

          // R_te += w P : grad v_te
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

  /**
   * @brief Momentum residual of the mixed u-p formulation. Bound to the
   * displacement test function @f$ \mathbf{v} @f$.
   *
   * @f[
   *   R^e_i = \int_{K} \big( \mathbf{P}_{\text{iso}}
   *     + p\,J\,\mathbf{F}^{-T} \big) : \nabla_0 \phi_i \, dX
   * @f]
   *
   * The constitutive law must supply only the isochoric part of the stress;
   * the volumetric response is carried by the pressure field.
   */
  template <class LawDerived, class TestFunctionType, class DisplacementType,
            class PressureType>
  class InternalVirtualWorkResidual<LawDerived, TestFunctionType,
                                    DisplacementType, PressureType> final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using LawType = LawDerived;
      using TestType = TestFunctionType;
      using StateType = DisplacementType;
      using PressureStateType = PressureType;

      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;
      using PressureFESType =
        typename FormLanguage::Traits<PressureStateType>::FESType;

      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::InternalVirtualWorkResidual expects a Rodin test function.");

      InternalVirtualWorkResidual(
          const LawDerived& law,
          const TestType& v,
          const StateType& displacement,
          const PressureStateType& pressure)
        : Parent(v),
          m_law(law),
          m_test(v),
          m_displacement(displacement),
          m_pressure(pressure),
          m_testfes(v.getFiniteElementSpace()),
          m_statefes(displacement.getFiniteElementSpace()),
          m_pressfes(pressure.getFiniteElementSpace()),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      InternalVirtualWorkResidual(const InternalVirtualWorkResidual& other)
        : Parent(other),
          m_law(other.m_law),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_pressure(other.m_pressure),
          m_testfes(other.m_testfes),
          m_statefes(other.m_statefes),
          m_pressfes(other.m_pressfes),
          m_quadOrder(other.m_quadOrder),
          m_input(other.m_input)
      {}

      InternalVirtualWorkResidual& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        m_statefes = std::cref(gf.getFiniteElementSpace());
        return *this;
      }

      InternalVirtualWorkResidual& setPressure(const PressureStateType& gf)
      {
        m_pressure = std::cref(gf);
        m_pressfes = std::cref(gf.getFiniteElementSpace());
        return *this;
      }

      InternalVirtualWorkResidual& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

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
        const auto& pressFES = m_pressfes.get();
        const size_t vdim = testFES.getVectorDimension();

        const auto& testFE = testFES.getFiniteElement(d, idx);
        const auto& stateFE = stateFES.getFiniteElement(d, idx);
        const auto& pressFE = pressFES.getFiniteElement(d, idx);
        const size_t testDofs = testFE.getCount();

        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * std::max({testFE.getOrder(), stateFE.getOrder(),
                          pressFE.getOrder()});
        const auto& qf = QF::PolytopeQuadratureFormula::get(effectiveOrder, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        auto stateGradient = Variational::Jacobian(m_displacement.get());
        auto testGradient = Variational::Jacobian(m_test.get());

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const ScalarType wq = qf.getWeight(q);

          const ScalarType distortion = pt.getDistortion();

          testGradient.setIntegrationPoint(ip);
          const auto H = stateGradient.getValue(ip);

          KinematicState state(d);
          state.setDisplacementGradient(H);

          ConstitutivePoint cp(pt, state);
          cp.set<Tags::CellIndex>(idx);
          cp.set<Tags::QuadraturePointIndex>(q);

          if (m_input)
            m_input(cp);

          typename LawType::Cache cache;
          m_law.setCache(cache, cp);

          Math::SpatialMatrix<ScalarType> P;
          m_law.getFirstPiolaKirchhoffStress(P, cache, cp);

          const ScalarType p = m_pressure.get().getValue(ip);
          const ScalarType J = state.getJacobian();
          const auto& FinvT = state.getDeformationGradientInverseTranspose();

          // R_te += w (P_iso + p J F^{-T}) : grad v_te
          for (size_t te = 0; te < testDofs; ++te)
          {
            const auto gradTest = testGradient.getBasis(te);
            ScalarType val = 0;
            for (size_t c = 0; c < vdim; ++c)
              for (size_t k = 0; k < d; ++k)
                val += (P(c, k) + p * J * FinvT(c, k)) * gradTest(c, k);
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

      const LawType& getLaw() const { return m_law; }

    private:
      void checkCompatibility(const StateType& displacement) const
      {
        const auto& testFES = m_testfes.get();
        const auto& stateFES = displacement.getFiniteElementSpace();
        const auto& pressFES = m_pressfes.get();

        assert(&stateFES.getMesh() == &testFES.getMesh());
        assert(&pressFES.getMesh() == &testFES.getMesh());
        assert(stateFES.getVectorDimension() == testFES.getVectorDimension());
        (void) testFES;
        (void) stateFES;
        (void) pressFES;
      }

      LawType m_law;
      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      std::reference_wrapper<const PressureStateType> m_pressure;
      std::reference_wrapper<const TestFESType> m_testfes;
      std::reference_wrapper<const StateFESType> m_statefes;
      std::reference_wrapper<const PressureFESType> m_pressfes;
      size_t m_quadOrder;
      InputFunction m_input;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /**
   * @brief Incompressibility constraint residual of the mixed u-p
   * formulation. Bound to the pressure test function @f$ q @f$.
   *
   * @f[
   *   R^e_i = \int_{K} (J - 1)\, \psi_i \, dX
   * @f]
   */
  template <class TestPressFunctionType, class DisplacementType>
  class InternalVirtualWorkResidualP final
    : public Variational::LinearFormIntegratorBase<Real>
  {
    public:
      using ScalarType = Real;
      using Parent = Variational::LinearFormIntegratorBase<ScalarType>;
      using TestType = TestPressFunctionType;
      using StateType = DisplacementType;

      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

      static_assert(Variational::IsTestFunction<TestType>::Value,
        "Solid::InternalVirtualWorkResidualP expects a Rodin test function.");

      InternalVirtualWorkResidualP(
          const TestType& q,
          const StateType& displacement)
        : Parent(q),
          m_test(q),
          m_displacement(displacement),
          m_testfes(q.getFiniteElementSpace()),
          m_statefes(displacement.getFiniteElementSpace()),
          m_quadOrder(0)
      {
        checkCompatibility(displacement);
      }

      InternalVirtualWorkResidualP(const InternalVirtualWorkResidualP& other)
        : Parent(other),
          m_test(other.m_test),
          m_displacement(other.m_displacement),
          m_testfes(other.m_testfes),
          m_statefes(other.m_statefes),
          m_quadOrder(other.m_quadOrder)
      {}

      InternalVirtualWorkResidualP& setDisplacement(const StateType& gf)
      {
        checkCompatibility(gf);
        m_displacement = std::cref(gf);
        m_statefes = std::cref(gf.getFiniteElementSpace());
        return *this;
      }

      InternalVirtualWorkResidualP& setQuadratureOrder(size_t order)
      {
        m_quadOrder = order;
        return *this;
      }

      InternalVirtualWorkResidualP& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;

        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& testFES = m_testfes.get();
        const auto& stateFES = m_statefes.get();

        const auto& testFE = testFES.getFiniteElement(d, idx);
        const auto& stateFE = stateFES.getFiniteElement(d, idx);
        const size_t testDofs = testFE.getCount();

        const size_t effectiveOrder = (m_quadOrder > 0)
          ? m_quadOrder
          : 2 * std::max(testFE.getOrder(), stateFE.getOrder());
        const auto& qf = QF::PolytopeQuadratureFormula::get(effectiveOrder, polytope.getGeometry());
        const auto& quadrature = polytope.getQuadrature(qf);
        const size_t nqp = quadrature.getSize();

        m_elemVec.resize(testDofs);
        m_elemVec.setZero();

        auto stateGradient = Variational::Jacobian(m_displacement.get());

        for (size_t q = 0; q < nqp; ++q)
        {
          const auto& pt = quadrature.getPoint(q);
          const Variational::IntegrationPoint ip(pt, &qf, q);
          const ScalarType wq = qf.getWeight(q);

          const ScalarType distortion = pt.getDistortion();

          const auto H = stateGradient.getValue(ip);

          KinematicState state(d);
          state.setDisplacementGradient(H);

          const ScalarType J = state.getJacobian();
          const auto& rc = pt.getReferenceCoordinates();

          // R_te += w (J - 1) q_te
          for (size_t te = 0; te < testDofs; ++te)
          {
            const ScalarType qv = testFE.getBasis(te)(rc);
            m_elemVec(te) += wq * distortion * (J - 1) * qv;
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

      InternalVirtualWorkResidualP* copy() const noexcept final override
      {
        return new InternalVirtualWorkResidualP(*this);
      }

    private:
      void checkCompatibility(const StateType& displacement) const
      {
        const auto& testFES = m_testfes.get();
        const auto& stateFES = displacement.getFiniteElementSpace();

        assert(&stateFES.getMesh() == &testFES.getMesh());
        (void) testFES;
        (void) stateFES;
      }

      std::reference_wrapper<const TestType> m_test;
      std::reference_wrapper<const StateType> m_displacement;
      std::reference_wrapper<const TestFESType> m_testfes;
      std::reference_wrapper<const StateFESType> m_statefes;
      size_t m_quadOrder;

      Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
      Math::Vector<ScalarType> m_elemVec;
  };

  /// CTAD deduction guide for the displacement-only residual
  template <class LawDerived, class TestFunctionType, class DisplacementType>
  InternalVirtualWorkResidual(const LawDerived&,
                               const TestFunctionType&,
                               const DisplacementType&)
    -> InternalVirtualWorkResidual<
         LawDerived,
         std::decay_t<TestFunctionType>,
         std::decay_t<DisplacementType>>;

  /// CTAD deduction guide for the mixed u-p momentum residual
  template <class LawDerived, class TestFunctionType, class DisplacementType,
            class PressureType>
  InternalVirtualWorkResidual(const LawDerived&,
                               const TestFunctionType&,
                               const DisplacementType&,
                               const PressureType&)
    -> InternalVirtualWorkResidual<
         LawDerived,
         std::decay_t<TestFunctionType>,
         std::decay_t<DisplacementType>,
         std::decay_t<PressureType>>;

  /// CTAD deduction guide for the incompressibility constraint residual
  template <class TestPressFunctionType, class DisplacementType>
  InternalVirtualWorkResidualP(const TestPressFunctionType&,
                               const DisplacementType&)
    -> InternalVirtualWorkResidualP<
         std::decay_t<TestPressFunctionType>,
         std::decay_t<DisplacementType>>;
}

#endif
