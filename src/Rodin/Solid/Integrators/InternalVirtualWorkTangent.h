/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalVirtualWorkTangent.h
 * @brief Tangent integrators for the internal virtual work in hyperelastic
 * formulations (pure displacement and mixed displacement-pressure), written
 * in terms of the first Piola-Kirchhoff stress.
 *
 * Pure displacement block:
 * @f[
 *   D(\delta W^{\text{int}})[\Delta\mathbf{u}, \mathbf{v}]
 *     = \int_{\Omega_0} D\mathbf{P}[\nabla_0 \Delta\mathbf{u}]
 *       : \nabla_0 \mathbf{v} \, dX
 * @f]
 *
 * Mixed u-p formulation with @f$ \mathbf{P} = \mathbf{P}_{\text{iso}} +
 * p\,J\,\mathbf{F}^{-T} @f$. Each block of the Newton matrix
 * @f[
 *   \begin{pmatrix} K_{uu} & K_{up} \\ K_{pu} & 0 \end{pmatrix}
 * @f]
 * is a separate integrator bound to its own (trial, test) pair, so that the
 * blockwise Problem assembly routes it to the correct rows and columns:
 * - @f$ K_{uu} @f$: InternalVirtualWorkTangent(law, u, v, displacement, pressure)
 * - @f$ K_{up} @f$: InternalVirtualWorkTangentUP(p, v, displacement)
 * - @f$ K_{pu} @f$: InternalVirtualWorkTangentPU(u, q, displacement)
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORKTANGENT_H
#define RODIN_SOLID_INTEGRATORS_INTERNALVIRTUALWORKTANGENT_H

#include <algorithm>
#include <cassert>
#include <functional>
#include <type_traits>
#include <vector>

#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Types.h"
#include "Rodin/Variational/BilinearFormIntegrator.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "Rodin/Variational/Jacobian.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/TrialFunction.h"

#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"
#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"
#include "Rodin/Solid/Local/Input.h"

<<<<<<< HEAD
namespace Rodin::Solid {

template <class... Args>
class InternalVirtualWorkTangent;
=======
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
      /// @brief Constitutive law type.
      using LawType = LawDerived;
      /// @brief Trial function type.
      using TrialType = TrialFunctionType;
      /// @brief Test function type.
      using TestType = TestFunctionType;
      /// @brief Current displacement state type.
      using StateType = DisplacementType;

      /// @brief Trial finite element space type.
      using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
      /// @brief Test finite element space type.
      using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
      /// @brief State finite element space type.
      using StateFESType = typename FormLanguage::Traits<StateType>::FESType;
>>>>>>> origin/develop

/**
 * @brief Displacement tangent (pure displacement formulation).
 *
 * @f[
 *   K^e_{ij} = \int_{K} D\mathbf{P}[\nabla_0 \phi_j] : \nabla_0 \phi_i \, dX
 * @f]
 */
template <class LawDerived, class TrialFunctionType, class TestFunctionType,
          class DisplacementType>
class InternalVirtualWorkTangent<LawDerived, TrialFunctionType,
                                 TestFunctionType, DisplacementType> final
    : public Variational::LocalBilinearFormIntegratorBase<Real> {
public:
  using ScalarType = Real;
  using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
  using LawType = LawDerived;
  using TrialType = TrialFunctionType;
  using TestType = TestFunctionType;
  using StateType = DisplacementType;

  using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
  using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
  using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

  static_assert(
      Variational::IsTrialFunction<TrialType>::Value,
      "Solid::InternalVirtualWorkTangent expects a Rodin trial function.");
  static_assert(
      Variational::IsTestFunction<TestType>::Value,
      "Solid::InternalVirtualWorkTangent expects a Rodin test function.");

  InternalVirtualWorkTangent(const LawDerived &law, const TrialType &u,
                             const TestType &v, const StateType &displacement)
      : Parent(u, v), m_law(law), m_trial(u), m_test(v),
        m_displacement(displacement), m_trialfes(u.getFiniteElementSpace()),
        m_testfes(v.getFiniteElementSpace()),
        m_statefes(displacement.getFiniteElementSpace()), m_quadOrder(0) {
    checkCompatibility(displacement);
  }

  InternalVirtualWorkTangent(const InternalVirtualWorkTangent &other)
      : Parent(other), m_law(other.m_law), m_trial(other.m_trial),
        m_test(other.m_test), m_displacement(other.m_displacement),
        m_trialfes(other.m_trialfes), m_testfes(other.m_testfes),
        m_statefes(other.m_statefes), m_quadOrder(other.m_quadOrder),
        m_input(other.m_input) {}

  InternalVirtualWorkTangent &setDisplacement(const StateType &gf) {
    checkCompatibility(gf);
    m_displacement = std::cref(gf);
    m_statefes = std::cref(gf.getFiniteElementSpace());
    return *this;
  }

  InternalVirtualWorkTangent &setQuadratureOrder(size_t order) {
    m_quadOrder = order;
    return *this;
  }

  InternalVirtualWorkTangent &setInput(InputFunction input) {
    m_input = std::move(input);
    return *this;
  }

  InternalVirtualWorkTangent &
  setPolytope(const Geometry::Polytope &polytope) final override {
    m_polytope = polytope;

    const size_t d = polytope.getDimension();
    const Index idx = polytope.getIndex();
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = m_statefes.get();
    const size_t vdim = testFES.getVectorDimension();

    const auto &trialFE = trialFES.getFiniteElement(d, idx);
    const auto &testFE = testFES.getFiniteElement(d, idx);
    const auto &stateFE = stateFES.getFiniteElement(d, idx);
    const size_t trialDofs = trialFE.getCount();
    const size_t testDofs = testFE.getCount();

    const size_t effectiveOrder =
        (m_quadOrder > 0) ? m_quadOrder
                          : 2 * std::max({trialFE.getOrder(), testFE.getOrder(),
                                          stateFE.getOrder()});
    const auto &qf = QF::PolytopeQuadratureFormula::get(effectiveOrder,
                                                        polytope.getGeometry());
    const auto &quadrature = polytope.getQuadrature(qf);
    const size_t nqp = quadrature.getSize();

    m_matrix.resize(testDofs, trialDofs);
    m_matrix.setZero();

    auto stateGradient = Variational::Jacobian(m_displacement.get());
    auto trialGradient = Variational::Jacobian(m_trial.get());
    auto testGradient = Variational::Jacobian(m_test.get());

    for (size_t q = 0; q < nqp; ++q) {
      const auto &pt = quadrature.getPoint(q);
      const Variational::IntegrationPoint ip(pt, &qf, q);
      const ScalarType wq = qf.getWeight(q);

      const ScalarType distortion = pt.getDistortion();

      trialGradient.setIntegrationPoint(ip);
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

      for (size_t tr = 0; tr < trialDofs; ++tr) {
        const auto dF = trialGradient.getBasis(tr);

        Math::SpatialMatrix<ScalarType> dP;
        m_law.getMaterialTangent(dP, cache, cp, dF);

        // K_{te,tr} += w (A : grad du_tr) : grad v_te
        for (size_t te = 0; te < testDofs; ++te) {
          const auto gradTest = testGradient.getBasis(te);
          ScalarType val = 0;
          for (size_t c = 0; c < vdim; ++c)
            for (size_t k = 0; k < d; ++k)
              val += dP(c, k) * gradTest(c, k);
          m_matrix(te, tr) += wq * distortion * val;
        }
      }
    }

<<<<<<< HEAD
    return *this;
  }
=======
      /// @brief Copy constructor.
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
>>>>>>> origin/develop

  ScalarType integrate(size_t tr, size_t te) final override {
    return m_matrix(te, tr);
  }

  const Geometry::Polytope &getPolytope() const final override {
    assert(m_polytope);
    return m_polytope->get();
  }

  Geometry::Region getRegion() const final override {
    return Geometry::Region::Cells;
  }

<<<<<<< HEAD
  InternalVirtualWorkTangent *copy() const noexcept final override {
    return new InternalVirtualWorkTangent(*this);
  }
=======
      /// @brief Sets the current polytope and assembles the element tangent matrix.
      InternalVirtualWorkTangent& setPolytope(const Geometry::Polytope& polytope) final override
      {
        m_polytope = polytope;
>>>>>>> origin/develop

  const LawType &getLaw() const { return m_law; }

private:
  void checkCompatibility(const StateType &displacement) const {
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = displacement.getFiniteElementSpace();

    assert(&trialFES.getMesh() == &testFES.getMesh());
    assert(&stateFES.getMesh() == &testFES.getMesh());
    assert(trialFES.getVectorDimension() == testFES.getVectorDimension());
    assert(stateFES.getVectorDimension() == testFES.getVectorDimension());
    (void)trialFES;
    (void)testFES;
    (void)stateFES;
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

/**
 * @brief Displacement-displacement block @f$ K_{uu} @f$ of the mixed u-p
 * formulation. Bound to the (u, v) pair.
 *
 * @f[
 *   K^e_{ij} = \int_K \Big( D\mathbf{P}_{\text{iso}}[\nabla_0\phi_j]
 *     + p\,J\big[(\mathbf{F}^{-T} : \nabla_0\phi_j)\,\mathbf{F}^{-T}
 *     - \mathbf{F}^{-T} (\nabla_0\phi_j)^T \mathbf{F}^{-T}\big]\Big)
 *     : \nabla_0\phi_i \, dX
 * @f]
 *
 * The constitutive law must supply only the isochoric part of the stress;
 * the volumetric response is carried by the pressure field.
 */
template <class LawDerived, class TrialFunctionType, class TestFunctionType,
          class DisplacementType, class PressureType>
class InternalVirtualWorkTangent<LawDerived, TrialFunctionType,
                                 TestFunctionType, DisplacementType,
                                 PressureType> final
    : public Variational::LocalBilinearFormIntegratorBase<Real> {
public:
  using ScalarType = Real;
  using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
  using LawType = LawDerived;
  using TrialType = TrialFunctionType;
  using TestType = TestFunctionType;
  using StateType = DisplacementType;
  using PressureStateType = PressureType;

  using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
  using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
  using StateFESType = typename FormLanguage::Traits<StateType>::FESType;
  using PressureFESType =
      typename FormLanguage::Traits<PressureStateType>::FESType;

  static_assert(
      Variational::IsTrialFunction<TrialType>::Value,
      "Solid::InternalVirtualWorkTangent expects a Rodin trial function.");
  static_assert(
      Variational::IsTestFunction<TestType>::Value,
      "Solid::InternalVirtualWorkTangent expects a Rodin test function.");

  InternalVirtualWorkTangent(const LawDerived &law, const TrialType &u,
                             const TestType &v, const StateType &displacement,
                             const PressureStateType &pressure)
      : Parent(u, v), m_law(law), m_trial(u), m_test(v),
        m_displacement(displacement), m_pressure(pressure),
        m_trialfes(u.getFiniteElementSpace()),
        m_testfes(v.getFiniteElementSpace()),
        m_statefes(displacement.getFiniteElementSpace()),
        m_pressfes(pressure.getFiniteElementSpace()), m_quadOrder(0) {
    checkCompatibility(displacement);
  }

  InternalVirtualWorkTangent(const InternalVirtualWorkTangent &other)
      : Parent(other), m_law(other.m_law), m_trial(other.m_trial),
        m_test(other.m_test), m_displacement(other.m_displacement),
        m_pressure(other.m_pressure), m_trialfes(other.m_trialfes),
        m_testfes(other.m_testfes), m_statefes(other.m_statefes),
        m_pressfes(other.m_pressfes), m_quadOrder(other.m_quadOrder),
        m_input(other.m_input) {}

  InternalVirtualWorkTangent &setDisplacement(const StateType &gf) {
    checkCompatibility(gf);
    m_displacement = std::cref(gf);
    m_statefes = std::cref(gf.getFiniteElementSpace());
    return *this;
  }

  InternalVirtualWorkTangent &setPressure(const PressureStateType &gf) {
    m_pressure = std::cref(gf);
    m_pressfes = std::cref(gf.getFiniteElementSpace());
    return *this;
  }

  InternalVirtualWorkTangent &setQuadratureOrder(size_t order) {
    m_quadOrder = order;
    return *this;
  }

  InternalVirtualWorkTangent &setInput(InputFunction input) {
    m_input = std::move(input);
    return *this;
  }

  InternalVirtualWorkTangent &
  setPolytope(const Geometry::Polytope &polytope) final override {
    m_polytope = polytope;

    const size_t d = polytope.getDimension();
    const Index idx = polytope.getIndex();
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = m_statefes.get();
    const auto &pressFES = m_pressfes.get();
    const size_t vdim = testFES.getVectorDimension();

    const auto &trialFE = trialFES.getFiniteElement(d, idx);
    const auto &testFE = testFES.getFiniteElement(d, idx);
    const auto &stateFE = stateFES.getFiniteElement(d, idx);
    const auto &pressFE = pressFES.getFiniteElement(d, idx);
    const size_t trialDofs = trialFE.getCount();
    const size_t testDofs = testFE.getCount();

    const size_t effectiveOrder =
        (m_quadOrder > 0) ? m_quadOrder
                          : 2 * std::max({trialFE.getOrder(), testFE.getOrder(),
                                          stateFE.getOrder(),
                                          pressFE.getOrder()});
    const auto &qf = QF::PolytopeQuadratureFormula::get(effectiveOrder,
                                                        polytope.getGeometry());
    const auto &quadrature = polytope.getQuadrature(qf);
    const size_t nqp = quadrature.getSize();

    m_matrix.resize(testDofs, trialDofs);
    m_matrix.setZero();

    auto stateGradient = Variational::Jacobian(m_displacement.get());
    auto trialGradient = Variational::Jacobian(m_trial.get());
    auto testGradient = Variational::Jacobian(m_test.get());

    for (size_t q = 0; q < nqp; ++q) {
      const auto &pt = quadrature.getPoint(q);
      const Variational::IntegrationPoint ip(pt, &qf, q);
      const ScalarType wq = qf.getWeight(q);

      const ScalarType distortion = pt.getDistortion();

      trialGradient.setIntegrationPoint(ip);
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

      const ScalarType p = m_pressure.get().getValue(ip);
      const ScalarType J = state.getJacobian();
      const auto &FinvT = state.getDeformationGradientInverseTranspose();

      for (size_t tr = 0; tr < trialDofs; ++tr) {
        const auto dF = trialGradient.getBasis(tr);

        Math::SpatialMatrix<ScalarType> dP;
        m_law.getMaterialTangent(dP, cache, cp, dF);

        // DJ[du] = J (F^{-T} : grad du)
        ScalarType FinvTdF = 0;
        for (size_t c = 0; c < vdim; ++c)
          for (size_t k = 0; k < d; ++k)
            FinvTdF += FinvT(c, k) * dF(c, k);

        // D(F^{-T})[du] = -F^{-T} (grad du)^T F^{-T}
        const Math::SpatialMatrix<ScalarType> M =
            FinvT * dF.transpose() * FinvT;

        // K_{te,tr} += w [ (A : grad du) : grad v
        //   + p J ( (F^{-T} : grad du)(F^{-T} : grad v)
        //           - (F^{-T} (grad du)^T F^{-T}) : grad v ) ]
        for (size_t te = 0; te < testDofs; ++te) {
          const auto gradTest = testGradient.getBasis(te);
          ScalarType val = 0;
          ScalarType FinvTG = 0;
          ScalarType MG = 0;
          for (size_t c = 0; c < vdim; ++c) {
            for (size_t k = 0; k < d; ++k) {
              val += dP(c, k) * gradTest(c, k);
              FinvTG += FinvT(c, k) * gradTest(c, k);
              MG += M(c, k) * gradTest(c, k);
            }
          }
          val += p * J * (FinvTdF * FinvTG - MG);
          m_matrix(te, tr) += wq * distortion * val;
        }
      }
    }

<<<<<<< HEAD
    return *this;
  }

  ScalarType integrate(size_t tr, size_t te) final override {
    return m_matrix(te, tr);
  }

  const Geometry::Polytope &getPolytope() const final override {
    assert(m_polytope);
    return m_polytope->get();
  }

  Geometry::Region getRegion() const final override {
    return Geometry::Region::Cells;
  }

  InternalVirtualWorkTangent *copy() const noexcept final override {
    return new InternalVirtualWorkTangent(*this);
  }

  const LawType &getLaw() const { return m_law; }

private:
  void checkCompatibility(const StateType &displacement) const {
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = displacement.getFiniteElementSpace();
    const auto &pressFES = m_pressfes.get();

    assert(&trialFES.getMesh() == &testFES.getMesh());
    assert(&stateFES.getMesh() == &testFES.getMesh());
    assert(&pressFES.getMesh() == &testFES.getMesh());
    assert(trialFES.getVectorDimension() == testFES.getVectorDimension());
    assert(stateFES.getVectorDimension() == testFES.getVectorDimension());
    (void)trialFES;
    (void)testFES;
    (void)stateFES;
    (void)pressFES;
  }

  LawType m_law;
  std::reference_wrapper<const TrialType> m_trial;
  std::reference_wrapper<const TestType> m_test;
  std::reference_wrapper<const StateType> m_displacement;
  std::reference_wrapper<const PressureStateType> m_pressure;
  std::reference_wrapper<const TrialFESType> m_trialfes;
  std::reference_wrapper<const TestFESType> m_testfes;
  std::reference_wrapper<const StateFESType> m_statefes;
  std::reference_wrapper<const PressureFESType> m_pressfes;
  size_t m_quadOrder;
  InputFunction m_input;

  Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
  Math::Matrix<ScalarType> m_matrix;
};

/**
 * @brief Displacement-pressure coupling block @f$ K_{up} @f$. Bound to the
 * (p, v) pair: columns are pressure DOFs, rows are displacement test DOFs.
 *
 * @f[
 *   K^e_{ij} = \int_K \psi_j \, J\,\mathbf{F}^{-T} : \nabla_0\phi_i \, dX
 * @f]
 */
template <class TrialPressFunctionType, class TestFunctionType,
          class DisplacementType>
class InternalVirtualWorkTangentUP final
    : public Variational::LocalBilinearFormIntegratorBase<Real> {
public:
  using ScalarType = Real;
  using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
  using TrialType = TrialPressFunctionType;
  using TestType = TestFunctionType;
  using StateType = DisplacementType;

  using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
  using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
  using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

  static_assert(
      Variational::IsTrialFunction<TrialType>::Value,
      "Solid::InternalVirtualWorkTangentUP expects a Rodin trial function.");
  static_assert(
      Variational::IsTestFunction<TestType>::Value,
      "Solid::InternalVirtualWorkTangentUP expects a Rodin test function.");

  InternalVirtualWorkTangentUP(const TrialType &p, const TestType &v,
                               const StateType &displacement)
      : Parent(p, v), m_trial(p), m_test(v), m_displacement(displacement),
        m_trialfes(p.getFiniteElementSpace()),
        m_testfes(v.getFiniteElementSpace()),
        m_statefes(displacement.getFiniteElementSpace()), m_quadOrder(0) {
    checkCompatibility(displacement);
  }

  InternalVirtualWorkTangentUP(const InternalVirtualWorkTangentUP &other)
      : Parent(other), m_trial(other.m_trial), m_test(other.m_test),
        m_displacement(other.m_displacement), m_trialfes(other.m_trialfes),
        m_testfes(other.m_testfes), m_statefes(other.m_statefes),
        m_quadOrder(other.m_quadOrder) {}

  InternalVirtualWorkTangentUP &setDisplacement(const StateType &gf) {
    checkCompatibility(gf);
    m_displacement = std::cref(gf);
    m_statefes = std::cref(gf.getFiniteElementSpace());
    return *this;
  }

  InternalVirtualWorkTangentUP &setQuadratureOrder(size_t order) {
    m_quadOrder = order;
    return *this;
  }

  InternalVirtualWorkTangentUP &
  setPolytope(const Geometry::Polytope &polytope) final override {
    m_polytope = polytope;

    const size_t d = polytope.getDimension();
    const Index idx = polytope.getIndex();
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = m_statefes.get();
    const size_t vdim = testFES.getVectorDimension();

    const auto &trialFE = trialFES.getFiniteElement(d, idx);
    const auto &testFE = testFES.getFiniteElement(d, idx);
    const auto &stateFE = stateFES.getFiniteElement(d, idx);
    const size_t trialDofs = trialFE.getCount();
    const size_t testDofs = testFE.getCount();

    const size_t effectiveOrder =
        (m_quadOrder > 0) ? m_quadOrder
                          : 2 * std::max({trialFE.getOrder(), testFE.getOrder(),
                                          stateFE.getOrder()});
    const auto &qf = QF::PolytopeQuadratureFormula::get(effectiveOrder,
                                                        polytope.getGeometry());
    const auto &quadrature = polytope.getQuadrature(qf);
    const size_t nqp = quadrature.getSize();

    m_matrix.resize(testDofs, trialDofs);
    m_matrix.setZero();

    auto stateGradient = Variational::Jacobian(m_displacement.get());
    auto testGradient = Variational::Jacobian(m_test.get());

    for (size_t q = 0; q < nqp; ++q) {
      const auto &pt = quadrature.getPoint(q);
      const Variational::IntegrationPoint ip(pt, &qf, q);
      const ScalarType wq = qf.getWeight(q);

      const ScalarType distortion = pt.getDistortion();

      testGradient.setIntegrationPoint(ip);
      const auto H = stateGradient.getValue(ip);

      KinematicState state(d);
      state.setDisplacementGradient(H);

      const ScalarType J = state.getJacobian();
      const auto &FinvT = state.getDeformationGradientInverseTranspose();
      const auto &rc = pt.getReferenceCoordinates();

      // K_{te,tr} += w dp_tr J (F^{-T} : grad v_te)
      for (size_t te = 0; te < testDofs; ++te) {
        const auto gradTest = testGradient.getBasis(te);
        ScalarType FinvTG = 0;
        for (size_t c = 0; c < vdim; ++c)
          for (size_t k = 0; k < d; ++k)
            FinvTG += FinvT(c, k) * gradTest(c, k);
        for (size_t tr = 0; tr < trialDofs; ++tr) {
          const ScalarType dp = trialFE.getBasis(tr)(rc);
          m_matrix(te, tr) += wq * distortion * dp * J * FinvTG;
        }
=======
      /// @brief Returns an entry of the current element tangent matrix.
      ScalarType integrate(size_t tr, size_t te) final override
      {
        return m_matrix(te, tr);
>>>>>>> origin/develop
      }
    }

<<<<<<< HEAD
    return *this;
  }

  ScalarType integrate(size_t tr, size_t te) final override {
    return m_matrix(te, tr);
  }

  const Geometry::Polytope &getPolytope() const final override {
    assert(m_polytope);
    return m_polytope->get();
  }

  Geometry::Region getRegion() const final override {
    return Geometry::Region::Cells;
  }

  InternalVirtualWorkTangentUP *copy() const noexcept final override {
    return new InternalVirtualWorkTangentUP(*this);
  }

private:
  void checkCompatibility(const StateType &displacement) const {
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = displacement.getFiniteElementSpace();

    assert(&trialFES.getMesh() == &testFES.getMesh());
    assert(&stateFES.getMesh() == &testFES.getMesh());
    assert(stateFES.getVectorDimension() == testFES.getVectorDimension());
    (void)trialFES;
    (void)testFES;
    (void)stateFES;
  }

  std::reference_wrapper<const TrialType> m_trial;
  std::reference_wrapper<const TestType> m_test;
  std::reference_wrapper<const StateType> m_displacement;
  std::reference_wrapper<const TrialFESType> m_trialfes;
  std::reference_wrapper<const TestFESType> m_testfes;
  std::reference_wrapper<const StateFESType> m_statefes;
  size_t m_quadOrder;

  Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
  Math::Matrix<ScalarType> m_matrix;
};

/**
 * @brief Pressure-displacement coupling block @f$ K_{pu} @f$ (linearized
 * incompressibility constraint). Bound to the (u, q) pair: columns are
 * displacement DOFs, rows are pressure test DOFs.
 *
 * @f[
 *   K^e_{ij} = \int_K \psi_i \, J\,\mathbf{F}^{-T} : \nabla_0\phi_j \, dX
 * @f]
 */
template <class TrialFunctionType, class TestPressFunctionType,
          class DisplacementType>
class InternalVirtualWorkTangentPU final
    : public Variational::LocalBilinearFormIntegratorBase<Real> {
public:
  using ScalarType = Real;
  using Parent = Variational::LocalBilinearFormIntegratorBase<ScalarType>;
  using TrialType = TrialFunctionType;
  using TestType = TestPressFunctionType;
  using StateType = DisplacementType;

  using TrialFESType = typename FormLanguage::Traits<TrialType>::FESType;
  using TestFESType = typename FormLanguage::Traits<TestType>::FESType;
  using StateFESType = typename FormLanguage::Traits<StateType>::FESType;

  static_assert(
      Variational::IsTrialFunction<TrialType>::Value,
      "Solid::InternalVirtualWorkTangentPU expects a Rodin trial function.");
  static_assert(
      Variational::IsTestFunction<TestType>::Value,
      "Solid::InternalVirtualWorkTangentPU expects a Rodin test function.");

  InternalVirtualWorkTangentPU(const TrialType &u, const TestType &q,
                               const StateType &displacement)
      : Parent(u, q), m_trial(u), m_test(q), m_displacement(displacement),
        m_trialfes(u.getFiniteElementSpace()),
        m_testfes(q.getFiniteElementSpace()),
        m_statefes(displacement.getFiniteElementSpace()), m_quadOrder(0) {
    checkCompatibility(displacement);
  }

  InternalVirtualWorkTangentPU(const InternalVirtualWorkTangentPU &other)
      : Parent(other), m_trial(other.m_trial), m_test(other.m_test),
        m_displacement(other.m_displacement), m_trialfes(other.m_trialfes),
        m_testfes(other.m_testfes), m_statefes(other.m_statefes),
        m_quadOrder(other.m_quadOrder) {}

  InternalVirtualWorkTangentPU &setDisplacement(const StateType &gf) {
    checkCompatibility(gf);
    m_displacement = std::cref(gf);
    m_statefes = std::cref(gf.getFiniteElementSpace());
    return *this;
  }

  InternalVirtualWorkTangentPU &setQuadratureOrder(size_t order) {
    m_quadOrder = order;
    return *this;
  }

  InternalVirtualWorkTangentPU &
  setPolytope(const Geometry::Polytope &polytope) final override {
    m_polytope = polytope;

    const size_t d = polytope.getDimension();
    const Index idx = polytope.getIndex();
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = m_statefes.get();
    const size_t vdim = trialFES.getVectorDimension();

    const auto &trialFE = trialFES.getFiniteElement(d, idx);
    const auto &testFE = testFES.getFiniteElement(d, idx);
    const auto &stateFE = stateFES.getFiniteElement(d, idx);
    const size_t trialDofs = trialFE.getCount();
    const size_t testDofs = testFE.getCount();

    const size_t effectiveOrder =
        (m_quadOrder > 0) ? m_quadOrder
                          : 2 * std::max({trialFE.getOrder(), testFE.getOrder(),
                                          stateFE.getOrder()});
    const auto &qf = QF::PolytopeQuadratureFormula::get(effectiveOrder,
                                                        polytope.getGeometry());
    const auto &quadrature = polytope.getQuadrature(qf);
    const size_t nqp = quadrature.getSize();

    m_matrix.resize(testDofs, trialDofs);
    m_matrix.setZero();

    auto stateGradient = Variational::Jacobian(m_displacement.get());
    auto trialGradient = Variational::Jacobian(m_trial.get());

    for (size_t q = 0; q < nqp; ++q) {
      const auto &pt = quadrature.getPoint(q);
      const Variational::IntegrationPoint ip(pt, &qf, q);
      const ScalarType wq = qf.getWeight(q);

      const ScalarType distortion = pt.getDistortion();

      trialGradient.setIntegrationPoint(ip);
      const auto H = stateGradient.getValue(ip);

      KinematicState state(d);
      state.setDisplacementGradient(H);

      const ScalarType J = state.getJacobian();
      const auto &FinvT = state.getDeformationGradientInverseTranspose();
      const auto &rc = pt.getReferenceCoordinates();

      // K_{te,tr} += w q_te DJ[du_tr],  DJ[du] = J (F^{-T} : grad du)
      for (size_t tr = 0; tr < trialDofs; ++tr) {
        const auto dF = trialGradient.getBasis(tr);
        ScalarType FinvTdF = 0;
        for (size_t c = 0; c < vdim; ++c)
          for (size_t k = 0; k < d; ++k)
            FinvTdF += FinvT(c, k) * dF(c, k);
        const ScalarType dJ = J * FinvTdF;
        for (size_t te = 0; te < testDofs; ++te) {
          const ScalarType qv = testFE.getBasis(te)(rc);
          m_matrix(te, tr) += wq * distortion * qv * dJ;
        }
=======
      /// @brief Returns the current polytope.
      const Geometry::Polytope& getPolytope() const final override
      {
        assert(m_polytope);
        return m_polytope->get();
>>>>>>> origin/develop
      }
    }

<<<<<<< HEAD
    return *this;
  }

  ScalarType integrate(size_t tr, size_t te) final override {
    return m_matrix(te, tr);
  }
=======
      /// @brief Returns the integration region.
      Geometry::Region getRegion() const final override
      {
        return Geometry::Region::Cells;
      }

      /// @brief Polymorphically copies this tangent integrator.
      InternalVirtualWorkTangent* copy() const noexcept final override
      {
        return new InternalVirtualWorkTangent(*this);
      }
>>>>>>> origin/develop

  const Geometry::Polytope &getPolytope() const final override {
    assert(m_polytope);
    return m_polytope->get();
  }

  Geometry::Region getRegion() const final override {
    return Geometry::Region::Cells;
  }

  InternalVirtualWorkTangentPU *copy() const noexcept final override {
    return new InternalVirtualWorkTangentPU(*this);
  }

private:
  void checkCompatibility(const StateType &displacement) const {
    const auto &trialFES = m_trialfes.get();
    const auto &testFES = m_testfes.get();
    const auto &stateFES = displacement.getFiniteElementSpace();

    assert(&trialFES.getMesh() == &testFES.getMesh());
    assert(&stateFES.getMesh() == &testFES.getMesh());
    assert(stateFES.getVectorDimension() == trialFES.getVectorDimension());
    (void)trialFES;
    (void)testFES;
    (void)stateFES;
  }

  std::reference_wrapper<const TrialType> m_trial;
  std::reference_wrapper<const TestType> m_test;
  std::reference_wrapper<const StateType> m_displacement;
  std::reference_wrapper<const TrialFESType> m_trialfes;
  std::reference_wrapper<const TestFESType> m_testfes;
  std::reference_wrapper<const StateFESType> m_statefes;
  size_t m_quadOrder;

  Optional<std::reference_wrapper<const Geometry::Polytope>> m_polytope;
  Math::Matrix<ScalarType> m_matrix;
};

/// CTAD deduction guide for the displacement-only tangent
template <class LawDerived, class TrialFunctionType, class TestFunctionType,
          class DisplacementType>
InternalVirtualWorkTangent(const LawDerived &, const TrialFunctionType &,
                           const TestFunctionType &, const DisplacementType &)
    -> InternalVirtualWorkTangent<LawDerived, std::decay_t<TrialFunctionType>,
                                  std::decay_t<TestFunctionType>,
                                  std::decay_t<DisplacementType>>;

/// CTAD deduction guide for the mixed u-p displacement block
template <class LawDerived, class TrialFunctionType, class TestFunctionType,
          class DisplacementType, class PressureType>
InternalVirtualWorkTangent(const LawDerived &, const TrialFunctionType &,
                           const TestFunctionType &, const DisplacementType &,
                           const PressureType &)
    -> InternalVirtualWorkTangent<LawDerived, std::decay_t<TrialFunctionType>,
                                  std::decay_t<TestFunctionType>,
                                  std::decay_t<DisplacementType>,
                                  std::decay_t<PressureType>>;

/// CTAD deduction guide for the K_up block
template <class TrialPressFunctionType, class TestFunctionType,
          class DisplacementType>
InternalVirtualWorkTangentUP(const TrialPressFunctionType &,
                             const TestFunctionType &,
                             const DisplacementType &)
    -> InternalVirtualWorkTangentUP<std::decay_t<TrialPressFunctionType>,
                                    std::decay_t<TestFunctionType>,
                                    std::decay_t<DisplacementType>>;

/// CTAD deduction guide for the K_pu block
template <class TrialFunctionType, class TestPressFunctionType,
          class DisplacementType>
InternalVirtualWorkTangentPU(const TrialFunctionType &,
                             const TestPressFunctionType &,
                             const DisplacementType &)
    -> InternalVirtualWorkTangentPU<std::decay_t<TrialFunctionType>,
                                    std::decay_t<TestPressFunctionType>,
                                    std::decay_t<DisplacementType>>;
} // namespace Rodin::Solid

#endif
