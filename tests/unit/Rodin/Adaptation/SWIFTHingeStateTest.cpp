/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Assembly.h"
#include "Rodin/Adaptation/SWIFT/Hinge.h"
#include "Rodin/Adaptation/SWIFT/HingeForce.h"
#include "Rodin/Adaptation/SWIFT/HingeMetric.h"

using namespace Rodin;
using namespace Rodin::Adaptation;

namespace Rodin::Tests::Unit
{

  TEST(Rodin_Adaptation_SWIFTHingeState, AffineHingeDerivatives2D3D)
  {
    for (const std::size_t dimension : {2u, 3u})
      for (const Real amplitude : {Real(0), Real(1)})
      {
        CellDeformation deformation(dimension);
        Math::SpatialMatrix<Real> F =
          Math::SpatialMatrix<Real>::Identity(dimension, dimension);
        F(0, 0) = Real(1.4);
        F(1, 1) = Real(0.7);
        deformation.setDeformationGradient(F);
        Math::SpatialMatrix<Real> inner(dimension, dimension);
        inner.setZero();
        inner(0, 0) = amplitude * Real(0.5);
        inner(1, 1) = amplitude * Real(-0.9);
        Math::SpatialMatrix<Real> direction =
          Math::SpatialMatrix<Real>::Identity(dimension, dimension);
        direction(0, 1) = Real(0.2);
        SWIFT::Parameters p;
        p.model.distortion = Real(1.4);
        constexpr Real mu = Real(0.02), eps = Real(1e-6);
        const Real rowJ = -deformation.getJacobianAction(direction);
        const Real rowQ = deformation.getRelativeDistortionAction(direction);
        const auto scalarEnergy = [&](Real slack, Real delta) {
          const Real violation = std::max(Real(0), Real(1) - slack / delta);
          return mu * Real(0.5) * violation * violation;
        };
        const auto energy = [&](const Math::SpatialMatrix<Real>& v) {
          SWIFT::HingeState s(deformation, v, p, mu);
          const Real expected = scalarEnergy(s.getJacobianSlack(),
                                  p.model.qualityGuard * (Real(1) - p.model.jacobian)) +
            scalarEnergy(s.getDistortionSlack(),
              p.model.qualityGuard * (p.model.distortion - Real(1)));
          EXPECT_NEAR(s.getEnergy(p, mu), expected, Real(1e-12));
          return s.getEnergy(p, mu);
        };
        const auto derivative = [&](const Math::SpatialMatrix<Real>& v) {
          SWIFT::HingeState s(deformation, v, p, mu);
          return (s.getJacobianHessian() * s.getJacobianAction() - s.getJacobianForce()) *
            rowJ +
            (s.getDistortionHessian() * s.getDistortionAction() -
              s.getDistortionForce()) *
            rowQ;
        };
        SWIFT::HingeState s(deformation, inner, p, mu);
        ASSERT_TRUE(s.isAdmissible());
        const Math::SpatialMatrix<Real> plus(inner + eps * direction),
          minus(inner - eps * direction);
        EXPECT_NEAR(derivative(inner), (energy(plus) - energy(minus)) / (Real(2) * eps),
          Real(1e-7));
        EXPECT_NEAR(
          s.getJacobianHessian() * rowJ * rowJ + s.getDistortionHessian() * rowQ * rowQ,
          (derivative(plus) - derivative(minus)) / (Real(2) * eps), Real(1e-7));
      }
  }

  TEST(Rodin_Adaptation_SWIFTHingeState, AffineHingeAssembledTangentP1P2)
  {
    using namespace Geometry;
    using namespace Variational;
    const auto check = []<std::size_t Order>() {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
      mesh.scale(Real(0.5));
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);
      auto fes = [&]() {
        if constexpr (Order == 1)
          return P1<Math::SpatialVector<Real>, LocalMesh>(mesh, 2);
        else
          return H1(std::integral_constant<std::size_t, Order>{}, mesh, std::size_t{2});
      }();
      TrialFunction trial(fes);
      TestFunction test(fes);
      GridFunction current(fes), inner(fes), direction(fes);
      current = VectorFunction(std::size_t{2}, [](const Point& p) {
        return Math::SpatialVector<Real>{
          Real(0.4) * p.getCoordinates()(0), Real(-0.3) * p.getCoordinates()(1)};
      });
      inner = VectorFunction(std::size_t{2}, [](const Point& p) {
        return Math::SpatialVector<Real>{
          Real(0.5) * p.getCoordinates()(0), Real(-0.9) * p.getCoordinates()(1)};
      });
      direction = VectorFunction(std::size_t{2}, [](const Point& p) {
        return Math::SpatialVector<Real>{
          p.getCoordinates()(0) + Real(0.2) * p.getCoordinates()(1),
          p.getCoordinates()(1)};
      });
      SWIFT::Parameters p;
      p.model.distortion = Real(1.4);
      BilinearForm metric(trial, test);
      LinearForm force(test);
      const auto residual = [&]() {
        metric = SWIFT::HingeMetric(trial, test, current, inner, p, Real(0.02));
        force = SWIFT::HingeForce(test, current, inner, p, Real(0.02));
        metric.assemble();
        force.assemble();
        return Math::Vector<Real>(
          metric.getOperator() * inner.getData() - force.getVector());
      };
      residual();
      const Math::Vector<Real> tangent(metric.getOperator() * direction.getData());
      const Math::Vector<Real> original(inner.getData());
      constexpr Real eps = Real(1e-6);
      inner.getData() = original + eps * direction.getData();
      const auto plus = residual();
      inner.getData() = original - eps * direction.getData();
      const auto minus = residual();
      const Math::Vector<Real> difference((plus - minus) / (Real(2) * eps));
      ASSERT_GT(tangent.norm(), Real(0));
      EXPECT_LT((difference - tangent).norm() / tangent.norm(), Real(1e-6));
    };
    check.template operator()<1>();
    check.template operator()<2>();
  }

  TEST(Rodin_Adaptation_SWIFTHingeState, NegativeAffineSlackRemainsFinite)
  {
    CellDeformation outer(2);
    Math::SpatialMatrix<Real> increment = -Math::SpatialMatrix<Real>::Identity(2, 2);
    SWIFT::Parameters parameters;
    SWIFT::HingeState state(outer, increment, parameters, Real(0.02));
    EXPECT_TRUE(state.isAdmissible());
    EXPECT_LT(state.getJacobianSlack(), Real(0));
    EXPECT_TRUE(std::isfinite(state.getEnergy(parameters, Real(0.02))));
    EXPECT_GT(state.getJacobianHessian(), Real(0));
    EXPECT_EQ(state.getEnergy(parameters, Real(0)), Real(0));
    increment(0, 0) = Real(-1);
    increment(1, 1) = Real(1);
    outer.setDeformationGradient(increment);
    SWIFT::HingeState inverted(outer, increment, parameters, Real(0.02));
    EXPECT_FALSE(inverted.isAdmissible());
    EXPECT_TRUE(std::isinf(inverted.getEnergy(parameters, Real(0.02))));
  }

}
