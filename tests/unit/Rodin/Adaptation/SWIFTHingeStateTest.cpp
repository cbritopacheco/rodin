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
  TEST(Rodin_Adaptation_SWIFTHingeState, AdaptiveMeasuresPreserveMassAcrossGeometries)
  {
    using namespace Geometry;
    struct Gradient
    {
        bool zero = false;
        Math::SpatialMatrix<Real> getValue(const Variational::IntegrationPoint& ip) const
        {
          const auto dim = ip.getPoint().getPolytope().getDimension();
          Math::SpatialMatrix<Real> gradient(dim, dim);
          gradient.setZero();
          if (!zero)
          {
            gradient(0, 0) = Real(-1.5) * ip.getPoint().getReferenceCoordinates()(0);
            if (dim > 1)
              gradient(0, 1) = Real(3) * ip.getPoint().getReferenceCoordinates()(1);
          }
          return gradient;
        }
    };
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
    {
      const Polytope::Traits traits(geometry);
      Array<size_t> shape(traits.getDimension());
      shape.setConstant(2);
      auto mesh = LocalMesh::UniformGrid(geometry, shape);
      const auto cell = mesh.getCell(0);
      SWIFT::Parameters p;
      p.model.distortion = 2;
      p.quadrature.quality = 4;
      const SWIFT::QualitySamples samples(*cell, 2, p);
      Real measure = 0;
      size_t count = 0;
      samples.forEach([&](const Variational::IntegrationPoint&, Real weight) {
        measure += weight;
        ++count;
      });
      for (const bool zero : {false, true})
      {
        Real massJ = 0, massQ = 0;
        bool different = false;
        samples.forEachHinge(Gradient{true}, Gradient{zero},
          [&](const Variational::IntegrationPoint&, Real weightJ, Real weightQ) {
            EXPECT_TRUE(std::isfinite(weightJ));
            EXPECT_TRUE(std::isfinite(weightQ));
            EXPECT_GE(weightJ, Real(0.5) * measure / count);
            EXPECT_GE(weightQ, Real(0.5) * measure / count);
            massJ += weightJ;
            massQ += weightQ;
            different = different || std::abs(weightJ - weightQ) > Real(1e-12);
            if (zero)
            {
              EXPECT_DOUBLE_EQ(weightJ, measure / count);
              EXPECT_DOUBLE_EQ(weightQ, measure / count);
            }
          });
        EXPECT_NEAR(massJ, measure, Real(1e-12));
        EXPECT_NEAR(massQ, measure, Real(1e-12));
        if (!zero && traits.getDimension() > 1)
          EXPECT_TRUE(different);
      }
    }
  }


  TEST(Rodin_Adaptation_SWIFTHingeState, QualityQuadratureAcrossGeometries)
  {
    using namespace Geometry;
    for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge, Polytope::Type::Pyramid})
    {
      const Polytope::Traits traits(geometry);
      Array<size_t> shape(traits.getDimension());
      shape.setConstant(2);
      auto mesh = LocalMesh::UniformGrid(geometry, shape);
      SWIFT::Parameters parameters;
      parameters.quadrature.quality = 4;
      Real mass = 0, moment = 0;
      for (Index index = 0; index < mesh.getCellCount(); ++index)
      {
        const auto cell = mesh.getCell(index);
        const SWIFT::QualitySamples samples(*cell, 2, parameters);
        std::vector<bool> vertices(traits.getVertexCount(), false);
        samples.forEach([&](const Variational::IntegrationPoint& ip, Real weight) {
          EXPECT_GT(weight, Real(0));
          mass += weight;
          const Real x = ip.getPoint().getCoordinates()(0);
          moment += weight * x * x;
          for (size_t vertex = 0; vertex < vertices.size(); ++vertex)
          {
            const Point point(*cell, traits.getVertex(vertex));
            vertices[vertex] = vertices[vertex] ||
              (point.getCoordinates() - ip.getPoint().getCoordinates()).norm() < Real(1e-14);
          }
        });
        for (const bool included : vertices)
          EXPECT_TRUE(included);
      }
      EXPECT_NEAR(mass, Real(1), Real(1e-12)) << static_cast<int>(geometry);
      EXPECT_NEAR(moment, Real(1) / 3, Real(1e-12)) << static_cast<int>(geometry);
    }
  }

  TEST(Rodin_Adaptation_SWIFTHingeState, AffineHingeDerivatives2D3D)
  {
    for (const std::size_t dimension : {2u, 3u})
    {
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
      GridFunction current(fes), inner(fes), direction(fes), predictor(fes);
      current = VectorFunction(std::size_t{2}, [](const Point& p) {
        return Math::SpatialVector<Real>{
          Real(0.4) * p.getCoordinates()(0) +
            (Order > 1 ? Real(0.05) * p.getCoordinates()(0) * p.getCoordinates()(0)
                       : Real(0)),
          Real(-0.3) * p.getCoordinates()(1)};
      });
      inner = VectorFunction(std::size_t{2}, [](const Point& p) {
        return Math::SpatialVector<Real>{
          Real(0.5) * p.getCoordinates()(0), Real(-0.9) * p.getCoordinates()(1) +
            (Order > 1 ? Real(0.03) * p.getCoordinates()(1) * p.getCoordinates()(1)
                       : Real(0))};
      });
      direction = VectorFunction(std::size_t{2}, [](const Point& p) {
        return Math::SpatialVector<Real>{
          p.getCoordinates()(0) + Real(0.2) * p.getCoordinates()(1),
          p.getCoordinates()(1)};
      });
      predictor = inner;
      SWIFT::Parameters p;
      p.model.distortion = Real(1.4);
      BilinearForm metric(trial, test);
      LinearForm force(test);
      const auto residual = [&]() {
        metric = SWIFT::HingeMetric(trial, test, current, inner, predictor, p, Real(0.02));
        force = SWIFT::HingeForce(test, current, inner, predictor, p, Real(0.02));
        metric.assemble();
        force.assemble();
        return Math::Vector<Real>(
          metric.getOperator() * inner.getData() - force.getVector());
      };
      residual();
      const Math::Vector<Real> baseResidual(residual());
      const auto energy = [&]() {
        Real value = 0;
        auto currentJacobian = Jacobian(current);
        auto innerJacobian = Jacobian(inner);
        auto predictorJacobian = Jacobian(predictor);
        for (Index index = 0; index < mesh.getCellCount(); ++index)
        {
          const auto cell = mesh.getCell(index);
          const SWIFT::QualitySamples samples(*cell, Order, p);
          samples.forEachHinge(currentJacobian, predictorJacobian,
            [&](const IntegrationPoint& ip, Real weightJ, Real weightQ) {
            CellDeformation deformation(2);
            deformation.setDisplacementGradient(currentJacobian.getValue(ip));
            const SWIFT::HingeState state(
              deformation, innerJacobian.getValue(ip), p, Real(0.02));
            value += state.getEnergy(p, Real(0.02), weightJ, weightQ);
          });
        }
        return value;
      };
      const Math::Vector<Real> tangent(metric.getOperator() * direction.getData());
      const Math::Vector<Real> original(inner.getData());
      constexpr Real eps = Real(1e-6);
      inner.getData() = original + eps * direction.getData();
      const Real energyPlus = energy();
      const auto plus = residual();
      inner.getData() = original - eps * direction.getData();
      const Real energyMinus = energy();
      const auto minus = residual();
      const Math::Vector<Real> difference((plus - minus) / (Real(2) * eps));
      ASSERT_GT(tangent.norm(), Real(0));
      EXPECT_LT((difference - tangent).norm() / tangent.norm(), Real(1e-6));
      EXPECT_NEAR((energyPlus - energyMinus) / (Real(2) * eps),
        baseResidual.dot(direction.getData()), Real(1e-7));
    };
    check.template operator()<1>();
    check.template operator()<2>();
  }

  TEST(Rodin_Adaptation_SWIFTHingeState, P2VertexOnlyHingeIsAssembled)
  {
    using namespace Geometry;
    using namespace Variational;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(1, 2);
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, size_t{2});
    TrialFunction trial(fes);
    TestFunction test(fes);
    GridFunction current(fes), inner(fes);
    current = Math::SpatialVector<Real>{0, 0};
    inner = VectorFunction(size_t{2}, [](const Point& point) {
      const Real y = point.getCoordinates()(1);
      return Math::SpatialVector<Real>{0, Real(-0.48) * y * y};
    });
    SWIFT::Parameters parameters;
    parameters.quadrature.quality = 1;
    parameters.quadrature.volume = 1;
    auto jacobian = Jacobian(inner);
    Real mass = 0;
    bool vertexActive = false;
    for (Index index = 0; index < mesh.getCellCount(); ++index)
    {
      const auto cell = mesh.getCell(index);
      const SWIFT::QualitySamples samples(*cell, 2, parameters);
      samples.forEach([&](const IntegrationPoint& ip, Real weight) {
        EXPECT_GT(weight, Real(0));
        mass += weight;
        CellDeformation deformation(2);
        const SWIFT::HingeState state(deformation, jacobian.getValue(ip), parameters, 1);
        if (ip.getPoint().getCoordinates()(1) == Real(1))
          vertexActive = vertexActive || state.getJacobianHessian() > 0;
        else
          EXPECT_EQ(state.getJacobianHessian(), Real(0));
      });
    }
    EXPECT_TRUE(vertexActive);
    EXPECT_NEAR(mass, Real(1), Real(1e-12));
    BilinearForm metric(trial, test);
    LinearForm force(test);
    metric = SWIFT::HingeMetric(trial, test, current, inner, inner, parameters, 1);
    force = SWIFT::HingeForce(test, current, inner, inner, parameters, 1);
    metric.assemble();
    force.assemble();
    EXPECT_GT(metric.getOperator().norm(), Real(0));
    EXPECT_GT(force.getVector().norm(), Real(0));
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
