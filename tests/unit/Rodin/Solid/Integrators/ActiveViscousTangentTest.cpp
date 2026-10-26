/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Coupled active viscous tangent checks on nonaffine wedges.
#include <gtest/gtest.h>
#include <array>
#include "Rodin/Solid.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /**
   * @brief Checks the full four-field derivative, including rate and active couplings.
   * @param space Displacement finite-element space on nonaffine wedges.
   */
  template <class Space>
  void checkCoupledDerivative(const Space& space)
  {
    const auto& mesh = space.getMesh();
    P0 activeSpace(mesh);
    GridFunction u(space);
    GridFunction c(activeSpace), g(activeSpace), b(activeSpace);
    u = VectorFunction(3, [](const Point& x) {
      Math::SpatialVector<Real> value(3);
      value(0) = .02 * x(0) + .001 * x(0) * x(1);
      value(1) = -.01 * x(1);
      value(2) = .005 * x(2);
      return value;
    });
    c = -.03;
    g = 1.4;
    b = .8;
    TrialFunction du(space);
    TestFunction v(space);
    TrialFunction dc(activeSpace), dg(activeSpace), db(activeSpace);
    TestFunction qc(activeSpace), qg(activeSpace), qb(activeSpace);
    Solid::ActiveFiberLaw::Parameters active;
    active.stiffness = 20;
    active.damping = 3;
    active.destructionRate = .4;
    active.crossBridgeStiffness = 100;
    active.contractility = 20;
    Solid::ActiveContraction law(
      Solid::KelvinVoigt(Solid::DemirayGasser(Solid::DemirayGasser::Parameters{}), .2),
      Solid::ActiveFiberLaw(active));
    constexpr Real velocityCoefficient = 20;
    const auto input = [&](Solid::ConstitutivePoint& point) {
      const auto& x = point.getPoint()->get();
      Math::SpatialVector<Real> fiber(3);
      fiber.setZero();
      fiber(0) = 1;
      point.set<Solid::Tags::FiberDirection>(fiber);
      point.set<Solid::Tags::VelocityGradient>(
        Math::SpatialMatrix<Real>(.03 * Math::SpatialMatrix<Real>::Identity(3, 3) +
          velocityCoefficient * point.getKinematicState().getDisplacementGradient()));
      point.set<Solid::Tags::VelocityGradientDerivative>(velocityCoefficient);
      point.set<Solid::Tags::TimeStep>(.05);
      point.set<Solid::Tags::PreviousActiveExtension>(-.02);
      point.set<Solid::Tags::PreviousActiveGamma>(1.2);
      point.set<Solid::Tags::PreviousActiveBeta>(.6);
      point.set<Solid::Tags::PreviousFiberStrain>(.01);
      point.set<Solid::Tags::ElectricalActivation>(2);
      point.set<Solid::Tags::ActiveExtension>(c(x));
      point.set<Solid::Tags::ActiveGamma>(g(x));
      point.set<Solid::Tags::ActiveBeta>(b(x));
    };
    const auto evaluate = [&](const Point& point) {
      Solid::KinematicState state(3);
      state.setDisplacementGradient(Jacobian(u).getValue(point));
      Solid::ConstitutivePoint cp(point, state);
      input(cp);
      auto prepared = law;
      prepared.setConstitutivePoint(cp);
      return prepared;
    };
    constexpr size_t order = 8;
    auto work = Solid::InternalVirtualWork(law, u).setInput(input);
    work.setQuadratureOrder(order);
    Problem problem(du, dc, dg, db, v, qc, qg, qb);
    problem = work(du, v) + 2 * Integral(du, v) + 2 * Integral(u, v);
    const auto coupling = [&](bool momentum) {
      return MatrixFunction(3, 3, [&, momentum](const Point& point) {
        const auto prepared = evaluate(point);
        const Real factor = momentum ? .5 * prepared.getActiveResponse().dStressDc
                                     : prepared.getActiveResidual().strainDerivative[0];
        return Math::SpatialMatrix<Real>(factor *
          prepared.getKinematicState().getDeformationGradient() *
          prepared.getFiberKinematics().tensor());
      });
    };
    auto uc = Integral(dc * coupling(true), Jacobian(v));
    uc.setOrder(order);
    problem += uc;
    auto cu = Integral(Dot(coupling(false), Jacobian(du)), qc);
    cu.setOrder(order);
    problem += cu;
    const std::array unknowns{std::ref(dc), std::ref(dg), std::ref(db)};
    const std::array tests{std::ref(qc), std::ref(qg), std::ref(qb)};
    for (size_t row = 0; row < 3; ++row)
    {
      auto load = Integral(RealFunction([&, row](const Point& x) {
        return evaluate(x).getActiveResidual().value[row];
      }),
        tests[row].get());
      load.setOrder(order);
      problem += load;
      for (size_t column = 0; column < 3; ++column)
      {
        auto coefficient = RealFunction([&, row, column](const Point& x) {
          return evaluate(x).getActiveResidual().derivative(row, column);
        });
        auto block = Integral(coefficient * unknowns[column].get(), tests[row].get());
        block.setOrder(order);
        problem += block;
      }
    }
    const Index mechanics = u.getData().size(), size = c.getData().size();
    Math::Vector<Real> saved(mechanics + 3 * size);
    saved << u.getData(), c.getData(), g.getData(), b.getData();
    Math::Vector<Real> direction(saved.size());
    for (Index i = 0; i < static_cast<Index>(direction.size()); ++i)
      direction(i) = std::sin(Real(i) + 1);
    direction.normalize();
    problem.assemble();
    const Math::Vector<Real> expected =
      problem.getLinearSystem().getOperator() * direction;
    const auto scatter = [&](const Math::Vector<Real>& x) {
      u.getData() = x.head(mechanics);
      c.getData() = x.segment(mechanics, size);
      g.getData() = x.segment(mechanics + size, size);
      b.getData() = x.tail(size);
    };
    constexpr Real step = 1e-6;
    scatter(saved + step * direction);
    problem.assemble();
    const Math::Vector<Real> plus = problem.getLinearSystem().getVector();
    scatter(saved - step * direction);
    problem.assemble();
    const Math::Vector<Real> minus = problem.getLinearSystem().getVector();
    const Real relative =
      (expected + (plus - minus) / (2 * step)).norm() / expected.norm();
    EXPECT_LT(relative, 2e-7);
  }
}

TEST(Rodin_Solid_ActiveViscousTangent, P1OnNonaffineWedges)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, {2, 2, 2});
  mesh.setVertexCoordinates(mesh.getVertexCount() - 1, 1.07, 0);
  P1 space(mesh, 3);
  checkCoupledDerivative(space);
}

TEST(Rodin_Solid_ActiveViscousTangent, QuadraticOnNonaffineWedges)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, {2, 2, 2});
  mesh.setVertexCoordinates(mesh.getVertexCount() - 1, 1.07, 0);
  for (size_t dimension = 3; dimension > 0; --dimension)
    mesh.getConnectivity().compute(dimension, dimension - 1);
  H1<2, Math::SpatialVector<Real>> space(std::integral_constant<size_t, 2>{}, mesh, 3);
  checkCoupledDerivative(space);
}
