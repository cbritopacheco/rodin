/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <limits>
#include <utility>
#include "Rodin/Adaptation/SWIFT/Admissibility.h"
#include "Rodin/Adaptation/SWIFT/HingeProblem.h"
#include "Rodin/Variational.h"
#include "../../../../examples/Adaptation/SWIFT/SWIFTExampleParameters.h"

using namespace Rodin;

TEST(Rodin_Adaptation_SWIFTExampleParameters, CalibratedDefaultsAreInherited)
{
  std::string name = "example";
  char* argv[] = {name.data()};
  const auto p = Examples::makeSWIFTParameters(1, argv, Real(0.1), 10);
  EXPECT_EQ(p.model.fit, Real(1));
  EXPECT_EQ(p.model.distribution.deviatoric, Real(1e-4));
  EXPECT_EQ(p.model.distribution.divergence, Real(1e-2));
  EXPECT_EQ(p.model.hinge, Real(10));
  EXPECT_EQ(p.convergence.iterations.outer, 30u);
  EXPECT_EQ(p.convergence.iterations.inner, 15u);
  EXPECT_EQ(p.quadrature.getSurfaceOrder(1), 8u);
  EXPECT_EQ(p.quadrature.getVolumeOrder(1), 2u);
  EXPECT_EQ(p.quadrature.getQualityOrder(1), 2u);
  EXPECT_EQ(p.quadrature.getSurfaceOrder(2), 12u);
  EXPECT_EQ(p.quadrature.getVolumeOrder(2), 8u);
  EXPECT_EQ(p.quadrature.getQualityOrder(2), 16u);
  EXPECT_EQ(Adaptation::SWIFT::Parameters::Quadrature::getValidationOrder(1), 32u);
}

TEST(Rodin_Adaptation_SWIFTExampleParameters, CanonicalNamesMapToHierarchicalParameters)
{
  std::vector<std::string> arguments = {"example", "--swift-fit=2",
    "--swift-distribution-deviatoric=0.0003", "--swift-distribution-divergence=0.0004",
    "--swift-hinge=10", "--swift-jacobian=0.02", "--swift-distortion=5",
    "--swift-outer-iterations=20", "--swift-inner-iterations=10",
    "--swift-linear-solver=cg", "--swift-linear-relative-tolerance=1e-9",
    "--quad-order=12", "--surface-quadrature-order=8", "--volume-quadrature-order=2",
    "--quality-validation-order=16", "--geometric-validation-order=32"};
  std::vector<char*> argv;
  for (auto& argument : arguments)
    argv.push_back(argument.data());
  const auto p = Examples::makeSWIFTParameters(argv.size(), argv.data(), Real(0.1), 10);
  EXPECT_EQ(p.model.h, Real(0.1));
  EXPECT_EQ(p.model.fit, Real(2));
  EXPECT_EQ(p.model.distribution.deviatoric, Real(0.0003));
  EXPECT_EQ(p.model.distribution.divergence, Real(0.0004));
  EXPECT_EQ(p.model.hinge, Real(10));
  EXPECT_EQ(p.model.jacobian, Real(0.02));
  EXPECT_EQ(p.model.distortion, Real(5));
  EXPECT_EQ(p.convergence.iterations.outer, 20u);
  EXPECT_EQ(p.convergence.iterations.inner, 10u);
  EXPECT_EQ(p.linear.solver, Adaptation::SWIFT::Parameters::LinearSolver::CG);
  EXPECT_EQ(p.convergence.tolerance.linearRelative, Real(1e-9));
  EXPECT_EQ(p.quadrature.order, 12u);
  EXPECT_EQ(p.quadrature.surface, 8u);
  EXPECT_EQ(p.quadrature.volume, 2u);
  EXPECT_EQ(p.quadrature.quality, 16u);
  EXPECT_EQ(p.quadrature.validation, 32u);
}

TEST(Rodin_Adaptation_SWIFTExampleParameters, RemovedOptionsAreRejected)
{
  for (const char* option : {"--swift-volume-gauge=1", "--swift-omega-min=0.1",
         "--swift-jls=0.01", "--swift-kappa-d=1", "--swift-primal-barrier-iterations=15"})
  {
    std::string name = "example", argument = option;
    char* argv[] = {name.data(), argument.data()};
    EXPECT_THROW(Examples::makeSWIFTParameters(2, argv, Real(0.1), 10), Alert::Exception);
  }
}

TEST(Rodin_Adaptation_SWIFTExampleParameters, PreviousMethodOptionsAreRejected)
{
  std::string name = "example", argument = "--wngir-fit=2";
  char* argv[] = {name.data(), argument.data()};
  EXPECT_THROW(Examples::makeSWIFTParameters(2, argv, Real(0.1), 10), Alert::Exception);
}

TEST(Rodin_Adaptation_SWIFTHingeProblem, CenteringEntersTheStationarityResidual)
{
  using namespace Geometry;
  using namespace Variational;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 0);
  P1<Math::SpatialVector<Real>, LocalMesh> space(mesh, 2);
  TrialFunction trial(space);
  TestFunction test(space);
  GridFunction state(space);
  state = VectorFunction(Real(0.2), Real(0.4));
  Math::Matrix<Real> centering =
    Math::Matrix<Real>::Constant(space.getSize(), 1, Real(0.01));
  Adaptation::SWIFT::HingeProblem problem(trial, test);
  problem = Integral(Dot(trial, test));
  problem.setState(state).setCentering(centering).assemble();
  const auto& system = problem.getLinearSystem();
  const Math::Vector<Real> expected = -system.getOperator() * state.getData() +
    centering * (centering.transpose() * state.getData());
  EXPECT_LT((system.getVector() - expected).norm(), Real(1e-12));
}

TEST(Rodin_Adaptation_SWIFTAdmissibility, SamplingIsReadOnly)
{
  using namespace Geometry;
  using namespace Variational;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 0);
  P1<Math::SpatialVector<Real>, LocalMesh> space(mesh, 2);
  GridFunction displacement(space);
  displacement = VectorFunction(size_t(2), [](const Point& point) {
    Math::SpatialVector<Real> value(2);
    value(0) = Real(0.2) - Real(0.1) * point.x();
    value(1) = Real(0.3) - Real(0.1) * point.y();
    return value;
  });
  const Math::Vector<Real> before = displacement.getData();
  const auto report = Adaptation::SWIFT::evaluateAdmissibility(
    std::as_const(displacement), Real(0.01));
  EXPECT_NEAR(report.minJ, Real(0.81), Real(1e-12));
  EXPECT_NEAR(report.maxQRel, Real(1), Real(1e-12));
  EXPECT_EQ(report.inadmissibleCount, 0u);
  EXPECT_EQ((displacement.getData() - before).norm(), Real(0));

  P1<Math::SpatialVector<Real>, LocalMesh> wrongSpace(mesh, 3);
  GridFunction wrongDimension(wrongSpace);
  EXPECT_THROW(Adaptation::SWIFT::evaluateAdmissibility(
    wrongDimension, Real(0.01)), Alert::Exception);
}

TEST(Rodin_Adaptation_SWIFTAdmissibility, RejectsNonfiniteGeometry)
{
  using namespace Geometry;
  using namespace Variational;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 0);
  P1<Math::SpatialVector<Real>, LocalMesh> space(mesh, 2);
  GridFunction displacement(space);
  for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
         std::numeric_limits<Real>::infinity()})
  {
    displacement.setData(Math::Vector<Real>::Constant(space.getSize(), invalid));
    const auto report = Adaptation::SWIFT::evaluateAdmissibility(
      std::as_const(displacement), Real(0.01));
    EXPECT_GT(report.inadmissibleCount, 0u);
  }
}

TEST(Rodin_Adaptation_SWIFTReport, InnerFailuresUseCanonicalNames)
{
  Adaptation::SWIFT::Report report;
  report.reason = Adaptation::SWIFT::Report::Reason::InnerLineSearchFailure;
  EXPECT_STREQ(report.getReasonString(), "inner-line-search-failure");
  report.reason = Adaptation::SWIFT::Report::Reason::InnerIterationLimit;
  EXPECT_STREQ(report.getReasonString(), "inner-iteration-limit");
}
