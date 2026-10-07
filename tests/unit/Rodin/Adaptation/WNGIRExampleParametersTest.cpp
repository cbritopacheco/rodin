/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <utility>
#include "Rodin/Adaptation/WNGIR/Admissibility.h"
#include "Rodin/Adaptation/WNGIR/HingeProblem.h"
#include "Rodin/Variational.h"
#include "../../../../examples/WNGIRExampleParameters.h"

using namespace Rodin;

TEST(Rodin_Adaptation_WNGIRExampleParameters, CanonicalNamesMapToHierarchicalParameters)
{
  std::vector<std::string> arguments = {"example", "--wngir-fit=2",
    "--wngir-distribution-deviatoric=0.0003", "--wngir-distribution-divergence=0.0004",
    "--wngir-hinge=10", "--wngir-jacobian=0.02", "--wngir-distortion=5",
    "--wngir-outer-iterations=20", "--wngir-inner-iterations=10",
    "--wngir-linear-solver=cg", "--wngir-linear-relative-tolerance=1e-9"};
  std::vector<char*> argv;
  for (auto& argument : arguments)
    argv.push_back(argument.data());
  const auto p = Examples::makeWNGIRParameters(argv.size(), argv.data(), Real(0.1), 10);
  EXPECT_EQ(p.model.h, Real(0.1));
  EXPECT_EQ(p.model.fit, Real(2));
  EXPECT_EQ(p.model.distribution.deviatoric, Real(0.0003));
  EXPECT_EQ(p.model.distribution.divergence, Real(0.0004));
  EXPECT_EQ(p.model.hinge, Real(10));
  EXPECT_EQ(p.model.jacobian, Real(0.02));
  EXPECT_EQ(p.model.distortion, Real(5));
  EXPECT_EQ(p.convergence.iterations.outer, 20u);
  EXPECT_EQ(p.convergence.iterations.inner, 10u);
  EXPECT_EQ(p.linear.solver, Adaptation::WNGIRParameters::LinearSolver::CG);
  EXPECT_EQ(p.convergence.tolerance.linearRelative, Real(1e-9));
}

TEST(Rodin_Adaptation_WNGIRExampleParameters, RemovedOptionsAreRejected)
{
  for (const char* option : {"--wngir-volume-gauge=1", "--wngir-omega-min=0.1",
         "--wngir-jls=0.01", "--wngir-kappa-d=1", "--wngir-primal-barrier-iterations=15"})
  {
    std::string name = "example", argument = option;
    char* argv[] = {name.data(), argument.data()};
    EXPECT_THROW(Examples::makeWNGIRParameters(2, argv, Real(0.1), 10), Alert::Exception);
  }
}

TEST(Rodin_Adaptation_WNGIRHingeProblem, CenteringEntersTheStationarityResidual)
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
  Adaptation::WNGIRHingeProblem problem(trial, test);
  problem = Integral(Dot(trial, test));
  problem.setState(state).setCentering(centering).assemble();
  const auto& system = problem.getLinearSystem();
  const Math::Vector<Real> expected = -system.getOperator() * state.getData() +
    centering * (centering.transpose() * state.getData());
  EXPECT_LT((system.getVector() - expected).norm(), Real(1e-12));
}

TEST(Rodin_Adaptation_WNGIRAdmissibility, SamplingIsReadOnly)
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
  const auto report = Adaptation::evaluateWNGIRAdmissibilitySampled(
    std::as_const(displacement), Real(0.01));
  EXPECT_NEAR(report.minJ, Real(0.81), Real(1e-12));
  EXPECT_NEAR(report.maxQRel, Real(1), Real(1e-12));
  EXPECT_EQ(report.inadmissibleCount, 0u);
  EXPECT_EQ((displacement.getData() - before).norm(), Real(0));

  P1<Math::SpatialVector<Real>, LocalMesh> wrongSpace(mesh, 3);
  GridFunction wrongDimension(wrongSpace);
  EXPECT_THROW(Adaptation::evaluateWNGIRAdmissibilitySampled(
    wrongDimension, Real(0.01)), Alert::Exception);
}

TEST(Rodin_Adaptation_WNGIRReport, InnerFailuresUseCanonicalNames)
{
  Adaptation::WNGIRReport report;
  report.reason = Adaptation::WNGIRReport::Reason::InnerLineSearchFailure;
  EXPECT_STREQ(report.getReasonString(), "inner-line-search-failure");
  report.reason = Adaptation::WNGIRReport::Reason::InnerIterationLimit;
  EXPECT_STREQ(report.getReasonString(), "inner-iteration-limit");
}
