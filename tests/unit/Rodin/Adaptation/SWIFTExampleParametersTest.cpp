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
#include "../../../../experiments/swift_calibration/Parameters.h"
#include "../../../../examples/Adaptation/SWIFT/Options.h"

using namespace Rodin;

TEST(Rodin_Adaptation_SWIFTExampleParameters, ReconstructionNamedOptions)
{
  std::vector<std::string> arguments = {"SWIFT_ReconstructionP2", "--n", "8",
    "--dimension=3", "--lobes=4", "--amp=0.03", "--R0=0.2", "--phase=0.4", "--cx=0.4",
    "--cy=0.6", "--cz=0.3", "--output=results/test", "--swift-fit=2",
    "--swift-distribution-deviatoric=0.001", "--swift-distribution-divergence=0.01",
    "--swift-hinge=100", "--swift-jacobian=0.02", "--swift-distortion=8",
    "--swift-outer-iterations=20", "--swift-inner-iterations=10",
    "--swift-linear-solver=cg", "--swift-linear-threads=2",
    "--swift-linear-relative-tolerance=1e-9", "--swift-linear-iterations=500",
    "--quad-order=6", "--surface-quadrature-order=8", "--volume-quadrature-order=4",
    "--quality-validation-order=16", "--geometric-validation-order=32",
    "--swift-directional-newton=0", "--swift-trace", "--swift-quality-witness=1"};
  std::vector<char*> argv;
  for (auto& argument : arguments)
    argv.push_back(argument.data());
  const Examples::ReconstructionOptions options(argv.size(), argv.data());
  const auto& p = options.parameters;
  EXPECT_EQ(options.n, 8u);
  EXPECT_EQ(options.dimension, 3u);
  EXPECT_EQ(options.lobes, 4);
  EXPECT_EQ(options.amplitude, Real(0.03));
  EXPECT_EQ(options.radius, Real(0.2));
  EXPECT_EQ(options.phase, Real(0.4));
  EXPECT_EQ(options.cx, Real(0.4));
  EXPECT_EQ(options.cy, Real(0.6));
  EXPECT_EQ(options.cz, Real(0.3));
  EXPECT_EQ(options.output, "results/test");
  EXPECT_EQ(p.model.h, Real(1) / Real(7));
  EXPECT_EQ(p.model.fit, 2);
  EXPECT_EQ(p.model.distribution.deviatoric, Real(0.001));
  EXPECT_EQ(p.model.distribution.divergence, Real(0.01));
  EXPECT_EQ(p.model.hinge, 100);
  EXPECT_EQ(p.model.jacobian, Real(0.02));
  EXPECT_EQ(p.model.distortion, 8);
  EXPECT_EQ(p.convergence.iterations.outer, 20u);
  EXPECT_EQ(p.convergence.iterations.inner, 10u);
  EXPECT_EQ(p.linear.solver, Adaptation::SWIFT::Parameters::LinearSolver::CG);
  EXPECT_EQ(p.linear.threads, 2u);
  EXPECT_EQ(p.convergence.tolerance.linearRelative, Real(1e-9));
  EXPECT_EQ(p.convergence.iterations.linear, 500u);
  EXPECT_EQ(p.quadrature.order, 6u);
  EXPECT_EQ(p.quadrature.surface, 8u);
  EXPECT_EQ(p.quadrature.volume, 4u);
  EXPECT_EQ(p.quadrature.quality, 16u);
  EXPECT_EQ(p.quadrature.validation, 32u);
  EXPECT_FALSE(p.globalization.directionalNewton);
  EXPECT_TRUE(p.trace);
  EXPECT_TRUE(p.traceQualityWitness);
}

TEST(Rodin_Adaptation_SWIFTExampleParameters, ReconstructionRejectsMalformedOptions)
{
  for (const char* invalid :
    {"--n", "--n=", "--n=-1", "--n=3", "--n=8junk", "--n=184467440737095516160",
      "--dimension=4", "--swift-fit=nan", "--swift-fit=inf", "--swift-fit=1bad",
      "--swift-trace=2", "--swift-linear-solver=unknown", "--wngir-fit=1", "--unknown=1",
      "--lobes=2.5", "--amp=-0.1", "--R0=0.01", "16",
      "--swift-equal-hinge-weights", "--swift-stratified-hinge-weights",
      "--swift-critical-fraction=0.5", "--swift-adaptive-hinge-weights",
      "--swift-nonlinear-hinges"})
  {
    SCOPED_TRACE(invalid);
    std::string executable = "SWIFT_ReconstructionP1", argument = invalid;
    char* argv[] = {executable.data(), argument.data()};
    EXPECT_THROW(Examples::ReconstructionOptions(2, argv), Alert::Exception);
  }
}

TEST(Rodin_Adaptation_SWIFTExampleParameters, ReconstructionInheritsProductionDefaults)
{
  std::string executable = "SWIFT_ReconstructionP3";
  char* argv[] = {executable.data()};
  const Examples::ReconstructionOptions options(1, argv);
  const Adaptation::SWIFT::Parameters defaults;
  EXPECT_EQ(options.n, 16u);
  EXPECT_EQ(options.dimension, 2u);
  EXPECT_EQ(options.lobes, 0);
  EXPECT_EQ(options.parameters.model.h, Real(1) / Real(15));
  EXPECT_EQ(options.parameters.model.fit, defaults.model.fit);
  EXPECT_EQ(options.parameters.model.distribution.deviatoric,
    defaults.model.distribution.deviatoric);
  EXPECT_EQ(options.parameters.model.distribution.divergence,
    defaults.model.distribution.divergence);
  EXPECT_EQ(options.parameters.model.hinge, defaults.model.hinge);
  EXPECT_EQ(options.parameters.convergence.iterations.outer,
    defaults.convergence.iterations.outer);
  EXPECT_EQ(options.parameters.convergence.iterations.inner,
    defaults.convergence.iterations.inner);
}

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
  const auto report =
    Adaptation::SWIFT::evaluateAdmissibility(std::as_const(displacement), Real(0.01));
  EXPECT_NEAR(report.minJ, Real(0.81), Real(1e-12));
  EXPECT_NEAR(report.maxQRel, Real(1), Real(1e-12));
  EXPECT_EQ(report.inadmissibleCount, 0u);
  EXPECT_EQ((displacement.getData() - before).norm(), Real(0));

  P1<Math::SpatialVector<Real>, LocalMesh> wrongSpace(mesh, 3);
  GridFunction wrongDimension(wrongSpace);
  EXPECT_THROW(Adaptation::SWIFT::evaluateAdmissibility(wrongDimension, Real(0.01)),
    Alert::Exception);
}

TEST(Rodin_Adaptation_SWIFTAdmissibility, RejectsNonfiniteGeometry)
{
  using namespace Geometry;
  using namespace Variational;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 0);
  P1<Math::SpatialVector<Real>, LocalMesh> space(mesh, 2);
  GridFunction displacement(space);
  for (const Real invalid :
    {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
  {
    displacement.setData(Math::Vector<Real>::Constant(space.getSize(), invalid));
    const auto report =
      Adaptation::SWIFT::evaluateAdmissibility(std::as_const(displacement), Real(0.01));
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
