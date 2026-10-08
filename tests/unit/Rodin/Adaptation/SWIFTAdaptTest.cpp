/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include "Rodin/Adaptation.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace SWIFT = Rodin::Adaptation::SWIFT;

namespace Rodin::Tests::Unit
{
  namespace
  {
    constexpr Attribute Interface = 10;
    constexpr Attribute Fixed = 20;

    LocalMesh makeMesh()
    {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {5, 5});
      mesh.scale(Real(0.25));
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.getConnectivity().compute(1, 2);
      for (auto face = mesh.getFace(); face; ++face)
      {
        bool interface = true, left = true, right = true;
        for (const Index vertex : face->getVertices())
        {
          const Real x = mesh.getVertexCoordinates(vertex)(0);
          interface &= x == Real(0.5);
          left &= x == Real(0);
          right &= x == Real(1);
        }
        if (interface)
          mesh.setAttribute({1, face->getIndex()}, Interface);
        if (left || right)
          mesh.setAttribute({1, face->getIndex()}, Fixed);
      }
      return mesh;
    }

    struct MPIMesh
    {
      const Rodin::Context::MPI& getContext() const;
    };

    TEST(Rodin_Adaptation_SWIFTAdapt, ContextSelection)
    {
      static_assert(SWIFT::Adapt<LocalMesh>::isSupported);
      static_assert(!SWIFT::Adapt<MPIMesh>::isSupported);
      static_assert(!std::is_constructible_v<SWIFT::Adapt<MPIMesh>, MPIMesh&>);
    }

    TEST(Rodin_Adaptation_SWIFTAdapt, MatchesExplicitProblemAndDoesNotApplyTwice)
    {
      auto mesh = makeMesh();
      auto reference = makeMesh();
      SWIFT::Adapt adapt(mesh);
      SWIFT::Parameters parameters;
      parameters.model.h = Real(0.25);
      parameters.linear.solver = SWIFT::Parameters::LinearSolver::SparseLU;
      parameters.convergence.tolerance.geometric = Real(1e-7);
      adapt.setParameters(parameters).setInterfaceAttribute(Interface);
      auto& trial = adapt.getTrialFunction();
      adapt.getProblem() += DirichletBC(trial, VectorFunction(Real(0), Real(0))).on(Fixed);

      P1<Math::SpatialVector<Real>, LocalMesh> space(reference, 2);
      TrialFunction u(space);
      TestFunction v(space);
      SWIFT::Problem problem(u, v);
      problem.setParameters(parameters).setInterfaceAttribute(Interface);
      problem += DirichletBC(u, VectorFunction(Real(0), Real(0))).on(Fixed);
      RealFunction phi([](const Point& p) { return p.x() - Real(0.55); });
      VectorFunction gradient(Real(1), Real(0));
      const auto expected = problem.solve(phi, gradient);
      const auto report = adapt.execute(phi, gradient);
      ASSERT_TRUE(report.qualityBudgetSatisfied);
      EXPECT_EQ(report.reason, expected.reason);
      EXPECT_EQ(report.iterations, expected.iterations);
      EXPECT_NEAR(report.geometricSup, expected.geometricSup, Real(1e-12));
      reference.displace(u.getSolution());
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        EXPECT_LT((mesh.getVertexCoordinates(vertex) -
          reference.getVertexCoordinates(vertex)).norm(), Real(1e-12));

      const auto again = adapt.execute(phi, gradient);
      EXPECT_TRUE(again.geometricTargetReached);
      EXPECT_EQ(again.iterations, 0u);
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        EXPECT_LT((mesh.getVertexCoordinates(vertex) -
          reference.getVertexCoordinates(vertex)).norm(), Real(1e-12));
    }

    TEST(Rodin_Adaptation_SWIFTAdapt, EmptyInterfaceDoesNotMoveMesh)
    {
      auto mesh = makeMesh();
      const auto reference = makeMesh();
      SWIFT::Adapt adapt(mesh);
      SWIFT::Parameters parameters;
      parameters.model.h = Real(0.25);
      adapt.setParameters(parameters).setInterfaceAttribute(999);
      RealFunction phi([](const Point& p) { return p.x() - Real(0.55); });
      VectorFunction gradient(Real(1), Real(0));
      const auto report = adapt.execute(phi, gradient);
      EXPECT_FALSE(report.qualityBudgetSatisfied);
      EXPECT_STREQ(report.getReasonString(), "empty-interface");
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        EXPECT_EQ((mesh.getVertexCoordinates(vertex) -
          reference.getVertexCoordinates(vertex)).norm(), Real(0));
    }

    TEST(Rodin_Adaptation_SWIFTAdapt, AppliesQualityValidBestEffort)
    {
      auto mesh = makeMesh();
      const auto reference = makeMesh();
      SWIFT::Adapt adapt(mesh);
      SWIFT::Parameters parameters;
      parameters.model.h = Real(0.25);
      parameters.linear.solver = SWIFT::Parameters::LinearSolver::SparseLU;
      parameters.convergence.iterations.outer = 1;
      parameters.convergence.tolerance.geometric = Real(1e-30);
      adapt.setParameters(parameters).setInterfaceAttribute(Interface);
      adapt.getProblem() += DirichletBC(adapt.getTrialFunction(),
        VectorFunction(Real(0), Real(0))).on(Fixed);
      RealFunction phi([](const Point& p) { return p.x() - Real(0.55); });
      VectorFunction gradient(Real(1), Real(0));
      const auto report = adapt.execute(phi, gradient);
      ASSERT_TRUE(report.qualityBudgetSatisfied);
      EXPECT_FALSE(report.geometricTargetReached);
      Real motion = 0;
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        motion = std::max(motion, (mesh.getVertexCoordinates(vertex) -
          reference.getVertexCoordinates(vertex)).norm());
      EXPECT_GT(motion, Real(0));
    }

    TEST(Rodin_Adaptation_SWIFTAdapt, ThreeDimensionalInitialFit)
    {
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
      mesh.scale(Real(0.5));
      mesh.getConnectivity().compute(3, 2);
      mesh.getConnectivity().compute(2, 0);
      mesh.getConnectivity().compute(2, 3);
      for (auto face = mesh.getFace(); face; ++face)
      {
        bool marked = true;
        for (const Index vertex : face->getVertices())
          marked &= mesh.getVertexCoordinates(vertex)(0) == Real(0.5);
        if (marked)
          mesh.setAttribute({2, face->getIndex()}, Interface);
      }
      SWIFT::Adapt adapt(mesh);
      SWIFT::Parameters parameters;
      parameters.model.h = Real(0.5);
      adapt.setParameters(parameters).setInterfaceAttribute(Interface);
      RealFunction phi([](const Point& p) { return p.x() - Real(0.5); });
      VectorFunction gradient(Real(1), Real(0), Real(0));
      const auto report = adapt.execute(phi, gradient);
      EXPECT_TRUE(report.geometricTargetReached);
      EXPECT_EQ(report.iterations, 0u);
      EXPECT_NEAR(report.geometricSup, Real(0), Real(1e-14));
    }

    TEST(Rodin_Adaptation_SWIFTAdapt, RejectsHigherOrderTransformations)
    {
      auto mesh = makeMesh();
      RealH1Element<2> element(Polytope::Type::Triangle);
      PointCloud points(2, element.getCount());
      const auto& transformation = mesh.getPolytopeTransformation(2, 0);
      for (size_t node = 0; node < element.getCount(); ++node)
      {
        Math::SpatialPoint point;
        transformation.transform(point, element.getNode(node));
        points.col(node) = point;
      }
      mesh.setPolytopeTransformation({2, 0},
        new ParametricTransformation(points, element));
      EXPECT_THROW((SWIFT::Adapt(mesh)), Alert::Exception);
    }
  }
}
