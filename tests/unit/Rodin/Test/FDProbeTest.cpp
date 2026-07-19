/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Test/FDProbe.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  const char* getName(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Point:
        return "Point";
      case Polytope::Type::Segment:
        return "Segment";
      case Polytope::Type::Triangle:
        return "Triangle";
      case Polytope::Type::Quadrilateral:
        return "Quadrilateral";
      case Polytope::Type::Tetrahedron:
        return "Tetrahedron";
      case Polytope::Type::Pyramid:
        return "Pyramid";
      case Polytope::Type::Hexahedron:
        return "Hexahedron";
      case Polytope::Type::Wedge:
        return "Wedge";
    }
    return "Unknown";
  }

  LocalMesh makeUnitMesh(Polytope::Type type)
  {
    LocalMesh mesh;
    switch (type)
    {
      case Polytope::Type::Segment:
        mesh = LocalMesh::UniformGrid(type, {5});
        mesh.scale(1.0 / 4.0);
        mesh.getConnectivity().compute(0, 1);
        break;
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        mesh = LocalMesh::UniformGrid(type, {4, 4});
        mesh.scale(1.0 / 3.0);
        mesh.getConnectivity().compute(1, 2);
        break;
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
        mesh = LocalMesh::UniformGrid(type, {3, 3, 3});
        mesh.scale(1.0 / 2.0);
        mesh.getConnectivity().compute(2, 3);
        break;
      case Polytope::Type::Point:
        mesh = LocalMesh::Builder()
                 .initialize(1)
                 .nodes(1)
                 .vertex({0.0})
                 .polytope(Polytope::Type::Point, {0})
                 .finalize();
        break;
    }
    return mesh;
  }

  LocalMesh makeUnitSquareMesh()
  {
    return makeUnitMesh(Polytope::Type::Triangle);
  }

  template <class State>
  void setSmoothState(State& state)
  {
    state = [](const Geometry::Point& p) {
      const auto& x = p.getPhysicalCoordinates();
      Real value = 0.3;
      for (Index d = 0; d < x.size(); ++d)
      {
        const Real scale = static_cast<Real>(d + 1);
        value += 0.04 * scale * x(d);
        value += 0.01 * std::sin(scale * Math::Constants::pi() * x(d));
      }
      return value;
    };
  }

  template <class ProblemType>
  Rodin::Test::FDProbeReport runNonlinearScalarProbe(ProblemType& problem)
  {
    Rodin::Test::FDProbe probe(problem);
    return probe.test(1e-6);
  }

  template <class ProblemType>
  auto makeInteriorVertexDirection(const LocalMesh& mesh, ProblemType& problem)
  {
    auto direction = problem.getTrialFunction().getSolution().getData();
    direction.setZero();

    const auto& fes = problem.getTrialFunction().getFiniteElementSpace();
    for (auto it = mesh.getVertex(); it; ++it)
    {
      const auto& x = it->getCoordinates();
      if (x(0) > 0.0 && x(0) < 1.0 && x(1) > 0.0 && x(1) < 1.0)
      {
        const auto& dofs = fes.getDOFs(0, it->getIndex());
        for (decltype(dofs.size()) i = 0; i < dofs.size(); ++i)
          direction[dofs[i]] = 1.0;
      }
    }

    const Real norm = direction.norm();
    if (norm > 0.0)
      direction *= 1.0 / norm;
    return direction;
  }
}

TEST(Rodin_Test_FDProbe, NonlinearScalarProblemPasses)
{
  LocalMesh mesh = makeUnitSquareMesh();
  P1 Vh(mesh);
  TrialFunction du(Vh);
  TestFunction v(Vh);

  setSmoothState(du.getSolution());

  const auto& u = du.getSolution();
  Problem problem(du, v);
  problem = 2.0 * Integral(u * du, v) + Integral(u * u, v);

  Rodin::Test::FDProbe probe(problem);
  const auto report = probe.test(1e-6);

  EXPECT_LT(report.relativeError, 1e-6)
    << "absoluteError = " << report.absoluteError
    << ", tangentNorm = " << report.tangentNorm
    << ", finiteDifferenceNorm = " << report.finiteDifferenceNorm;
}

TEST(Rodin_Test_FDProbe, NonlinearScalarProblemPassesAcrossGeometries)
{
  for (const Polytope::Type type : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    SCOPED_TRACE(getName(type));

    LocalMesh mesh = makeUnitMesh(type);
    P1 Vh(mesh);
    TrialFunction du(Vh);
    TestFunction v(Vh);

    setSmoothState(du.getSolution());

    const auto& u = du.getSolution();
    Problem problem(du, v);
    problem = 2.0 * Integral(u * du, v) + Integral(u * u, v);

    const auto report = runNonlinearScalarProbe(problem);
    EXPECT_LT(report.relativeError, 1e-6)
      << "absoluteError = " << report.absoluteError
      << ", tangentNorm = " << report.tangentNorm
      << ", finiteDifferenceNorm = " << report.finiteDifferenceNorm;
  }
}

TEST(Rodin_Test_FDProbe, LinearSystemSpecializationPasses)
{
  Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>> system;

  std::vector<Eigen::Triplet<Real>> triplets;
  triplets.emplace_back(0, 0, 4.0);
  triplets.emplace_back(0, 1, -1.0);
  triplets.emplace_back(1, 0, 2.0);
  triplets.emplace_back(1, 1, 3.0);
  triplets.emplace_back(1, 2, 0.5);
  triplets.emplace_back(2, 1, -2.0);
  triplets.emplace_back(2, 2, 5.0);
  system.getOperator().resize(3, 3);
  system.getOperator().setFromTriplets(triplets.begin(), triplets.end());

  system.getVector().resize(3);
  system.getVector() << 1.0, -2.0, 0.5;

  system.getSolution().resize(3);
  system.getSolution() << 0.1, -0.2, 0.3;

  Rodin::Test::FDProbe probe(system);
  const auto report = probe.test(1e-6);

  EXPECT_LT(report.relativeError, 1e-10)
    << "absoluteError = " << report.absoluteError
    << ", tangentNorm = " << report.tangentNorm
    << ", finiteDifferenceNorm = " << report.finiteDifferenceNorm;
}

TEST(Rodin_Test_FDProbe, WrongTangentFails)
{
  LocalMesh mesh = makeUnitSquareMesh();
  P1 Vh(mesh);
  TrialFunction du(Vh);
  TestFunction v(Vh);

  setSmoothState(du.getSolution());

  const auto& u = du.getSolution();
  Problem problem(du, v);
  problem = Integral(u * du, v) + Integral(u * u, v);

  Rodin::Test::FDProbe probe(problem);
  const auto report = probe.test(1e-6);

  EXPECT_GT(report.relativeError, 1e-3);
}

TEST(Rodin_Test_FDProbe, DirichletProblemPassesWithAdmissibleDirection)
{
  LocalMesh mesh = makeUnitSquareMesh();
  P1 Vh(mesh);
  TrialFunction du(Vh);
  TestFunction v(Vh);

  setSmoothState(du.getSolution());

  const auto& u = du.getSolution();
  Problem problem(du, v);
  problem = 2.0 * Integral(u * du, v) + Integral(u * u, v) + DirichletBC(du, Zero());

  Rodin::Test::FDProbe probe(problem);
  const auto direction = makeInteriorVertexDirection(mesh, problem);
  ASSERT_GT(direction.norm(), 0.0);

  const auto report = probe.test(direction, 1e-6);
  EXPECT_LT(report.relativeError, 1e-6)
    << "absoluteError = " << report.absoluteError
    << ", tangentNorm = " << report.tangentNorm
    << ", finiteDifferenceNorm = " << report.finiteDifferenceNorm;
}
