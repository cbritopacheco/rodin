/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Degree improvement for mixed Stokes velocity and pressure. */

#include "../../StokesProblem.h"
#include "../../CurvedGeometry.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::P::Stokes
{
  TEST(StokesSpaceTest, CoarseMeshesHaveTooFewFreeVelocityDOFs)
  {
    using namespace Variational;
    for (const auto geometry :
      {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
        Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(UniformGrid::getGeometryName(geometry));
      for (bool curved : {false, true})
      {
        SCOPED_TRACE(::testing::Message() << "curved=" << curved);
        auto mesh = UniformGrid(geometry).makeMesh(2);
        if (curved)
        {
          CurvedGeometry mapping(mesh);
          mapping.install<2>();
        }
        const auto dimension = mesh.getDimension();
        H1 velocitySpace(std::integral_constant<size_t, 2>{}, mesh, dimension);
        H1 pressureSpace(std::integral_constant<size_t, 1>{}, mesh);
        TrialFunction u(velocitySpace);
        auto boundary = DirichletBC(u, Zero());
        boundary.assemble();
        const auto& constrained = std::get<IndexMap<Real>>(boundary.getDOFs());
        ASSERT_LE(constrained.size(), velocitySpace.getSize());
        for (const auto& [index, value] : constrained)
          EXPECT_LT(index, velocitySpace.getSize());
        const auto freeVelocity = velocitySpace.getSize() - constrained.size();
        ASSERT_EQ(freeVelocity, dimension);
        ASSERT_EQ(pressureSpace.getSize(), dimension == 2 ? 4u : 8u);
        const auto zeroMeanPressure = pressureSpace.getSize() - 1;
        // Rank(B) <= dim(V_h^0) < dim(Q_h^0), independently of the map:
        // a nonconstant pressure null mode remains after fixing the mean.
        EXPECT_LT(freeVelocity, zeroMeanPressure);
      }
    }
  }

  class StokesPTest : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(StokesPTest, BothFieldsImproveAtEveryDegree)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const StokesProblem problem(
      mesh, StokesData(mesh.getDimension(), StokesData::Field::Smooth));
    const auto e2 = problem.solve<2>(), e3 = problem.solve<3>(), e4 = problem.solve<4>();
    ErrorHistory velocity, pressure;
    velocity.append(2, e2.velocity).append(3, e3.velocity).append(4, e4.velocity);
    pressure.append(2, e2.pressure).append(3, e3.pressure).append(4, e4.pressure);
    for (bool isVelocity : {true, false})
    {
      const auto& history = isVelocity ? velocity : pressure;
      for (size_t i = 1; i < history.getSize(); ++i)
      {
        const auto& coarse = history.getSample(i - 1).error;
        const auto& fine = history.getSample(i).error;
        SCOPED_TRACE(::testing::Message()
          << "velocity=" << isVelocity << " degree=" << i + 2 << " L2 " << coarse.getL2()
          << " -> " << fine.getL2() << " H1 " << coarse.getH1Seminorm() << " -> "
          << fine.getH1Seminorm());
        ASSERT_TRUE(coarse.isFinite());
        ASSERT_TRUE(fine.isFinite());
        ASSERT_GT(fine.getL2(), 0);
        ASSERT_GT(fine.getH1Seminorm(), 0);
        ASSERT_GT(coarse.getL2(), fine.getL2());
        ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
        const auto decay = history.getExponentialRates(i);
        EXPECT_GT(decay.getL2(), 0.1);
        EXPECT_GT(decay.getH1Seminorm(), 0.1);
      }
    }
  }

  TEST_P(StokesPTest, EveryPairReproducesPolynomialPatches)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const StokesProblem quadratic(
      mesh, StokesData(mesh.getDimension(), StokesData::Field::Quadratic));
    const StokesProblem cubic(
      mesh, StokesData(mesh.getDimension(), StokesData::Field::Cubic));
    const StokesProblem quartic(
      mesh, StokesData(mesh.getDimension(), StokesData::Field::Quartic));
    for (const auto& errors :
      {quadratic.solve<2>(), cubic.solve<3>(), quartic.solve<4>()})
    {
      EXPECT_LT(errors.velocity.getL2(), 1e-9);
      EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
      EXPECT_LT(errors.pressure.getL2(), 1e-9);
      EXPECT_LT(errors.pressure.getH1Seminorm(), 1e-9);
      EXPECT_LT(errors.divergence, 1e-9);
    }
  }

  TEST_P(StokesPTest, EveryPairRejectsIncorrectViscosity)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const StokesProblem problem(
      mesh, StokesData(mesh.getDimension(), StokesData::Field::Quadratic));
    for (const auto& errors :
      {problem.solve<2>(2), problem.solve<3>(2), problem.solve<4>(2)})
    {
      // The wrong viscosity is absorbed by pressure, not necessarily velocity.
      EXPECT_LT(errors.velocity.getL2(), 1e-9);
      EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
      EXPECT_GT(errors.pressure.getL2(), 0.1);
      EXPECT_GT(errors.pressure.getH1Seminorm(), 1);
    }
  }

  TEST_P(StokesPTest, P4QuadratureSensitivity)
  {
    auto mesh = UniformGrid(GetParam()).makeMesh(3);
    const StokesData data(mesh.getDimension(), StokesData::Field::Smooth);
    const auto baseline = StokesProblem(mesh, data).solve<4>();
    const auto refined = StokesProblem(mesh, data, 18).solve<4>();
    for (const auto& pair : {std::pair{baseline.velocity, refined.velocity},
           std::pair{baseline.pressure, refined.pressure}})
    {
      ASSERT_GT(pair.first.getL2(), 0);
      ASSERT_GT(pair.first.getH1Seminorm(), 0);
      EXPECT_LT(std::abs(pair.second.getL2() / pair.first.getL2() - 1), 1e-6);
      EXPECT_LT(
        std::abs(pair.second.getH1Seminorm() / pair.first.getH1Seminorm() - 1), 1e-6);
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, StokesPTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
