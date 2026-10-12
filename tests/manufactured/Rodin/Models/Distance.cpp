/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Exact interval solutions for distance-model wrappers. */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Distance/Eikonal.h"
#include "Rodin/Distance/Poisson.h"
#include "Rodin/Distance/SignedPoisson.h"
#include "Rodin/Distance/Rvachev.h"
#include "Rodin/Distance/SpaldingTucker.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured
{
  /**
   * @brief Builds a unit interval with the incidences required by distance solvers.
   * @param count Number of vertices.
   * @returns Uniform interval mesh.
   */
  LocalMesh distanceInterval(size_t count)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {count});
    mesh.scale(Real(1) / (count - 1));
    mesh.getConnectivity().compute(0, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(0, 0);
    return mesh;
  }

  /// @brief Exact marching distances and signs survive repeated solves and shared vertices.
  TEST(Rodin_Manufactured_Distance, EikonalSignedInterval)
  {
    auto mesh = distanceInterval(5);
    mesh.setAttribute({0, 2}, 9);
    for (Index cell = 0; cell < mesh.getCellCount(); ++cell)
      mesh.setAttribute({1, cell}, cell < 2 ? 7 : 8);
    P1 fes(mesh);
    GridFunction solution(fes);
    Distance::Eikonal model(solution);
    model.setInterface(Attribute(9), Attribute(99))
      .setInterior(Attribute(7), Attribute(70));
    EXPECT_EQ(model.getInterface(), (FlatSet<Attribute>{9, 99}));
    EXPECT_EQ(model.getInterior(), (FlatSet<Attribute>{7, 70}));
    for (size_t repeat = 0; repeat < 2; ++repeat)
    {
      model.solve();
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        const Real x = mesh.getVertexCoordinates(vertex)[0];
        EXPECT_NEAR(solution[vertex], std::abs(x - Real(0.5)), 1e-12);
      }
      model.sign();
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        const Real x = mesh.getVertexCoordinates(vertex)[0];
        EXPECT_NEAR(solution[vertex], x - Real(0.5), 1e-12);
      }
    }
    const FlatSet<Attribute> interface{9}, interior{7};
    model.setInterface(interface).setInterior(interior);
    EXPECT_EQ(model.getInterface(), interface);
    EXPECT_EQ(model.getInterior(), interior);
  }

  /// @brief A missing interface cannot silently produce an infinite distance field.
  TEST(Rodin_Manufactured_Distance, EikonalMissingInterface)
  {
    auto mesh = distanceInterval(5);
    P1 fes(mesh);
    GridFunction solution(fes);
    Distance::Eikonal model(solution);
    model.setInterface(Attribute(999));
    EXPECT_THROW(model.solve(), Alert::Exception);
  }

  /// @brief P1 Poisson nodal values are exact and midpoint error decreases quadratically.
  TEST(Rodin_Manufactured_Distance, PoissonP1Refinement)
  {
    Real previousError = 0;
    for (size_t count : {5u, 9u})
    {
      auto mesh = distanceInterval(count);
      P1 fes(mesh);
      const auto solution = Distance::Poisson()(fes);
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        const Real x = mesh.getVertexCoordinates(vertex)[0];
        EXPECT_NEAR(solution[vertex], x * (1 - x) / 2, 1e-12);
      }
      Real error = 0;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        const Real x = point.getCoordinates()[0];
        error = std::max(error, std::abs(solution(point) - x * (1 - x) / 2));
      }
      const Real h = Real(1) / (count - 1);
      EXPECT_NEAR(error, h * h / 8, 1e-12);
      if (previousError > 0)
      {
        EXPECT_NEAR(previousError / error, 4, 1e-10);
      }
      previousError = error;
    }
  }

  /// @brief A quadratic H1 space reproduces both signed and unsigned Poisson fields exactly.
  TEST(Rodin_Manufactured_Distance, PoissonH1Quadratic)
  {
    auto mesh = distanceInterval(5);
    H1<2, Real> fes(std::integral_constant<size_t, 2>{}, mesh);
    const auto positive = Distance::Poisson()(fes);
    const auto negative = Distance::SignedPoisson()(fes);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      for (Real reference : {Real(0.2), Real(0.5), Real(0.8)})
      {
        Point point(*cell, Math::SpatialVector<Real>{reference});
        const Real x = point.getCoordinates()[0];
        const Real exact = x * (1 - x) / 2;
        EXPECT_NEAR(positive(point), exact, 1e-12);
        EXPECT_NEAR(negative(point), -exact, 1e-12);
      }
    }
  }

  /// @brief The signed full-domain overload solves the negative unit-source problem.
  TEST(Rodin_Manufactured_Distance, SignedPoissonFullDomain)
  {
    auto mesh = distanceInterval(9);
    P1 fes(mesh);
    const auto solution = Distance::SignedPoisson()(fes);
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const Real x = mesh.getVertexCoordinates(vertex)[0];
      EXPECT_NEAR(solution[vertex], -x * (1 - x) / 2, 1e-12);
    }
  }

  /// @brief Region forcing and an interior zero interface give exact piecewise quadratics.
  TEST(Rodin_Manufactured_Distance, SignedPoissonRegionOverloads)
  {
    auto mesh = distanceInterval(9);
    mesh.setAttribute({0, 4}, 9);
    for (Index cell = 0; cell < mesh.getCellCount(); ++cell)
      mesh.setAttribute({1, cell}, cell < 4 ? 7 : 8);
    P1 fes(mesh);
    const auto single = Distance::SignedPoisson()(9, 7, fes);
    const auto sets =
      Distance::SignedPoisson()(FlatSet<Attribute>{9}, FlatSet<Attribute>{7}, fes);
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const Real x = mesh.getVertexCoordinates(vertex)[0];
      const Real exact =
        x <= Real(0.5) ? (x * x - Real(0.25)) / 2 : (Real(0.25) - (1 - x) * (1 - x)) / 2;
      EXPECT_NEAR(single[vertex], exact, 1e-12);
      EXPECT_NEAR(sets[vertex], exact, 1e-12);
    }
  }

  /// @brief Affine level sets have analytically known normalization and unchanged inputs.
  TEST(Rodin_Manufactured_Distance, AffineNormalizations)
  {
    auto mesh = distanceInterval(9);
    P1 fes(mesh);
    GridFunction input(fes);
    input = [](const Point& point) { return 2 * point.getCoordinates()[0] - 1; };
    const auto rvachev = Distance::Rvachev()(input);
    const auto spalding = Distance::SpaldingTucker()(input);
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const Real value = 2 * mesh.getVertexCoordinates(vertex)[0] - 1;
      EXPECT_NEAR(input[vertex], value, 1e-12);
      EXPECT_NEAR(rvachev[vertex], value / std::sqrt(value * value + 4), 1e-12);
      EXPECT_NEAR(
        spalding[vertex], 2 * value / (2 + std::sqrt(4 + 2 * std::abs(value))), 1e-12);
    }
  }
}
