/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Isoparametric P2 geometry and Poisson convergence validation.
 *
 * The global map @f$\Phi(\xi)_i=\xi_i@f$ for @f$i<d-1@f$ and
 * @f$\Phi(\xi)_{d-1}=\xi_{d-1}+a\xi_0^2@f$ is a regular quadratic map.
 * It is reproduced exactly by a P2 geometry element, while its P1 nodal
 * interpolant is only piecewise affine.  The Poisson problem is posed in
 * physical coordinates with @f$u(x)=\prod_i\sin(\pi x_i)@f$ and
 * @f$f=d\pi^2u@f$, so both its load and trace are known independently of the
 * mesh parametrization.
 */

#include <cstdint>
#include <initializer_list>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "Convergence.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/SparseLU.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::Poisson
{
  constexpr Real Warp = Real(0.1);

  Math::SpatialPoint mapToPhysical(Math::SpatialPoint point)
  {
    assert(point.size() >= 1);
    point(point.size() - 1) += Warp * point(0) * point(0);
    return point;
  }

  Math::SpatialPoint referencePosition(
    const Polytope& polytope,
    const std::vector<Math::SpatialPoint>& vertices,
    const Math::SpatialPoint& rc)
  {
    RealP1Element affine(polytope.getGeometry());
    Math::SpatialPoint point(vertices.front().size());
    point.setZero();
    const auto polytopeVertices = polytope.getVertices();
    for (size_t local = 0; local < affine.getCount(); ++local)
      point += vertices.at(polytopeVertices[local])
        * affine.getBasis(local)(rc);
    return point;
  }

  /**
   * @brief Curves every positive-dimensional entity consistently.
   *
   * The original vertex positions are retained while the mesh vertices are
   * moved. Control points are then evaluated from that original affine mesh;
   * this avoids interpolating the already-warped P1 geometry a second time.
   */
  template <size_t GeometryOrder>
  void installGeometry(LocalMesh& mesh)
  {
    const size_t dim = mesh.getSpaceDimension();
    std::vector<Math::SpatialPoint> vertices;
    vertices.reserve(mesh.getVertexCount());
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      vertices.push_back(mesh.getVertexCoordinates(vertex));

    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      mesh.setVertexCoordinates(vertex, mapToPhysical(vertices.at(vertex)));

    if constexpr (GeometryOrder == 1)
    {
      return;
    }
    else
    {
      for (size_t entityDimension = 1; entityDimension <= dim;
        ++entityDimension)
      {
        for (Index index = 0;
          index < mesh.getPolytopeCount(entityDimension); ++index)
        {
          const auto polytope = mesh.getPolytope(entityDimension, index);
          RealH1Element<GeometryOrder> geometryElement(
            polytope->getGeometry());
          PointCloud nodes(dim, geometryElement.getCount());
          for (size_t local = 0; local < geometryElement.getCount(); ++local)
          {
            const auto point = mapToPhysical(referencePosition(
              *polytope, vertices, geometryElement.getNode(local)));
            for (size_t coordinate = 0; coordinate < dim; ++coordinate)
              nodes(coordinate, local) = point(coordinate);
          }
          mesh.setPolytopeTransformation({entityDimension, index},
            new ParametricTransformation<RealH1Element<GeometryOrder>>(
              std::move(nodes), std::move(geometryElement)));
        }
      }
    }
  }

  ErrorNorms solveP2(const LocalMesh& mesh)
  {
    const size_t dim = mesh.getSpaceDimension();
    const Real pi = Math::Constants::pi();
    const RealFunction exact([dim, pi](const Point& p)
      {
        Real value = 1;
        for (size_t i = 0; i < dim; ++i)
          value *= std::sin(pi * p(i));
        return value;
      });
    const RealFunction forcing([dim, pi](const Point& p)
      {
        Real value = Real(dim) * pi * pi;
        for (size_t i = 0; i < dim; ++i)
          value *= std::sin(pi * p(i));
        return value;
      });
    const VectorFunction gradient(dim, [dim, pi](const Point& p)
      {
        Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
        for (size_t i = 0; i < dim; ++i)
        {
          value(i) = pi * std::cos(pi * p(i));
          for (size_t j = 0; j < dim; ++j)
            if (j != i)
              value(i) *= std::sin(pi * p(j));
        }
        return value;
      });

    H1 space(std::integral_constant<size_t, 2>{}, mesh);
    TrialFunction u(space);
    TestFunction v(space);
    auto stiffness = Integral(Grad(u), Grad(v));
    auto load = Integral(forcing, v);
    stiffness.setOrder(16);
    load.setOrder(16);
    Problem problem(u, v);
    problem = stiffness - load + DirichletBC(u, exact);

    SparseLU solver(problem);
    solver.solve();
    return ErrorNorm::compute(mesh, u.getSolution(), exact, gradient, 16);
  }

  void expectOptimalP2Rates(const ErrorHistory& history)
  {
    ASSERT_EQ(history.getSize(), 3);
    for (size_t i = 1; i < history.getSize(); ++i)
    {
      const auto& coarse = history.getSample(i - 1).error;
      const auto& fine = history.getSample(i).error;
      ASSERT_TRUE(coarse.isFinite());
      ASSERT_TRUE(fine.isFinite());
      ASSERT_GT(coarse.getL2(), fine.getL2());
      ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
      const auto rates = history.getAlgebraicRates(i);
      SCOPED_TRACE(::testing::Message()
        << "L2 " << coarse.getL2() << " -> " << fine.getL2()
        << ", rate " << rates.getL2() << "; H1 "
        << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
        << ", rate " << rates.getH1Seminorm());
      EXPECT_GT(rates.getL2(), 2.35);
      EXPECT_LT(rates.getL2(), 3.65);
      EXPECT_GT(rates.getH1Seminorm(), 1.45);
      EXPECT_LT(rates.getH1Seminorm(), 2.55);
    }
  }

  class PoissonIsoparametricTest
    : public ::testing::TestWithParam<Polytope::Type>
  {};

  /**
   * @brief P2 geometry represents the prescribed quadratic curved map exactly.
   *
   * Sampling an interior reference point of every cell catches wrong geometry
   * element selection, control-point ordering, and entity transformation
   * installation. This is a geometric patch test, independent of the PDE.
   */
  TEST_P(PoissonIsoparametricTest, P2GeometryReproducesQuadraticMap)
  {
    const auto geometry = GetParam();
    UniformGrid grid(geometry);
    auto mesh = grid.makeMesh(grid.getDimension() == 3 ? 3 : 5);
    std::vector<Math::SpatialPoint> vertices;
    vertices.reserve(mesh.getVertexCount());
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      vertices.push_back(mesh.getVertexCoordinates(vertex));

    installGeometry<2>(mesh);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& transformation = cell->getTransformation();
      RealH1Element<2> element(cell->getGeometry());
      Math::SpatialPoint rc(Polytope::Traits(cell->getGeometry()).getDimension());
      rc.setZero();
      for (size_t local = 0; local < element.getCount(); ++local)
        rc += element.getNode(local);
      rc /= Real(element.getCount());
      Math::SpatialPoint physical;
      transformation.transform(physical, rc);
      const auto expected = mapToPhysical(referencePosition(
        *cell, vertices, rc));
      EXPECT_LT((physical - expected).norm(), 1e-11);
    }
  }

  /**
   * @brief P2 Poisson fields retain their optimal h-rates on P2 curved maps.
   *
   * Isoparametric interpolation represents the quadratic geometry exactly.
   * Standard mapped-element estimates thus retain the P2 rates
   * @f$O(h^3)@f$ in L2 and @f$O(h^2)@f$ in the H1 seminorm.
   */
  TEST_P(PoissonIsoparametricTest, P2PoissonHasOptimalRatesOnCurvedGeometry)
  {
    const auto geometry = GetParam();
    UniformGrid grid(geometry);
    const std::initializer_list<size_t> levels =
      grid.getDimension() == 3 ? std::initializer_list<size_t>{3, 4, 5}
                               : std::initializer_list<size_t>{5, 9, 17};
    ErrorHistory history;
    for (const size_t pointsPerAxis : levels)
    {
      auto mesh = grid.makeMesh(pointsPerAxis);
      installGeometry<2>(mesh);
      history.append(Real(1) / Real(pointsPerAxis - 1), solveP2(mesh));
    }
    expectOptimalP2Rates(history);
  }

  INSTANTIATE_TEST_SUITE_P(AllUniformGridGeometries,
    PoissonIsoparametricTest,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info)
    {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
