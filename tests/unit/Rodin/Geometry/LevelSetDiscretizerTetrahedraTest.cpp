/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTetrahedra.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  static LocalMesh makeUnitTetrahedron()
  {
    return LocalMesh::Builder()
      .initialize(3)
      .nodes(4)
      .vertex({0, 0, 0})
      .vertex({1, 0, 0})
      .vertex({0, 1, 0})
      .vertex({0, 0, 1})
      .tetrahedron(IndexArray{{0, 1, 2, 3}})
      .finalize();
  }

  static long double quality(const LocalMesh& mesh, const Polytope& cell)
  {
    const auto& v = cell.getVertices();
    const auto& A = mesh.getVertexCoordinates(v(0));
    const auto& B = mesh.getVertexCoordinates(v(1));
    const auto& C = mesh.getVertexCoordinates(v(2));
    const auto& D = mesh.getVertexCoordinates(v(3));

    const auto AB = B - A;
    const auto AC = C - A;
    const auto AD = D - A;

    const long double vol6 =
      std::abs((long double) AB.dot(AC.cross(AD)));

    long double hmax2 = 0.0L;
    auto update = [&](const auto& X, const auto& Y)
    {
      const auto d = Y - X;
      hmax2 = std::max(hmax2, (long double) d.squaredNorm());
    };

    update(A, B); update(A, C); update(A, D);
    update(B, C); update(B, D); update(C, D);

    if (!(hmax2 > 0.0L))
      return 0.0L;

    const long double hmax = std::sqrt(hmax2);
    return vol6 / (hmax * hmax * hmax);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTetrahedra, SnapsNearVertexSplit22)
  {
    for (const Real eps : { 1e-12, 1e-9, 1e-6, 1e-3 })
    {
      LocalMesh mesh = makeUnitTetrahedron();
      auto& conn = mesh.getConnectivity();
      conn.compute(3, 2);
      conn.compute(3, 1);
      conn.compute(2, 1);
      conn.compute(1, 0);

      P1 vh(mesh);
      GridFunction phi(vh);
      phi[0] = -eps;
      phi[1] = -1.0;
      phi[2] =  1.0;
      phi[3] =  1.0;

      LevelSetDiscretizerTetrahedra lsd(phi);
      LocalMesh split = lsd.discretize();

      ASSERT_GT(split.getCellCount(), 0);

      long double qmin = std::numeric_limits<long double>::infinity();
      for (auto it = split.getCell(); !it.end(); ++it)
      {
        const long double q = quality(split, *it);
        EXPECT_TRUE(std::isfinite((double) q));
        qmin = std::min(qmin, q);
      }

      EXPECT_GE(qmin, lsd.getMinimumQuality());
    }
  }
}
