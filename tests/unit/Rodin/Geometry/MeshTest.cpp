/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Math/Vector.h"
#include <fstream>
#include <gtest/gtest.h>

#include <Rodin/IO.h>
#include <Rodin/IO/MEDIT.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Geometry_Mesh, 2D_Square_Build)
  {
    constexpr const size_t sdim = 2;
    // constexpr const size_t mdim = sdim;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 2);
  }

  TEST(Rodin_Geometry_Mesh, 2D_Square_Boundary)
  {
    constexpr const size_t sdim = 2;
    constexpr const size_t mdim = sdim;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 2);

    mesh.getConnectivity().compute(mdim - 1, mdim);

    size_t count = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      count += 1;
      ASSERT_TRUE(it->isBoundary());
      ASSERT_FALSE(it->isInterface());
    }

    EXPECT_EQ(count, 4);
  }

  TEST(Rodin_Geometry_Mesh, 2D_Square_Interface)
  {
    constexpr const size_t sdim = 2;
    constexpr const size_t mdim = sdim;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 2);

    mesh.getConnectivity().compute(mdim - 1, mdim);

    size_t count = 0;
    for (auto it = mesh.getInterface(); !it.end(); ++it)
    {
      count += 1;
      ASSERT_TRUE(it->isInterface());
      ASSERT_FALSE(it->isBoundary());
    }

    EXPECT_EQ(count, 1);
  }

  TEST(Rodin_Geometry_Mesh_FuzzyTest, 2D_Square_PolytopeTransformation_1)
  {
    constexpr const size_t rdim = 2;
    constexpr const size_t sdim = 2;
    constexpr const size_t meshDim = 2;

    Mesh mesh =
      Mesh<Rodin::Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    {
      const auto& trans = mesh.getPolytopeTransformation(meshDim, 0);
      Math::SpatialPoint rc(rdim);
      rc[0] = 0.5;
      rc[1] = 0.5;

      Math::SpatialPoint pc(sdim);
      pc[0] = 0.5;
      pc[1] = 0.5;

      Math::SpatialPoint res;
      trans.transform(res, rc);
      Math::SpatialPoint inv;
      trans.inverse(inv, res);
      EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
    }

    {
      const auto& trans = mesh.getPolytopeTransformation(meshDim, 0);

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 0.0;
        rc[1] = 0.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = 0.0;
        pc[1] = 0.0;

        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 1.0;
        rc[1] = 0.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = 1.0;
        pc[1] = 0.0;

        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 0.0;
        rc[1] = 1.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = 0.0;
        pc[1] = 1.0;

        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }

      {
        Math::SpatialPoint r(rdim);
        r[0] = 1.0 / 3.0;
        r[1] = 1.0 / 3.0;

        Math::SpatialPoint p(sdim);
        p[0] = (1.0 / 3.0);
        p[1] = (1.0 / 3.0);

        Math::SpatialPoint res;
        trans.transform(res, r);
        EXPECT_NEAR((res - p).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }

    {
      const auto& trans = mesh.getPolytopeTransformation(meshDim, 1);

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 0.0;
        rc[1] = 0.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = 1.0;
        pc[1] = 0.0;

        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 1.0;
        rc[1] = 0.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = 1.0;
        pc[1] = 1.0;

        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 0.0;
        rc[1] = 1.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = 0.0;
        pc[1] = 1.0;

        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }

      {
        Math::SpatialPoint rc(rdim);
        rc[0] = 1.0 / 3.0;
        rc[1] = 1.0 / 3.0;

        Math::SpatialPoint pc(sdim);
        pc[0] = (2.0 / 3.0);
        pc[1] = (2.0 / 3.0);


        Math::SpatialPoint res;
        trans.transform(res, rc);
        Math::SpatialPoint inv;
        trans.inverse(inv, res);
        EXPECT_NEAR((res - pc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
        EXPECT_NEAR((inv - rc).norm(), 0.0, RODIN_FUZZY_CONSTANT);
      }
    }
  }


  // Additional comprehensive tests

  TEST(Rodin_Geometry_Mesh, MultipleTriangles)
  {
    constexpr const size_t sdim = 2;

    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(sdim)
      .nodes(5)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .vertex({0.5, 0.5})
      .polytope(Polytope::Type::Triangle, {0, 1, 4})
      .polytope(Polytope::Type::Triangle, {1, 3, 4})
      .polytope(Polytope::Type::Triangle, {3, 2, 4})
      .polytope(Polytope::Type::Triangle, {2, 0, 4})
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 5);
    EXPECT_EQ(mesh.getCellCount(), 4);
  }

  TEST(Rodin_Geometry_Mesh, 3D_Tetrahedron)
  {
    constexpr const size_t sdim = 3;

    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0, 0})
      .vertex({1, 0, 0})
      .vertex({0, 1, 0})
      .vertex({0, 0, 1})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 1);
  }
}
