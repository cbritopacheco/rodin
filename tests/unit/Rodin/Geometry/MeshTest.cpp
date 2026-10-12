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
  /// @brief Attribute-restricted measures recover half-unit-box volumes and boundaries.
  TEST(Rodin_Geometry_Mesh, AttributeRestrictedUnitBoxMeasures)
  {
    for (const auto geometry : {Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
           Polytope::Type::Tetrahedron, Polytope::Type::Hexahedron,
           Polytope::Type::Pyramid, Polytope::Type::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      const size_t dimension = Polytope::Traits(geometry).getDimension();
      const Array<size_t> shape = Array<size_t>::Constant(dimension, 3);
      auto mesh = LocalMesh::UniformGrid(geometry, shape);
      mesh.scale(0.5);
      mesh.getConnectivity().compute(dimension - 1, dimension);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const Point center(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
        mesh.setAttribute({dimension, cell->getIndex()}, center.x() < 0.5 ? 17 : 23);
      }
      for (auto face = mesh.getBoundary(); face; ++face)
      {
        const Point center(*face, Polytope::Traits(face->getGeometry()).getCentroid());
        mesh.setAttribute({dimension - 1, face->getIndex()}, center.x() < 0.5 ? 31 : 37);
      }

      for (const Attribute attribute : {17, 23})
      {
        EXPECT_NEAR(mesh.getMeasure(dimension, attribute), 0.5, 1e-13);
        EXPECT_NEAR(
          mesh.getMeasure(dimension, FlatSet<Attribute>{attribute}), 0.5, 1e-13);
        if (dimension == 2)
        {
          EXPECT_NEAR(mesh.getArea(attribute), 0.5, 1e-13);
          EXPECT_NEAR(mesh.getArea(FlatSet<Attribute>{attribute}), 0.5, 1e-13);
        }
        else
        {
          EXPECT_NEAR(mesh.getVolume(attribute), 0.5, 1e-13);
          EXPECT_NEAR(mesh.getVolume(FlatSet<Attribute>{attribute}), 0.5, 1e-13);
        }
      }
      EXPECT_NEAR(mesh.getMeasure(dimension, FlatSet<Attribute>{17, 23, 999}), 1, 1e-13);
      EXPECT_EQ(mesh.getMeasure(dimension, Attribute(999)), 0);
      EXPECT_EQ(mesh.getMeasure(dimension, FlatSet<Attribute>{}), 0);
      EXPECT_EQ(mesh.getMeasure(dimension, FlatSet<Attribute>{999}), 0);
      // A half of the unit square/cube's exterior has measure 2/3.
      for (const Attribute attribute : {31, 37})
      {
        EXPECT_NEAR(mesh.getPerimeter(attribute), dimension, 1e-13);
        EXPECT_NEAR(mesh.getPerimeter(FlatSet<Attribute>{attribute}), dimension, 1e-13);
      }
      EXPECT_NEAR(
        mesh.getPerimeter(FlatSet<Attribute>{31, 37, 999}), 2 * dimension, 1e-13);
      EXPECT_EQ(mesh.getPerimeter(Attribute(999)), 0);
      EXPECT_EQ(mesh.getPerimeter(FlatSet<Attribute>{}), 0);
    }
  }

  /// @brief Verifies 2 D square build for geometry mesh by checking exact expected values.
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

  /// @brief Verifies 2 D square boundary for geometry mesh by checking exact expected values, true predicates, false predicates.
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

  /// @brief Verifies 2 D square interface for geometry mesh by checking exact expected values, true predicates, false predicates.
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

  /// @brief Verifies 2 D square polytope transformation 1 for geometry mesh fuzzy test by checking tolerance-based numerical results.
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

  /// @brief Verifies multiple triangles for geometry mesh by checking exact expected values.
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

  /// @brief Verifies 3 D tetrahedron for geometry mesh by checking exact expected values.
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
