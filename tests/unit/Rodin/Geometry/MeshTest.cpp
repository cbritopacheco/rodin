/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Math/Vector.h"
#include <fstream>
#include <array>
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include <Rodin/IO.h>
#include <Rodin/IO/MEDIT.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Closed surfaces and single cells report absent face selections without changing topology.
  TEST(Rodin_Geometry_Mesh, EmptyFaceSelectionsReportMemberContext)
  {
    for (bool closed : {false, true})
    {
      SCOPED_TRACE(closed);
      auto mesh = closed ? LocalMesh::Box(Polytope::Type::Triangle, {2, 2, 2})
                         : LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {2, 2});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      const auto vertices = mesh.getVertexCount();
      const auto faces = mesh.getFaceCount();
      const auto cells = mesh.getCellCount();
      ASSERT_GT(faces, 0u);
      for (Index face = 0; face < faces; ++face)
      {
        EXPECT_EQ(
          mesh.getConnectivity().getIncidence({1, 2}, face).size(), closed ? 2u : 1u);
        EXPECT_EQ(mesh.isBoundary(face), !closed);
        EXPECT_EQ(mesh.isInterface(face), closed);
      }
      try
      {
        if (closed)
          mesh.getBoundary();
        else
          mesh.getInterface();
        FAIL() << "An absent face selection was accepted";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find(closed ? "In member function getBoundary"
                                      : "In member function getInterface"),
          std::string::npos);
        EXPECT_NE(
          message.find(closed ? "empty boundary" : "empty interface"), std::string::npos);
      }
      EXPECT_EQ(mesh.getVertexCount(), vertices);
      EXPECT_EQ(mesh.getFaceCount(), faces);
      EXPECT_EQ(mesh.getCellCount(), cells);
      size_t selected = 0;
      for (auto face = closed ? mesh.getInterface() : mesh.getBoundary(); face; ++face)
        ++selected;
      EXPECT_EQ(selected, faces);
    }
  }

  /// @brief Public mesh file operations replace existing geometry and preserve physical cells and attributes.
  TEST(Rodin_Geometry_Mesh, PublicFileOperationsRoundTripAndReplaceMesh)
  {
    const auto directory = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-public-mesh-%%%%-%%%%-%%%%");
    ASSERT_TRUE(boost::filesystem::create_directories(directory));
    auto original = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    original.scale(0.5);
    original.setAttribute({1, 0}, 17);
    original.setAttribute({1, 1}, 23);
    const FlatSet<Real> expectedVertices{0, 0.5, 1};
    const FlatSet<std::array<Real, 2>> expectedCells{
      std::array<Real, 2>{0, 0.5}, std::array<Real, 2>{0.5, 1}};
    for (const auto format : {FileFormat::MFEM, FileFormat::MEDIT, FileFormat::HDF5})
    {
      SCOPED_TRACE(static_cast<int>(format));
      const auto filename = directory /
        (std::to_string(static_cast<int>(format)) +
          (format == FileFormat::HDF5 ? ".h5" : ".mesh"));
      original.save(filename, format);
      auto loaded = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, {2, 2});
      EXPECT_EQ(&loaded.load(filename, format), &loaded);
      EXPECT_EQ(loaded.getDimension(), 1u);
      EXPECT_EQ(loaded.getSpaceDimension(), 1u);
      ASSERT_EQ(loaded.getVertexCount(), 3u);
      ASSERT_EQ(loaded.getCellCount(), 2u);
      FlatSet<Real> actualVertices;
      for (Index vertex = 0; vertex < loaded.getVertexCount(); ++vertex)
        actualVertices.insert(loaded.getVertices()(0, vertex));
      EXPECT_EQ(actualVertices, expectedVertices);
      FlatSet<std::array<Real, 2>> actualCells;
      for (auto cell = loaded.getCell(); cell; ++cell)
      {
        EXPECT_EQ(cell->getGeometry(), Polytope::Type::Segment);
        const auto& vertices = cell->getVertices();
        ASSERT_EQ(vertices.size(), 2u);
        std::array<Real, 2> coordinates{
          loaded.getVertices()(0, vertices[0]), loaded.getVertices()(0, vertices[1])};
        std::sort(coordinates.begin(), coordinates.end());
        actualCells.insert(coordinates);
        EXPECT_EQ(cell->getAttribute(), coordinates[0] == 0 ? 17 : 23);
      }
      EXPECT_EQ(actualCells, expectedCells);
    }
    EXPECT_GT(boost::filesystem::remove_all(directory), 0u);
    EXPECT_FALSE(boost::filesystem::exists(directory));
  }

  /// @brief Unsupported mesh formats report member context before touching files or mesh state.
  TEST(Rodin_Geometry_Mesh, UnsupportedFileFormatsPreserveMeshAndFilesystem)
  {
    const auto directory = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-unsupported-mesh-%%%%-%%%%-%%%%");
    const auto filename = directory / "mesh";
    ASSERT_FALSE(boost::filesystem::exists(directory));
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.setAttribute({2, 0}, 17);
    const auto original = mesh.getVertices();
    for (int value : {-1, 255})
    {
      SCOPED_TRACE(value);
      const auto format = static_cast<FileFormat>(value);
      for (bool loading : {true, false})
      {
        SCOPED_TRACE(loading);
        try
        {
          if (loading)
            mesh.load(filename, format);
          else
            mesh.save(filename, format);
          FAIL() << "An unsupported mesh format was accepted";
        }
        catch (const Alert::Exception& error)
        {
          const std::string message = error.what();
          EXPECT_NE(
            message.find(loading ? "In member function load" : "In member function save"),
            std::string::npos);
          EXPECT_NE(message.find("format unsupported"), std::string::npos);
          EXPECT_NE(message.find(std::to_string(value)), std::string::npos);
        }
        EXPECT_EQ(mesh.getDimension(), 2u);
        EXPECT_EQ(mesh.getSpaceDimension(), 2u);
        EXPECT_EQ(mesh.getCellCount(), 2u);
        EXPECT_EQ(mesh.getAttribute(2, 0), 17);
        ASSERT_EQ(mesh.getVertexCount(), original.getCount());
        for (size_t vertex = 0; vertex < original.getCount(); ++vertex)
        {
          for (size_t coordinate = 0; coordinate < original.getDimension(); ++coordinate)
          {
            EXPECT_EQ(
              mesh.getVertices()(coordinate, vertex), original(coordinate, vertex));
          }
        }
        EXPECT_FALSE(boost::filesystem::exists(directory));
      }
    }
  }

  /// @brief Failed file opens preserve mesh state and create no output.
  TEST(Rodin_Geometry_Mesh, FileOpenFailuresPreserveMeshAndFilesystem)
  {
    const auto missingDirectory = boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-missing-mesh-%%%%-%%%%-%%%%");
    const auto missingFile = missingDirectory / "mesh";
    ASSERT_FALSE(boost::filesystem::exists(missingDirectory));
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.setAttribute({2, 0}, 17);
    const auto vertices = mesh.getVertices();
    for (const auto format : {FileFormat::MFEM, FileFormat::MEDIT})
    {
      SCOPED_TRACE(static_cast<int>(format));
      EXPECT_THROW(mesh.load(missingFile, format), Alert::Exception);
      EXPECT_THROW(mesh.save(missingFile, format), Alert::Exception);
      EXPECT_EQ(mesh.getCellCount(), 2u);
      EXPECT_EQ(mesh.getAttribute(2, 0), 17);
      ASSERT_EQ(mesh.getVertexCount(), vertices.getCount());
      for (size_t vertex = 0; vertex < vertices.getCount(); ++vertex)
      {
        for (size_t component = 0; component < vertices.getDimension(); ++component)
          EXPECT_EQ(mesh.getVertices()(component, vertex), vertices(component, vertex));
      }
      EXPECT_FALSE(boost::filesystem::exists(missingDirectory));
    }
  }

  /// @brief Boundary trace remaps selected exterior faces and preserves the interface.
  TEST(Rodin_Geometry_Mesh, BoundaryAttributeTraceHonorsSelection)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(2)
                  .nodes(4)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .vertex({1, 1})
                  .polytope(Polytope::Type::Triangle, {0, 1, 2})
                  .polytope(Polytope::Type::Triangle, {1, 3, 2})
                  .finalize();
    mesh.getConnectivity().compute(1, 2);
    mesh.setAttribute({2, 0}, 17);
    mesh.setAttribute({2, 1}, 23);
    for (auto face = mesh.getFace(); face; ++face)
    {
      const Point center(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const bool lowerOrLeft = center.x() == 0 || center.y() == 0;
      mesh.setAttribute(
        {1, face->getIndex()}, face->isBoundary() ? (lowerOrLeft ? 31 : 37) : 41);
    }
    const Map<Attribute, Attribute> mapping{{17, 71}, {23, 73}};
    EXPECT_EQ(&mesh.trace(mapping, FlatSet<Attribute>{31}), &mesh);
    for (auto face = mesh.getFace(); face; ++face)
    {
      const Point center(*face, Polytope::Traits(face->getGeometry()).getCentroid());
      const bool lowerOrLeft = center.x() == 0 || center.y() == 0;
      EXPECT_EQ(face->getAttribute(), face->isBoundary() ? (lowerOrLeft ? 71 : 37) : 41);
    }
    for (size_t repetition = 0; repetition < 2; ++repetition)
    {
      mesh.trace(mapping);
      mesh.trace(Map<Attribute, Attribute>{{999, 101}});
      for (auto face = mesh.getFace(); face; ++face)
      {
        const Point center(*face, Polytope::Traits(face->getGeometry()).getCentroid());
        const bool lowerOrLeft = center.x() == 0 || center.y() == 0;
        EXPECT_EQ(
          face->getAttribute(), face->isBoundary() ? (lowerOrLeft ? 71 : 73) : 41);
      }
    }
  }

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

    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    Mesh mesh = Mesh<Rodin::Context::Local>::Builder()
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

    Mesh mesh = Mesh<Context::Local>::Builder()
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

    Mesh mesh = Mesh<Context::Local>::Builder()
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
