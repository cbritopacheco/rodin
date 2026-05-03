/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <boost/filesystem/fstream.hpp>
#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include <Rodin/IO.h>
#include <Rodin/IO/MFEM.h>
#include <Rodin/IO/MEDIT.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin::IO;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  namespace
  {
    Mesh<Context::Local> makePointMesh()
    {
      Mesh<Context::Local>::Builder builder;
      return builder
        .initialize(1)
        .nodes(1)
        .vertex({ 0.0 })
        .attribute({ 0, 0 }, 7)
        .finalize();
    }

    Mesh<Context::Local> makeCoverageMesh(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Point:
          return makePointMesh();
        case Polytope::Type::Segment:
          return LocalMesh::UniformGrid(geometry, { 5 });
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          return LocalMesh::UniformGrid(geometry, { 3, 3 });
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Hexahedron:
        case Polytope::Type::Wedge:
          return LocalMesh::UniformGrid(geometry, { 2, 2, 2 });
      }
      assert(false);
      return {};
    }

    std::string geometryName(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Point:         return "Point";
        case Polytope::Type::Segment:       return "Segment";
        case Polytope::Type::Triangle:      return "Triangle";
        case Polytope::Type::Quadrilateral: return "Quadrilateral";
        case Polytope::Type::Tetrahedron:   return "Tetrahedron";
        case Polytope::Type::Hexahedron:    return "Hexahedron";
        case Polytope::Type::Wedge:         return "Wedge";
      }
      return "Unknown";
    }

    void computeAllIncidences(Mesh<Context::Local>& mesh)
    {
      const size_t D = mesh.getDimension();
      for (size_t d = 0; d <= D; ++d)
      {
        for (size_t dp = 0; dp <= D; ++dp)
          mesh.getConnectivity().compute(d, dp);
      }
    }

    void expectSameMeshShape(
        const Mesh<Context::Local>& actual,
        const Mesh<Context::Local>& expected)
    {
      ASSERT_EQ(actual.getSpaceDimension(), expected.getSpaceDimension());
      ASSERT_EQ(actual.getDimension(), expected.getDimension());
      ASSERT_EQ(actual.getVertexCount(), expected.getVertexCount());

      for (Index v = 0; v < static_cast<Index>(expected.getVertexCount()); ++v)
      {
        const auto ax = actual.getVertexCoordinates(v);
        const auto ex = expected.getVertexCoordinates(v);
        ASSERT_EQ(ax.size(), ex.size()) << "vertex " << v;
        for (Index d = 0; d < ex.size(); ++d)
          EXPECT_DOUBLE_EQ(ax(d), ex(d)) << "vertex " << v << ", coordinate " << d;
      }

      for (size_t d = 0; d <= expected.getDimension(); ++d)
        ASSERT_EQ(actual.getPolytopeCount(d), expected.getPolytopeCount(d)) << "dimension " << d;

      const size_t D = expected.getDimension();
      for (Index c = 0; c < static_cast<Index>(expected.getCellCount()); ++c)
      {
        EXPECT_EQ(actual.getGeometry(D, c), expected.getGeometry(D, c)) << "cell " << c;
        EXPECT_TRUE(Polytope::Key::SymmetricEquality()(
            actual.getConnectivity().getPolytope(D, c),
            expected.getConnectivity().getPolytope(D, c)))
          << "cell " << c;
      }
    }

    template <FileFormat Format>
    void expectMeshStringRoundTrip(Polytope::Type geometry)
    {
      Mesh mesh = makeCoverageMesh(geometry);
      computeAllIncidences(mesh);

      std::stringstream out;
      MeshPrinter<Format, Context::Local> printer(mesh);
      printer.print(out);

      Mesh loaded;
      MeshLoader<Format, Context::Local> loader(loaded);
      loader.load(out);
      computeAllIncidences(loaded);

      expectSameMeshShape(loaded, mesh);
    }
  }

  class MeshFormatCoverage : public ::testing::TestWithParam<Polytope::Type>
  {};

  TEST_P(MeshFormatCoverage, MEDITRoundTrip)
  {
    expectMeshStringRoundTrip<FileFormat::MEDIT>(GetParam());
  }

  TEST_P(MeshFormatCoverage, MFEMRoundTrip)
  {
    expectMeshStringRoundTrip<FileFormat::MFEM>(GetParam());
  }

  INSTANTIATE_TEST_SUITE_P(
      AllSupportedGeometries,
      MeshFormatCoverage,
      ::testing::Values(
        Polytope::Type::Point,
        Polytope::Type::Segment,
        Polytope::Type::Triangle,
        Polytope::Type::Quadrilateral,
        Polytope::Type::Tetrahedron,
        Polytope::Type::Hexahedron,
        Polytope::Type::Wedge),
      [](const ::testing::TestParamInfo<Polytope::Type>& info)
      {
        return geometryName(info.param);
      });

  TEST(Rodin_IO_MeshLoader, SanityTest_MEDIT_2D_Square)
  {
    static constexpr const char* filename = "medit/Square.medit.mesh";
    static constexpr size_t sdim = 2;
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);

    boost::filesystem::ifstream in(meshfile);

    Mesh mesh;
    MeshLoader<FileFormat::MEDIT, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    EXPECT_EQ(mesh.getSpaceDimension(), sdim);

    size_t d = 0;
    EXPECT_EQ(mesh.getPolytopeCount(d), 4);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getVertexCount());
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 1);
    EXPECT_EQ(mesh.getAttribute(d, 2), 1);
    EXPECT_EQ(mesh.getAttribute(d, 3), 1);

    d = 1;
    EXPECT_EQ(mesh.getPolytopeCount(d), 5);
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);

    d = 2;
    EXPECT_EQ(mesh.getPolytopeCount(d), 2);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getCellCount());
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);
  }

  TEST(Rodin_IO_MeshLoader, SanityTest_MFEM_2D_Square)
  {
    static constexpr const char* filename = "mfem/Square.mfem.mesh";
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);

    boost::filesystem::ifstream in(meshfile);

    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    size_t d = 0;
    EXPECT_EQ(mesh.getPolytopeCount(d), 4);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getVertexCount());

    // MFEM format does not support vertex attributes
    EXPECT_EQ(mesh.getAttribute(d, 0), std::nullopt);
    EXPECT_EQ(mesh.getAttribute(d, 1), std::nullopt);
    EXPECT_EQ(mesh.getAttribute(d, 2), std::nullopt);
    EXPECT_EQ(mesh.getAttribute(d, 3), std::nullopt);

    d = 1;
    EXPECT_EQ(mesh.getPolytopeCount(d), 5);
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);

    d = 2;
    EXPECT_EQ(mesh.getPolytopeCount(d), 2);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getCellCount());
    EXPECT_EQ(mesh.getAttribute(d, 0), 1);
    EXPECT_EQ(mesh.getAttribute(d, 1), 2);
  }

  TEST(Rodin_IO_MeshLoader, SanityTest_MFEM_2D_StarSquare)
  {
    static constexpr const char* filename = "mfem/StarSquare.mfem.mesh";
    static constexpr size_t vcount = 101;
    static constexpr size_t ecount = 80;
    static constexpr size_t sdim = 2;

    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);

    boost::filesystem::ifstream in(meshfile);

    Mesh mesh;
    MeshLoader<FileFormat::MFEM, Rodin::Context::Local> loader(mesh);
    loader.load(in);

    EXPECT_EQ(mesh.getSpaceDimension(), sdim);

    size_t d = 0;
    EXPECT_EQ(mesh.getPolytopeCount(d), vcount);
    EXPECT_EQ(mesh.getPolytopeCount(d), mesh.getVertexCount());

    for (size_t i = 0; i < vcount; i++)
      EXPECT_EQ(mesh.getAttribute(d, i), std::nullopt);

    d = 2;
    for (size_t i = 0; i < ecount; i++)
      EXPECT_EQ(mesh.getAttribute(d, i), 1);

    for (size_t i = 0; i < ecount; i++)
      EXPECT_EQ(mesh.getGeometry(d, i), Polytope::Type::Quadrilateral);
  }
}
