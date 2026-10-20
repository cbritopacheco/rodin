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
  /// @brief Member diagnostics report physical lines and restart on each load attempt.
  TEST(Rodin_IO_MeshLoader, MFEMDiagnosticLinesAndLoaderReuse)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
    size_t expectedLine = 6;
    for (const std::string input :
      {"# comment\n\nMFEM mesh v1.0\n# marker\n\nvertices\n2\n",
        "MFEM mesh v1.0\nvertices\n2\n"})
    {
      std::stringstream stream(input);
      try
      {
        loader.load(stream);
        FAIL() << "Missing dimension keyword was accepted";
      }
      catch (const Alert::Exception& error)
      {
        const std::string message = error.what();
        EXPECT_NE(message.find("In member function readDimension"), std::string::npos);
        EXPECT_NE(message.find("on line " + std::to_string(expectedLine) + "."),
          std::string::npos);
      }
      EXPECT_EQ(mesh.getVertexCount(), 4u);
      EXPECT_EQ(mesh.getCellCount(), 2u);
      expectedLine = 2;
    }
  }

  /// @brief A valid-looking section cannot replace the required dimension keyword.
  TEST(Rodin_IO_MeshLoader, MFEMRequiresDimensionKeyword)
  {
    for (const std::string keyword : {"elements", "boundary", "vertices"})
    {
      SCOPED_TRACE(keyword);
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
      std::stringstream stream("MFEM mesh v1.0\n" + keyword +
        "\n2\nelements\n1\n7 2 0 1 2\nboundary\n0\nvertices\n3\n2\n0 0\n1 0\n0 1\n");
      MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
      EXPECT_THROW(loader.load(stream), Alert::Exception);
      EXPECT_EQ(mesh.getVertexCount(), 4u);
      EXPECT_EQ(mesh.getCellCount(), 2u);
    }
  }

  /// @brief Straight MFEM meshes require a separate spatial-dimension record.
  TEST(Rodin_IO_MeshLoader, MFEMRequiresExplicitSpaceDimension)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    std::stringstream stream("MFEM mesh v1.0\ndimension\n2\nelements\n1\n"
                             "7 2 0 1 2\nboundary\n0\nvertices\n3\n0 0\n1 0\n0 1\n");
    MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
    EXPECT_THROW(loader.load(stream), Alert::Exception);
    EXPECT_EQ(mesh.getVertexCount(), 4u);
    EXPECT_EQ(mesh.getCellCount(), 2u);
  }

  /// @brief Malformed keywords and spatial dimensions are rejected before replacing the mesh.
  TEST(Rodin_IO_MeshLoader, MFEMRejectsInvalidDimensionRecords)
  {
    for (const std::string input :
      {"MFEM mesh v1.0\n!\n2\n", "MFEM mesh v1.0\n", "MFEM mesh v1.0\ndimension\n4\n",
        "MFEM mesh v1.0\ndimension\n2\nvertices\n1\n0\n",
        "MFEM mesh v1.0\ndimension\n2\nvertices\n1\n1\n",
        "MFEM mesh v1.0\ndimension\n2\nvertices\n1\n4\n"})
    {
      SCOPED_TRACE(input);
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
      std::stringstream stream(input);
      MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
      EXPECT_THROW(loader.load(stream), Alert::Exception);
      EXPECT_EQ(mesh.getVertexCount(), 4u);
      EXPECT_EQ(mesh.getCellCount(), 2u);
    }
  }

  /// @brief Topological and spatial dimensions remain distinct for embedded surfaces.
  TEST(Rodin_IO_MeshLoader, MFEMEmbeddedTriangleCoordinatesAndArea)
  {
    Mesh mesh;
    std::stringstream stream(
      "MFEM mesh v1.0\ndimension\n2\nelements\n1\n"
      "7 2 0 1 2\nboundary\n0\nvertices\n3\n3\n0 0 1\n1 0 1\n0 1 1\n");
    MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
    loader.load(stream);
    ASSERT_EQ(mesh.getDimension(), 2u);
    ASSERT_EQ(mesh.getSpaceDimension(), 3u);
    ASSERT_EQ(mesh.getVertexCount(), 3u);
    ASSERT_EQ(mesh.getCellCount(), 1u);
    for (Index vertex = 0; vertex < 3; ++vertex)
    {
      const auto& point = mesh.getVertexCoordinates(vertex);
      EXPECT_EQ(point[0], vertex == 1 ? 1 : 0);
      EXPECT_EQ(point[1], vertex == 2 ? 1 : 0);
      EXPECT_EQ(point[2], 1);
    }
    auto cell = mesh.getCell();
    EXPECT_EQ(cell->getAttribute(), 7);
    EXPECT_NEAR(cell->getMeasure(), 0.5, 1e-12);
  }

  /// @brief Rejects malformed MFEM headers and dimensions without replacing the destination.
  TEST(Rodin_IO_MeshLoader, MFEMRejectsMalformedHeader)
  {
    for (const std::string input :
      {"not an MFEM mesh\n", "MFEM mesh v1.0\ndimension\nnotanumber\n"})
    {
      SCOPED_TRACE(input);
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
      const auto vertexCount = mesh.getVertexCount();
      const auto cellCount = mesh.getCellCount();
      std::stringstream stream(input);
      MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
      EXPECT_THROW(loader.load(stream), Alert::Exception);
      EXPECT_EQ(mesh.getVertexCount(), vertexCount);
      EXPECT_EQ(mesh.getCellCount(), cellCount);
    }
  }

  /// @brief Rejects malformed and truncated MFEM sections before committing a partial mesh.
  TEST(Rodin_IO_MeshLoader, MFEMRejectsMalformedSections)
  {
    for (const std::string section : {"!\n", "unknown\n", "dimension\n2\n",
           "elements\ninvalid\n", "elements\n1\ninvalid\n", "elements\n1\n",
           "boundary\ninvalid\n", "boundary\n1\ninvalid\n", "boundary\n1\n",
           "vertices\ninvalid\n", "vertices\n1\n2\ninvalid\n", "vertices\n1\n2\n"})
    {
      SCOPED_TRACE(section);
      Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
      const auto vertexCount = mesh.getVertexCount();
      const auto cellCount = mesh.getCellCount();
      std::stringstream stream("MFEM mesh v1.0\ndimension\n2\n" + section);
      MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
      EXPECT_THROW(loader.load(stream), Alert::Exception);
      EXPECT_EQ(mesh.getVertexCount(), vertexCount);
      EXPECT_EQ(mesh.getCellCount(), cellCount);
    }
  }

  namespace
  {
    Mesh<Context::Local> makePointMesh()
    {
      Mesh<Context::Local>::Builder builder;
      return builder.initialize(1).nodes(1).vertex({0.0}).attribute({0, 0}, 7).finalize();
    }

    Mesh<Context::Local> makeCoverageMesh(Polytope::Type geometry)
    {
      switch (geometry)
      {
        case Polytope::Type::Point:
          return makePointMesh();
        case Polytope::Type::Segment:
          return LocalMesh::UniformGrid(geometry, {5});
        case Polytope::Type::Triangle:
        case Polytope::Type::Quadrilateral:
          return LocalMesh::UniformGrid(geometry, {3, 3});
        case Polytope::Type::Tetrahedron:
        case Polytope::Type::Pyramid:
        case Polytope::Type::Hexahedron:
        case Polytope::Type::Wedge:
          return LocalMesh::UniformGrid(geometry, {2, 2, 2});
      }
      assert(false);
      return {};
    }

    Mesh<Context::Local> makeAttributeRegressionMesh(size_t dim)
    {
      Mesh<Context::Local>::Builder builder;
      builder.initialize(dim);

      switch (dim)
      {
        case 1:
        {
          builder.nodes(3)
            .vertex({0.0})
            .vertex({1.0})
            .vertex({2.0})
            .polytope(Polytope::Type::Segment, {0, 1})
            .attribute({1, 0}, 11)
            .polytope(Polytope::Type::Segment, {1, 2})
            .attribute({1, 1}, 22);
          break;
        }
        case 2:
        {
          builder.nodes(5)
            .vertex({0.0, 0.0})
            .vertex({1.0, 0.0})
            .vertex({0.0, 1.0})
            .vertex({2.0, 0.0})
            .vertex({2.0, 1.0})
            .polytope(Polytope::Type::Triangle, {0, 1, 2})
            .attribute({2, 0}, 11)
            .polytope(Polytope::Type::Quadrilateral, {1, 3, 4, 2})
            .attribute({2, 1}, 22);
          break;
        }
        case 3:
        {
          builder.nodes(8)
            .vertex({0.0, 0.0, 0.0})
            .vertex({1.0, 0.0, 0.0})
            .vertex({1.0, 1.0, 0.0})
            .vertex({0.0, 1.0, 0.0})
            .vertex({0.0, 0.0, 1.0})
            .vertex({1.0, 0.0, 1.0})
            .vertex({1.0, 1.0, 1.0})
            .vertex({0.0, 1.0, 1.0})
            .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 4})
            .attribute({3, 0}, 11)
            .polytope(Polytope::Type::Hexahedron, {0, 1, 2, 3, 4, 5, 6, 7})
            .attribute({3, 1}, 22);
          break;
        }
        default:
          assert(false);
      }
      return builder.finalize();
    }

    std::string geometryName(Polytope::Type geometry)
    {
      switch (geometry)
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
      const Mesh<Context::Local>& actual, const Mesh<Context::Local>& expected)
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
      {
        ASSERT_EQ(actual.getPolytopeCount(d), expected.getPolytopeCount(d))
          << "dimension " << d;
      }

      const size_t D = expected.getDimension();
      for (Index c = 0; c < static_cast<Index>(expected.getCellCount()); ++c)
      {
        EXPECT_EQ(actual.getGeometry(D, c), expected.getGeometry(D, c)) << "cell " << c;
        EXPECT_TRUE(
          Polytope::Key::SymmetricEquality()(actual.getConnectivity().getPolytope(D, c),
            expected.getConnectivity().getPolytope(D, c)))
          << "cell " << c;
      }
    }

    void expectSameMeshWithAttributes(
      const Mesh<Context::Local>& actual, const Mesh<Context::Local>& expected)
    {
      expectSameMeshShape(actual, expected);

      for (size_t d = 0; d <= expected.getDimension(); ++d)
      {
        for (Index i = 0; i < static_cast<Index>(expected.getPolytopeCount(d)); ++i)
        {
          EXPECT_EQ(actual.getGeometry(d, i), expected.getGeometry(d, i))
            << "polytope (" << d << ", " << i << ")";
          EXPECT_TRUE(
            Polytope::Key::SymmetricEquality()(actual.getConnectivity().getPolytope(d, i),
              expected.getConnectivity().getPolytope(d, i)))
            << "polytope (" << d << ", " << i << ")";
          EXPECT_EQ(actual.getAttribute(d, i), expected.getAttribute(d, i))
            << "polytope (" << d << ", " << i << ")";
        }
      }
    }

    void expectAttributeRegressionMesh(const Mesh<Context::Local>& mesh, size_t dim)
    {
      ASSERT_EQ(mesh.getSpaceDimension(), dim);
      ASSERT_EQ(mesh.getDimension(), dim);
      ASSERT_EQ(mesh.getCellCount(), 2);
      ASSERT_EQ(mesh.getPolytopeCount(dim), 2);

      EXPECT_EQ(mesh.getAttribute(dim, 0), 11);
      EXPECT_EQ(mesh.getAttribute(dim, 1), 22);

      switch (dim)
      {
        case 1:
          EXPECT_EQ(mesh.getGeometry(1, 0), Polytope::Type::Segment);
          EXPECT_EQ(mesh.getGeometry(1, 1), Polytope::Type::Segment);
          break;
        case 2:
          EXPECT_EQ(mesh.getGeometry(2, 0), Polytope::Type::Triangle);
          EXPECT_EQ(mesh.getGeometry(2, 1), Polytope::Type::Quadrilateral);
          break;
        case 3:
          EXPECT_EQ(mesh.getGeometry(3, 0), Polytope::Type::Tetrahedron);
          EXPECT_EQ(mesh.getGeometry(3, 1), Polytope::Type::Hexahedron);
          break;
        default:
          FAIL() << "Unsupported dimension " << dim;
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

    template <FileFormat Format>
    void expectAttributeRegressionStringRoundTrip(size_t dim)
    {
      Mesh mesh = makeAttributeRegressionMesh(dim);

      std::stringstream out;
      MeshPrinter<Format, Context::Local> printer(mesh);
      printer.print(out);

      Mesh loaded;
      MeshLoader<Format, Context::Local> loader(loaded);
      loader.load(out);

      expectAttributeRegressionMesh(loaded, dim);
    }
  }

  namespace
  {
    /// @brief A failed MEDIT parse must not commit its partially constructed mesh.
    void expectPreservedMEDITFailure(const std::string& input)
    {
      auto expected = makeAttributeRegressionMesh(2);
      computeAllIncidences(expected);
      auto mesh = expected;
      std::stringstream stream(input);
      MeshLoader<FileFormat::MEDIT, Context::Local> loader(mesh);
      EXPECT_THROW(loader.load(stream), Alert::Exception);
      ASSERT_EQ(mesh.getSpaceDimension(), expected.getSpaceDimension());
      ASSERT_EQ(mesh.getDimension(), expected.getDimension());
      for (size_t d = 0; d <= expected.getDimension(); ++d)
        ASSERT_EQ(mesh.getPolytopeCount(d), expected.getPolytopeCount(d));
      expectSameMeshWithAttributes(mesh, expected);
    }
  }

  /// @brief Missing or malformed numeric headers leave the destination intact.
  TEST(Rodin_IO_MeshLoader, MEDITRejectsMalformedHeaders)
  {
    for (const std::string input :
      {"", "MeshVersionFormatted invalid\ninvalid\n", "MeshVersionFormatted 2\n",
        "MeshVersionFormatted 2\nDimension invalid\ninvalid\n"})
    {
      SCOPED_TRACE(input);
      expectPreservedMEDITFailure(input);
    }
  }

  /// @brief Every entity reader rejects malformed counts, invalid records and truncation.
  TEST(Rodin_IO_MeshLoader, MEDITRejectsMalformedAndTruncatedEntities)
  {
    const std::string header = "MeshVersionFormatted 2\nDimension 3\n";
    const std::string vertices = "Vertices\n8\n0 0 0 0\n1 0 0 0\n1 1 0 0\n0 1 0 0\n"
                                 "0 0 1 0\n1 0 1 0\n1 1 1 0\n0 1 1 0\n";
    for (const std::string entity : {"Vertices", "Edges", "Triangles", "Quadrilaterals",
           "Tetrahedra", "Hexahedra", "Wedges", "Pyramids"})
    {
      for (const std::string suffix : {"invalid\n", "1\ninvalid\n", "1\n"})
      {
        const std::string input =
          header + (entity == "Vertices" ? "" : vertices) + entity + "\n" + suffix;
        SCOPED_TRACE(input);
        expectPreservedMEDITFailure(input);
      }
    }
    expectPreservedMEDITFailure(header + "Vertices\n2\n0 0 0 0\n");
  }

  class MeshFormatCoverage : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief Verifies MEDIT round trip for mesh format coverage.
  TEST_P(MeshFormatCoverage, MEDITRoundTrip)
  {
    expectMeshStringRoundTrip<FileFormat::MEDIT>(GetParam());
  }

  /// @brief Verifies MFEM round trip for mesh format coverage.
  TEST_P(MeshFormatCoverage, MFEMRoundTrip)
  {
    expectMeshStringRoundTrip<FileFormat::MFEM>(GetParam());
  }

  class MeshAttributeRegression : public ::testing::TestWithParam<size_t>
  {};

  /// @brief Verifies MEDIT preserves dimension attributes for mesh attribute regression.
  TEST_P(MeshAttributeRegression, MEDITPreservesDimensionAttributes)
  {
    expectAttributeRegressionStringRoundTrip<FileFormat::MEDIT>(GetParam());
  }

  /// @brief Verifies MFEM preserves dimension attributes for mesh attribute regression.
  TEST_P(MeshAttributeRegression, MFEMPreservesDimensionAttributes)
  {
    expectAttributeRegressionStringRoundTrip<FileFormat::MFEM>(GetParam());
  }

  /// @brief Verifies MEDIT load save load preserves canonical attributes for IO mesh loader by checking exact expected values.
  TEST(Rodin_IO_MeshLoader, MEDITLoadSaveLoadPreservesCanonicalAttributes)
  {
    const std::string input = "MeshVersionFormatted 1\n"
                              "Dimension 3\n"
                              "\n"
                              "Vertices\n"
                              "4\n"
                              "0 0 0 0\n"
                              "1 0 0 0\n"
                              "0 1 0 0\n"
                              "0 0 1 0\n"
                              "\n"
                              "Tetrahedra\n"
                              "1\n"
                              "1 2 3 4 7\n"
                              "\n"
                              "Triangles\n"
                              "3\n"
                              "1 2 3 11\n"
                              "1 2 3 22\n"
                              "1 2 4 33\n"
                              "\n"
                              "End\n";

    Mesh first;
    std::stringstream in(input);
    MeshLoader<FileFormat::MEDIT, Context::Local> loader(first);
    loader.load(in);

    ASSERT_EQ(first.getPolytopeCount(2), 2);
    EXPECT_EQ(first.getAttribute(2, 0), 22);
    EXPECT_EQ(first.getAttribute(2, 1), 33);
    ASSERT_EQ(first.getPolytopeCount(3), 1);
    EXPECT_EQ(first.getAttribute(3, 0), 7);

    first.scale(1);

    std::stringstream out;
    MeshPrinter<FileFormat::MEDIT, Context::Local> printer(first);
    printer.print(out);

    Mesh second;
    MeshLoader<FileFormat::MEDIT, Context::Local> reloader(second);
    reloader.load(out);

    expectSameMeshWithAttributes(second, first);
  }

  /// @brief Instantiates Mesh Attribute Regression over the Dimensions parameter coverage.
  INSTANTIATE_TEST_SUITE_P(Dimensions, MeshAttributeRegression,
    ::testing::Values(1, 2, 3), [](const ::testing::TestParamInfo<size_t>& info) {
      return "Dim" + std::to_string(info.param) + "D";
    });

  /// @brief Instantiates Mesh Format Coverage over the All Supported Geometries parameter coverage.
  INSTANTIATE_TEST_SUITE_P(AllSupportedGeometries, MeshFormatCoverage,
    ::testing::Values(Polytope::Type::Point, Polytope::Type::Segment,
      Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const ::testing::TestParamInfo<Polytope::Type>& info) {
      return geometryName(info.param);
    });

  /// @brief Verifies sanity test MEDIT 2 D square for IO mesh loader by checking exact expected values.
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

  /// @brief Verifies sanity test MFEM 2 D square for IO mesh loader by checking exact expected values.
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

  /// @brief Verifies sanity test MFEM 2 D star square for IO mesh loader by checking exact expected values.
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
