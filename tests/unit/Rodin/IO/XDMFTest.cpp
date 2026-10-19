/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <array>
#include <cstdint>
#include <sstream>
#include <set>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include <hdf5.h>

#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::IO;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Isolates generated XDMF and HDF5 files for each test.
  class XDMFTest : public ::testing::Test
  {
    protected:
      /// @brief Creates an independent output directory.
      void SetUp() override
      {
        m_directory = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-xdmf-%%%%-%%%%");
        boost::filesystem::create_directory(m_directory);
      }

      /// @brief Removes all generated output files.
      void TearDown() override
      {
        boost::filesystem::remove_all(m_directory);
      }

      boost::filesystem::path m_directory; ///< Per-test output directory.
  };

  /// @brief Failed publication preserves open state and history for flush/close recovery.
  TEST_F(XDMFTest, OutputOpenFailurePreservesHistoryForRecovery)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    XDMF writer(m_directory / "series");
    writer.setMesh(mesh).write(0.25);
    const auto output = m_directory / "series.xdmf";
    boost::filesystem::create_directory(output);
    try
    {
      writer.close();
      FAIL() << "A directory was accepted as an XML output file";
    }
    catch (const Alert::Exception& error)
    {
      const std::string message = error.what();
      EXPECT_NE(message.find("In member function flush"), std::string::npos);
      EXPECT_NE(message.find("Failed to open XDMF output file"), std::string::npos);
      EXPECT_NE(message.find(output.string()), std::string::npos);
    }
    ASSERT_FALSE(writer.isClosed());
    EXPECT_EQ(writer.getSnapshotCount(), 1u);
    ASSERT_TRUE(boost::filesystem::remove(output));
    writer.flush();
    EXPECT_FALSE(writer.isClosed());
    writer.write(0.75);
    writer.flush();
    EXPECT_FALSE(writer.isClosed());
    EXPECT_EQ(writer.getSnapshotCount(), 2u);
    {
      std::ifstream input(output.string());
      ASSERT_TRUE(input.good());
      const std::string xml((std::istreambuf_iterator<char>(input)), {});
      EXPECT_NE(xml.find("Time Value=\"0.25\""), std::string::npos);
      EXPECT_NE(xml.find("Time Value=\"0.75\""), std::string::npos);
    }
    writer.close();
    EXPECT_TRUE(writer.isClosed());
    ASSERT_TRUE(boost::filesystem::remove(output));
    writer.close();
    EXPECT_FALSE(boost::filesystem::exists(output));
    EXPECT_THROW(writer.write(1.0), Alert::Exception);
    EXPECT_EQ(writer.getSnapshotCount(), 2u);
  }

  /// @brief Three-dimensional metadata and raw HDF5 coordinates/topology match an analytic tetrahedron.
  TEST_F(XDMFTest, ThreeDimensionalGeometryAndAffineSamples)
  {
    auto mesh = LocalMesh::Builder()
                  .initialize(3)
                  .nodes(4)
                  .vertex({0, 0, 0})
                  .vertex({2, 0, 0})
                  .vertex({0, 3, 0})
                  .vertex({0, 0, 4})
                  .polytope(Polytope::Type::Tetrahedron, {{0, 1, 2, 3}})
                  .finalize();
    P1 space(mesh);
    GridFunction field(space);
    field = [](const Point& point) {
      return 1 + point.x() + 2 * point.y() + 4 * point.z();
    };
    XDMF writer(m_directory / "volume");
    XDMF::FilePatterns patterns;
    patterns.staticMesh = "mesh.h5";
    patterns.attribute = "value.{index}.h5";
    writer.setFilePatterns(patterns).setPadding(1).setMesh(mesh).add("value", field);
    writer.write(0.0).close();
    std::ifstream input((m_directory / "volume.xdmf").string());
    ASSERT_TRUE(input.good());
    const std::string xml((std::istreambuf_iterator<char>(input)), {});
    EXPECT_NE(xml.find("GeometryType=\"XYZ\""), std::string::npos);
    EXPECT_NE(xml.find("Dimensions=\"4 3\""), std::string::npos);
    EXPECT_NE(
      xml.find("TopologyType=\"Mixed\" NumberOfElements=\"1\""), std::string::npos);
    EXPECT_NE(xml.find("Dimensions=\"5\""), std::string::npos);
    EXPECT_NE(xml.find("mesh.h5:/Mesh/Geometry/Vertices"), std::string::npos);
    EXPECT_NE(xml.find("value.0.h5:/GridFunction/Values/Data"), std::string::npos);
    const auto file =
      HDF5::File(H5Fopen((m_directory / "mesh.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    ASSERT_TRUE(file);
    const auto vertices =
      HDF5::DataSet(H5Dopen2(file.get(), "/Mesh/Geometry/Vertices", H5P_DEFAULT));
    ASSERT_TRUE(vertices);
    const auto vertexSpace = HDF5::Space(H5Dget_space(vertices.get()));
    ASSERT_TRUE(vertexSpace);
    ASSERT_EQ(H5Sget_simple_extent_ndims(vertexSpace.get()), 2);
    std::array<hsize_t, 2> shape{};
    ASSERT_EQ(H5Sget_simple_extent_dims(vertexSpace.get(), shape.data(), nullptr), 2);
    ASSERT_EQ(shape[0], 4u);
    ASSERT_EQ(shape[1], 3u);
    std::array<Real, 12> coordinates{};
    ASSERT_GE(H5Dread(vertices.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                coordinates.data()),
      0);
    std::set<std::array<Real, 3>> actualCoordinates;
    for (size_t i = 0; i < 4; ++i)
    {
      actualCoordinates.insert(
        {coordinates[3 * i], coordinates[3 * i + 1], coordinates[3 * i + 2]});
    }
    const std::set<std::array<Real, 3>> expectedCoordinates{
      {0, 0, 0}, {2, 0, 0}, {0, 3, 0}, {0, 0, 4}};
    EXPECT_EQ(actualCoordinates, expectedCoordinates);
    const auto topology =
      HDF5::DataSet(H5Dopen2(file.get(), "/Mesh/XDMF/Topology", H5P_DEFAULT));
    ASSERT_TRUE(topology);
    const auto topologySpace = HDF5::Space(H5Dget_space(topology.get()));
    ASSERT_TRUE(topologySpace);
    ASSERT_EQ(H5Sget_simple_extent_ndims(topologySpace.get()), 1);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(topologySpace.get(), &count, nullptr), 1);
    ASSERT_EQ(count, 5u);
    std::array<std::uint64_t, 5> indices{};
    ASSERT_GE(H5Dread(topology.get(), H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                indices.data()),
      0);
    EXPECT_EQ(indices[0], 6u); // External XDMF mixed-topology tetrahedron identifier.
    EXPECT_EQ((std::set<std::uint64_t>(indices.begin() + 1, indices.end())),
      (std::set<std::uint64_t>{0, 1, 2, 3}));
    const auto attributeFile = HDF5::File(
      H5Fopen((m_directory / "value.0.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    ASSERT_TRUE(attributeFile);
    const auto attribute = HDF5::DataSet(
      H5Dopen2(attributeFile.get(), "/GridFunction/Values/Data", H5P_DEFAULT));
    ASSERT_TRUE(attribute);
    const auto attributeSpace = HDF5::Space(H5Dget_space(attribute.get()));
    ASSERT_TRUE(attributeSpace);
    ASSERT_EQ(H5Sget_simple_extent_ndims(attributeSpace.get()), 1);
    ASSERT_EQ(H5Sget_simple_extent_dims(attributeSpace.get(), &count, nullptr), 1);
    ASSERT_EQ(count, 4u);
    std::array<Real, 4> values{};
    ASSERT_GE(H5Dread(attribute.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                values.data()),
      0);
    for (size_t i = 0; i < 4; ++i)
    {
      EXPECT_EQ(values[i],
        1 + coordinates[3 * i] + 2 * coordinates[3 * i + 1] + 4 * coordinates[3 * i + 2]);
    }
  }

  /// @brief Empty filename placeholders normalize repeated separators and preserve static data.
  TEST_F(XDMFTest, EmptyPlaceholdersNormalizeFileReferences)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1 space(mesh);
    GridFunction field(space);
    field = 4.0;
    boost::filesystem::create_directory(m_directory / "data");
    XDMF writer(m_directory / "series");
    XDMF::FilePatterns patterns;
    patterns.xdmf = "./{grid}/{rank}//{stem}.{rank}.xdmf";
    patterns.staticMesh = "data//./{grid}/{rank}/mesh.{index}.h5";
    patterns.staticAttribute = "data//./{grid}/{rank}/{name}.{index}.h5";
    writer.setFilePatterns(patterns).setMesh(mesh).add(
      "value", field, XDMF::Center::Node, XDMF::AttributePolicy::Static);
    writer.write(0.0);
    field = 99.0;
    writer.write(1.0).close();
    EXPECT_TRUE(boost::filesystem::exists(m_directory / "data/mesh.h5"));
    EXPECT_TRUE(boost::filesystem::exists(m_directory / "data/value.h5"));
    std::ifstream input((m_directory / "series.xdmf").string());
    ASSERT_TRUE(input.good());
    const std::string xml((std::istreambuf_iterator<char>(input)), {});
    EXPECT_NE(xml.find("data/mesh.h5:"), std::string::npos);
    EXPECT_NE(xml.find("data/value.h5:/GridFunction/Values/Data"), std::string::npos);
    EXPECT_EQ(xml.find("data//"), std::string::npos);
    EXPECT_EQ(xml.find("data/./"), std::string::npos);
    EXPECT_EQ(xml.find("..h5"), std::string::npos);
    const auto file = HDF5::File(
      H5Fopen((m_directory / "data/value.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    ASSERT_TRUE(file);
    const auto dataset =
      HDF5::DataSet(H5Dopen2(file.get(), "/GridFunction/Values/Data", H5P_DEFAULT));
    ASSERT_TRUE(dataset);
    const auto shape = HDF5::Space(H5Dget_space(dataset.get()));
    ASSERT_TRUE(shape);
    ASSERT_EQ(H5Sget_simple_extent_ndims(shape.get()), 1);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(shape.get(), &count, nullptr), 1);
    ASSERT_EQ(count, 3u);
    std::array<Real, 3> values{};
    ASSERT_GE(H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                values.data()),
      0);
    for (const Real value : values)
      EXPECT_EQ(value, 4.0);
  }

  /// @brief Vector XML shapes and static HDF5 samples match independent node and cell-center values.
  TEST_F(XDMFTest, VectorMetadataAndSamplesPreserveBothCenters)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(0.5);
    P1 space(mesh, 2);
    GridFunction field(space);
    field = [](const Point& point) {
      return Math::SpatialVector<Real>{1 + 2 * point.x(), 3 - 4 * point.x()};
    };
    XDMF writer(m_directory / "vectors");
    writer.setMesh(mesh)
      .add("nodal", field, XDMF::Center::Node, XDMF::AttributePolicy::Static)
      .add("cell", field, XDMF::Center::Cell, XDMF::AttributePolicy::Static);
    writer.write(0.0);
    // A later constant must not replace the static visualization samples.
    field = []([[maybe_unused]] const Point& point) {
      return Math::SpatialVector<Real>{99, 99};
    };
    writer.write(1.0).close();
    std::ifstream input((m_directory / "vectors.xdmf").string());
    ASSERT_TRUE(input.good());
    const std::string xml((std::istreambuf_iterator<char>(input)), {});
    for (const bool nodal : {false, true})
    {
      const std::string name = nodal ? "nodal" : "cell";
      const std::string center = nodal ? "Node" : "Cell";
      const size_t count = nodal ? 3 : 2;
      SCOPED_TRACE(name);
      const std::string filename = "vectors." + name + ".h5";
      const std::string reference = filename + ":/GridFunction/Values/Data";
      const std::string opening = "<Attribute Name=\"" + name +
        "\" AttributeType=\"Vector\" Center=\"" + center + "\">";
      size_t snapshots = 0, position = 0;
      while ((position = xml.find(opening, position)) != std::string::npos)
      {
        const auto end = xml.find("</Attribute>", position);
        ASSERT_NE(end, std::string::npos);
        const auto block = xml.substr(position, end - position);
        EXPECT_NE(block.find("Dimensions=\"" + std::to_string(count) + " 2\""),
          std::string::npos);
        EXPECT_NE(block.find(reference), std::string::npos);
        ++snapshots;
        position = end;
      }
      EXPECT_EQ(snapshots, 2u);
      const auto file = HDF5::File(
        H5Fopen((m_directory / filename).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      ASSERT_TRUE(file);
      const auto dataset =
        HDF5::DataSet(H5Dopen2(file.get(), "/GridFunction/Values/Data", H5P_DEFAULT));
      ASSERT_TRUE(dataset);
      const auto shape = HDF5::Space(H5Dget_space(dataset.get()));
      ASSERT_TRUE(shape);
      ASSERT_EQ(H5Sget_simple_extent_ndims(shape.get()), 2);
      std::array<hsize_t, 2> dimensions{};
      ASSERT_EQ(H5Sget_simple_extent_dims(shape.get(), dimensions.data(), nullptr), 2);
      ASSERT_EQ(dimensions[0], count);
      ASSERT_EQ(dimensions[1], 2u);
      std::vector<Real> values(2 * count);
      ASSERT_GE(H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  values.data()),
        0);
      for (size_t row = 0; row < count; ++row)
      {
        const Real x = nodal ? Real(row) / 2 : Real(row) / 2 + 0.25;
        EXPECT_EQ(values[2 * row], 1 + 2 * x);
        EXPECT_EQ(values[2 * row + 1], 3 - 4 * x);
      }
    }
  }

  /// @brief Clearing preserves mesh/options; resetting removes mesh and history.
  TEST_F(XDMFTest, ClearAndResetHaveDistinctPersistenceEffects)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1 fes(mesh);
    GridFunction field(fes);
    field = 2.0;
    XDMF writer(m_directory / "series");
    auto grid = writer.grid("domain");
    auto alias = writer.grid("domain");
    EXPECT_EQ(writer.getGridCount(), 1u);
    EXPECT_EQ(grid.getName(), "domain");
    EXPECT_FALSE(grid.hasMesh());
    XDMF::GridOptions options;
    options.exportMeshAttribute = false;
    EXPECT_EQ(&grid.setMesh(mesh).setOptions(options).add("old", field), &grid);
    writer.write(0.0);
    EXPECT_EQ(alias.getAttributeCount(), 1u);
    EXPECT_EQ(&alias.clear(), &alias);
    EXPECT_EQ(grid.getAttributeCount(), 0u);
    EXPECT_TRUE(grid.hasMesh());
    EXPECT_FALSE(grid.getOptions().exportMeshAttribute);
    writer.write(1.0);
    writer.flush();
    std::ifstream before((m_directory / "series.xdmf").string());
    std::ostringstream oldXML;
    oldXML << before.rdbuf();
    EXPECT_NE(oldXML.str().find("Name=\"old\""), std::string::npos);
    EXPECT_EQ(&grid.reset(), &grid);
    EXPECT_FALSE(alias.hasMesh());
    EXPECT_EQ(alias.getAttributeCount(), 0u);
    EXPECT_TRUE(alias.getOptions().exportMeshAttribute);
    EXPECT_EQ(alias.getOptions().meshPolicy, XDMF::MeshPolicy::Static);
    EXPECT_FALSE(alias.getOptions().patterns);
    EXPECT_THROW(writer.write(2.0), Alert::Exception);
    EXPECT_EQ(writer.getSnapshotCount(), 2u);
    grid.setMesh(mesh).add("new", field);
    writer.write(3.0).close();
    EXPECT_TRUE(writer.isClosed());
    EXPECT_NO_THROW(writer.close());
    EXPECT_THROW(writer.write(), Alert::Exception);
    std::ifstream after((m_directory / "series.xdmf").string());
    std::ostringstream newXML;
    newXML << after.rdbuf();
    EXPECT_EQ(newXML.str().find("Name=\"old\""), std::string::npos);
    EXPECT_EQ(newXML.str().find("Time Value=\"0\""), std::string::npos);
    EXPECT_NE(newXML.str().find("Name=\"new\""), std::string::npos);
    EXPECT_NE(newXML.str().find("Time Value=\"3\""), std::string::npos);
  }

  /// @brief Per-grid patterns and padding control transient file references.
  TEST_F(XDMFTest, GridPatternsOverrideWriterPatterns)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1 fes(mesh);
    GridFunction field(fes);
    field = 4.0;
    XDMF writer(m_directory / "series");
    XDMF::FilePatterns patterns;
    patterns.xdmf = "{stem}-catalog.xdmf";
    patterns.transientMesh = "global-{grid}-{index}.h5";
    writer.setFilePatterns(patterns).setPadding(3);
    EXPECT_EQ(writer.getPadding(), 3u);
    EXPECT_EQ(writer.getFilePatterns().xdmf, "{stem}-catalog.xdmf");
    EXPECT_EQ(writer.getStem(), m_directory / "series");
    XDMF::GridOptions options;
    options.meshPolicy = XDMF::MeshPolicy::Transient;
    options.exportMeshAttribute = false;
    patterns.transientMesh = "{stem}-{grid}-mesh-{index}.h5";
    patterns.attribute = "{stem}-{grid}-{name}-{index}.h5";
    options.patterns = patterns;
    auto grid = writer.grid("region");
    grid.setMesh(mesh).setOptions(options).add("value", field);
    writer.write(0.25).write(0.5).close();
    for (const auto& index : {"000", "001"})
    {
      EXPECT_TRUE(boost::filesystem::exists(
        m_directory / (std::string("series-region-mesh-") + index + ".h5")));
      EXPECT_TRUE(boost::filesystem::exists(
        m_directory / (std::string("series-region-value-") + index + ".h5")));
    }
    EXPECT_FALSE(boost::filesystem::exists(m_directory / "global-region-000.h5"));
    std::ifstream input((m_directory / "series-catalog.xdmf").string());
    ASSERT_TRUE(input.good());
    std::ostringstream xml;
    xml << input.rdbuf();
    EXPECT_NE(xml.str().find("series-region-mesh-000.h5:"), std::string::npos);
    EXPECT_NE(xml.str().find("series-region-value-001.h5:"), std::string::npos);
    EXPECT_EQ(xml.str().find("Name=\"Attribute\""), std::string::npos);
    EXPECT_EQ(writer.getSnapshotCount(), 2u);
    EXPECT_FALSE(writer.isDistributed());
    EXPECT_EQ(writer.getRank(), 0u);
    EXPECT_EQ(writer.getNumRanks(), 1u);
    EXPECT_EQ(writer.getRootRank(), 0u);
  }

  /// @brief Static attributes retain the first snapshot despite subsequent field changes.
  TEST_F(XDMFTest, StaticAttributeIsWrittenOnlyOnce)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1 fes(mesh);
    GridFunction field(fes);
    field = 4.0;
    XDMF writer(m_directory / "series");
    XDMF::FilePatterns patterns;
    patterns.staticAttribute = "{stem}-{name}-fixed.h5";
    writer.setFilePatterns(patterns);
    writer.setMesh(mesh).add(
      "value", field, XDMF::Center::Node, XDMF::AttributePolicy::Static);
    writer.write();
    field = 9.0;
    writer.write().close();
    EXPECT_EQ(writer.getSnapshotCount(), 2u);
    const auto path = m_directory / "series-value-fixed.h5";
    const hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0);
    const hid_t dataset = H5Dopen2(file, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dataset, 0);
    const hid_t space = H5Dget_space(dataset);
    ASSERT_GE(space, 0);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(space, &count, nullptr), 1);
    EXPECT_EQ(count, mesh.getVertexCount());
    std::vector<double> values(count);
    EXPECT_GE(
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
      0);
    for (const auto value : values)
      EXPECT_DOUBLE_EQ(value, 4.0);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
    std::ifstream input((m_directory / "series.xdmf").string());
    std::ostringstream xml;
    xml << input.rdbuf();
    const auto first = xml.str().find("series-value-fixed.h5:");
    ASSERT_NE(first, std::string::npos);
    EXPECT_NE(xml.str().find("series-value-fixed.h5:", first + 1), std::string::npos);
    EXPECT_NE(xml.str().find("Time Value=\"0\""), std::string::npos);
    EXPECT_NE(xml.str().find("Time Value=\"1\""), std::string::npos);
  }

  /// @brief Attribute registration rejects unnamed, duplicate and foreign fields.
  TEST_F(XDMFTest, InvalidAttributesLeaveRegistrationUnchanged)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    Mesh other = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    P1 fes(mesh);
    P1 otherFES(other);
    GridFunction field(fes);
    GridFunction foreign(otherFES);
    XDMF writer(m_directory / "series");
    auto grid = writer.grid();
    grid.setMesh(mesh);
    EXPECT_THROW(grid.add(field), Alert::Exception);
    EXPECT_THROW(grid.add("foreign", foreign), Alert::Exception);
    EXPECT_EQ(grid.getAttributeCount(), 0u);
    grid.add("value", field);
    EXPECT_THROW(grid.add("value", field), Alert::Exception);
    EXPECT_EQ(grid.getAttributeCount(), 1u);
    EXPECT_EQ(grid.getName(), "");
  }
}

TEST(XDMFGridHandle, UnregisteredIndexRaises)
{
  IO::XDMF writer("unregistered-grid");
  EXPECT_THROW((IO::XDMF::Grid(writer, 0)), Alert::Exception);
}
