/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <memory>
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

  /// @brief Cell geometry and the independently specified XDMF node counts.
  struct XDMFGeometryCase
  {
    Polytope::Type geometry;
    HDF5::U64 linearId;
    HDF5::U64 quadraticId;
    size_t quadraticNodes;
  };

  /// @brief Compares XML metadata, heavy data and evaluated fields across geometry orders.
  class XDMFGeometryOrderTest : public XDMFTest,
                                public ::testing::WithParamInterface<XDMFGeometryCase>
  {
    protected:
      /**
       * @brief Reads a visualization dataset after verifying its exact shape.
       * @param filename HDF5 file to inspect.
       * @param dataset Dataset path within the file.
       * @param shape Expected dimensions.
       * @returns Flattened values, or an empty vector after a failed shape check.
       */
      std::vector<HDF5::F64> readValues(const boost::filesystem::path& filename,
        const std::string& dataset, const std::vector<hsize_t>& shape)
      {
        const auto file =
          HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        EXPECT_TRUE(file);
        if (!file)
          return {};
        const auto actual = HDF5::readDatasetShape(file.get(), dataset);
        EXPECT_EQ(actual, shape);
        if (actual != shape)
          return {};
        size_t count = 1;
        for (const auto extent : shape)
          count *= extent;
        std::vector<HDF5::F64> values(count);
        const auto data =
          HDF5::DataSet(H5Dopen2(file.get(), dataset.c_str(), H5P_DEFAULT));
        EXPECT_TRUE(data);
        if (!data)
          return {};
        EXPECT_GE(H5Dread(data.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    values.data()),
          0);
        return values;
      }
  };

  /// @brief Linear product cells share vertices; mixed curved cells keep topology and field sampling synchronized.
  TEST_P(XDMFGeometryOrderTest, TopologyCoordinatesAndFieldsUseInterpolationOrder)
  {
    const auto parameter = GetParam();
    const Polytope::Traits traits(parameter.geometry);
    const size_t dimension = traits.getDimension();

    /// @brief Small quadratic displacement that vanishes at reference vertices.
    constexpr Real WarpAmplitude = 0.05;

    /// @brief Roundoff allowance for quadratic-basis interpolation of coordinate fields.
    constexpr Real SamplingTolerance = 1e-11;
    for (const bool curved : {false, true})
    {
      SCOPED_TRACE(curved);
      auto mesh =
        LocalMesh::UniformGrid(parameter.geometry, Array<size_t>::Constant(dimension, 3));
      if (parameter.geometry == Polytope::Type::Wedge)
      {
        const auto path = m_directory / "wedges.mesh";
        mesh.save(path, FileFormat::MEDIT);
        mesh.load(path, FileFormat::MEDIT);
      }
      for (size_t d = dimension; d > 0; --d)
        mesh.getConnectivity().compute(d, d - 1);
      const auto first = *mesh.getCell();
      std::unique_ptr<PolytopeTransformation> original(first.getTransformation().copy());
      // Retain the reference fields while changing only the geometric mapping.
      H1<2, Real> scalarSpace(std::integral_constant<size_t, 2>{}, mesh);
      H1<2, Math::SpatialVector<Real>> vectorSpace(
        std::integral_constant<size_t, 2>{}, mesh, 2);
      GridFunction scalar(scalarSpace);
      GridFunction vector(vectorSpace);
      scalar = [](const Point& point) { return 1 + point.x(); };
      vector = [](const Point& point) {
        return Math::SpatialVector<Real>{1 + point.x(), 5 + 2 * point.x()};
      };
      if (curved)
      {
        RealH1Element<2> element(parameter.geometry);
        PointCloud controls(dimension, element.getCount());
        for (size_t node = 0; node < element.getCount(); ++node)
        {
          const auto& reference = element.getNode(node);
          Math::SpatialPoint physical;
          original->transform(physical, reference);
          physical[0] += WarpAmplitude * reference[0] * (1 - reference[0]);
          for (size_t component = 0; component < dimension; ++component)
            controls(component, node) = physical[component];
        }
        mesh.setPolytopeTransformation({dimension, first.getIndex()},
          new ParametricTransformation<RealH1Element<2>>(std::move(controls), element));
      }
      EXPECT_EQ(HDF5::hasCurvedXDMFGeometry(mesh), curved);
      const size_t points = curved
        ? parameter.quadraticNodes + (mesh.getCellCount() - 1) * traits.getVertexCount()
        : mesh.getVertexCount();
      EXPECT_EQ(HDF5::getXDMFVisualizationVertexCount(mesh), points);
      const auto layout = HDF5::getXDMFTopologyLayout(mesh);
      ASSERT_FALSE(layout.isUniform);
      EXPECT_EQ(layout.topologyType, "Mixed");

      std::vector<HDF5::U64> topology;
      std::vector<Math::SpatialPoint> coordinates;
      std::vector<Real> referenceValues;
      if (!curved)
      {
        for (Index vertex = 0; vertex < static_cast<Index>(mesh.getVertexCount());
             ++vertex)
        {
          coordinates.push_back(mesh.getVertexCoordinates(vertex));
          referenceValues.push_back(1 + coordinates.back()[0]);
        }
      }
      size_t next = 0;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const bool quadratic = curved && cell->getIndex() == first.getIndex();
        topology.push_back(quadratic ? parameter.quadraticId : parameter.linearId);
        if (!quadratic && parameter.geometry == Polytope::Type::Segment)
          topology.push_back(2);
        if (curved)
        {
          const auto nodes =
            HDF5::getXDMFReferenceNodes(parameter.geometry, quadratic ? 2 : 1);
          for (const auto& reference : nodes)
          {
            topology.push_back(next++);
            Math::SpatialPoint physical;
            if (quadratic)
            {
              original->transform(physical, reference);
              referenceValues.push_back(1 + physical[0]);
              physical[0] += WarpAmplitude * reference[0] * (1 - reference[0]);
            }
            else
            {
              cell->getTransformation().transform(physical, reference);
              referenceValues.push_back(1 + physical[0]);
            }
            coordinates.push_back(physical);
          }
        }
        else
        {
          for (const Index vertex : cell->getVertices())
            topology.push_back(vertex);
        }
      }
      ASSERT_EQ(coordinates.size(), points);
      ASSERT_EQ(referenceValues.size(), points);
      EXPECT_EQ(layout.entryCount, topology.size());
      EXPECT_EQ(HDF5::getXDMFMixedTopologySize(mesh), topology.size());

      const std::string name = curved ? "curved" : "linear";
      XDMF writer(m_directory / name);
      writer.setMesh(mesh)
        .add("scalar", scalar, XDMF::Center::Node, XDMF::AttributePolicy::Static)
        .add("vector", vector, XDMF::Center::Node, XDMF::AttributePolicy::Static);
      writer.write(0).close();
      std::ifstream input((m_directory / (name + ".xdmf")).string());
      ASSERT_TRUE(input.good());
      const std::string xml((std::istreambuf_iterator<char>(input)), {});
      EXPECT_NE(xml.find("TopologyType=\"Mixed\""), std::string::npos);
      EXPECT_NE(
        xml.find("NumberOfElements=\"" + std::to_string(mesh.getCellCount()) + "\""),
        std::string::npos);
      EXPECT_NE(xml.find("Dimensions=\"" + std::to_string(points) + " " +
                  std::to_string(dimension) + "\""),
        std::string::npos);
      EXPECT_NE(
        xml.find("Dimensions=\"" + std::to_string(points) + "\""), std::string::npos);
      EXPECT_NE(
        xml.find("Dimensions=\"" + std::to_string(points) + " 2\""), std::string::npos);
      const auto meshPath = m_directory / (name + ".mesh.h5");
      const auto file =
        HDF5::File(H5Fopen(meshPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      ASSERT_TRUE(file);
      EXPECT_EQ(
        HDF5::readVectorDataset<HDF5::U64>(file.get(), HDF5::Path::MeshXDMFTopology),
        topology);
      const auto geometry = readValues(meshPath, HDF5::Path::MeshGeometryVertices,
        {static_cast<hsize_t>(points), static_cast<hsize_t>(dimension)});
      const auto scalars = readValues(m_directory / (name + ".scalar.h5"),
        HDF5::Path::GridFunctionValuesData, {static_cast<hsize_t>(points)});
      const auto vectors = readValues(m_directory / (name + ".vector.h5"),
        HDF5::Path::GridFunctionValuesData, {static_cast<hsize_t>(points), 2});
      ASSERT_EQ(geometry.size(), points * dimension);
      ASSERT_EQ(scalars.size(), points);
      ASSERT_EQ(vectors.size(), 2 * points);
      for (size_t point = 0; point < points; ++point)
      {
        for (size_t component = 0; component < dimension; ++component)
        {
          EXPECT_NEAR(geometry[point * dimension + component],
            coordinates[point][component], SamplingTolerance);
        }
        EXPECT_NEAR(scalars[point], referenceValues[point], SamplingTolerance);
        EXPECT_NEAR(vectors[2 * point], referenceValues[point], SamplingTolerance);
        EXPECT_NEAR(
          vectors[2 * point + 1], 3 + 2 * referenceValues[point], SamplingTolerance);
      }
    }
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, XDMFGeometryOrderTest,
    ::testing::Values(XDMFGeometryCase{Polytope::Type::Segment, 2, 34, 3},
      XDMFGeometryCase{Polytope::Type::Triangle, 4, 36, 6},
      XDMFGeometryCase{Polytope::Type::Quadrilateral, 5, 35, 9},
      XDMFGeometryCase{Polytope::Type::Tetrahedron, 6, 38, 10},
      XDMFGeometryCase{Polytope::Type::Pyramid, 7, 39, 13},
      XDMFGeometryCase{Polytope::Type::Wedge, 8, 41, 18},
      XDMFGeometryCase{Polytope::Type::Hexahedron, 9, 50, 27}),
    [](const ::testing::TestParamInfo<XDMFGeometryCase>& info) {
      std::ostringstream label;
      label << info.param.geometry;
      return label.str();
    });

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

namespace Rodin::Tests::Unit
{
  /// @brief Mesh compatibility is independent of attribute/mesh registration order.
  TEST_F(XDMFTest, AttributeMeshBindingIsOrderIndependent)
  {
    const auto small = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    const auto large = LocalMesh::UniformGrid(Polytope::Type::Segment, {5});
    const LocalMesh sameSizeDifferentMesh(small);
    P1<Real> space(small);
    GridFunction field(space);
    field = RealFunction(7);
    for (const bool attributeFirst : {false, true})
    {
      XDMF writer(m_directory / (attributeFirst ? "attribute-first" : "mesh-first"));
      auto grid = writer.grid();
      if (attributeFirst)
      {
        grid.add("value", field);
        EXPECT_THROW(grid.setMesh(large), Alert::Exception);
        EXPECT_THROW(grid.setMesh(sameSizeDifferentMesh), Alert::Exception);
        EXPECT_FALSE(grid.hasMesh());
        grid.setMesh(small);
      }
      else
      {
        grid.setMesh(large);
        EXPECT_THROW(grid.add("value", field), Alert::Exception);
        grid.setMesh(small).add("value", field);
      }
      EXPECT_THROW(grid.setMesh(large), Alert::Exception);
      EXPECT_THROW(grid.setMesh(sameSizeDifferentMesh), Alert::Exception);
      writer.write(0);
      grid.clear().setMesh(large);
      writer.write(1).close();
      EXPECT_EQ(writer.getSnapshotCount(), 2u);
    }
  }
}

namespace Rodin::Tests::Unit
{
  /// @brief Names are escaped in XML while generated filesystem paths retain their literal names.
  TEST_F(XDMFTest, EscapedNamesAndPaths)
  {
    const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    P1 fes(mesh);
    GridFunction field(fes);
    field = Real(2);
    XDMF writer(m_directory / "series&data");
    writer.grid("mesh&<\"'>").setMesh(mesh).add("value&<\"'>", field);
    writer.write().close();
    std::ifstream input((m_directory / "series&data.xdmf").string());
    ASSERT_TRUE(input.good());
    const std::string xml((std::istreambuf_iterator<char>(input)), {});
    EXPECT_NE(xml.find("mesh&amp;&lt;&quot;&apos;&gt;"), std::string::npos);
    EXPECT_NE(xml.find("value&amp;&lt;&quot;&apos;&gt;"), std::string::npos);
    EXPECT_NE(xml.find("series&amp;data"), std::string::npos);
    EXPECT_EQ(xml.find("Name=\"mesh&<"), std::string::npos);
    size_t files = 0;
    for (const auto& entry : boost::filesystem::directory_iterator(m_directory))
    {
      if (entry.path().extension() == ".h5")
      {
        ++files;
        EXPECT_NE(
          entry.path().filename().string().find("series&data"), std::string::npos);
        auto file =
          HDF5::File(H5Fopen(entry.path().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        EXPECT_TRUE(file);
      }
    }
    EXPECT_EQ(files, 2u);
  }
}
