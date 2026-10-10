/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <set>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <gtest/gtest.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/HDF5.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>

namespace Rodin::Tests::Unit
{
  /// @brief Checks actual example output, including coordinates rather than only dataset shapes.
  TEST(Adaptation_SWIFT_Output, ReconstructionGeometry)
  {
    const char* output = std::getenv("SWIFT_OUTPUT");
    const char* dimension = std::getenv("SWIFT_DIMENSION");
    ASSERT_NE(output, nullptr);
    ASSERT_NE(dimension, nullptr);
    const size_t d = std::stoul(dimension);
    ASSERT_TRUE(d == 2 || d == 3);
    // Match the fixture's background and independently check minimum-cut labels.
    constexpr size_t GridPoints = 4;
    const Real h = Real(1) / Real(GridPoints - 1);
    auto mesh = d == 2
      ? Geometry::LocalMesh::UniformGrid(Geometry::Polytope::Type::Triangle,
          {GridPoints, GridPoints})
      : Geometry::LocalMesh::UniformGrid(Geometry::Polytope::Type::Tetrahedron,
          {GridPoints, GridPoints, GridPoints});
    mesh.scale(h);
    mesh.getConnectivity().compute(d - 1, d);
    Geometry::MinSTCut classifier(mesh);
    decltype(classifier)::Parameters parameters;
    parameters.smoothing = [d](const Geometry::Polytope&) {
      return d == 2 ? Real(0.008) : Real(0.004);
    };
    classifier.setParameters(parameters);
    const auto expected = classifier.classify([&](const Geometry::Polytope& cell) {
      const auto& formula = QF::PolytopeQuadratureFormula::get(2, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(formula);
      Real integral = 0;
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& point = quadrature.getPoint(q);
        Math::SpatialVector<Real> radial(point.getPhysicalCoordinates());
        for (size_t component = 0; component < d; ++component)
          radial(component) -= Real(0.5);
        const Real distance = radial.norm() - Real(0.25);
        integral += formula.getWeight(q) * point.getDistortion() *
          std::tanh(distance / (Real(1.25) * h));
      }
      return integral / cell.getMeasure();
    });
    const std::string stem(output);
    ASSERT_TRUE(std::filesystem::exists(stem + ".xdmf"));
    boost::property_tree::ptree document;
    boost::property_tree::read_xml(stem + ".xdmf", document);
    std::set<std::string> blocks;
    for (const auto& [tag, collection] : document.get_child("Xdmf.Domain"))
    {
      if (tag != "Grid")
        continue;
      const auto block = collection.get<std::string>("<xmlattr>.Name");
      EXPECT_TRUE(block == "background" || block == "moved");
      EXPECT_TRUE(blocks.insert(block).second);
      const auto& grid = collection.get_child("Grid");
      const auto file = IO::HDF5::File(
        H5Fopen((stem + "." + block + ".mesh.h5").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      ASSERT_TRUE(file);
      const auto [count, sdim] =
        IO::HDF5::readMatrixShape(file.get(), IO::HDF5::Path::MeshGeometryVertices);
      ASSERT_EQ(sdim, d);
      ASSERT_GT(count, 0);
      const auto coordinates = IO::HDF5::readVectorDataset<IO::HDF5::F64>(
        file.get(), IO::HDF5::Path::MeshGeometryVertices);
      ASSERT_EQ(coordinates.size(), count * d);
      constexpr Real CoordinateTolerance = 1e-10;
      for (size_t component = 0; component < d; ++component)
      {
        Real minimum = std::numeric_limits<Real>::infinity();
        Real maximum = -std::numeric_limits<Real>::infinity();
        for (size_t node = 0; node < count; ++node)
        {
          const Real value = coordinates[node * d + component];
          ASSERT_TRUE(std::isfinite(value));
          minimum = std::min(minimum, value);
          maximum = std::max(maximum, value);
        }
        EXPECT_GT(maximum - minimum, CoordinateTolerance);
        if (block == "background")
        {
          EXPECT_NEAR(minimum, Real(0), CoordinateTolerance);
          EXPECT_NEAR(maximum, Real(1), CoordinateTolerance);
        }
      }
      const auto topology = IO::HDF5::readVectorDataset<IO::HDF5::U64>(
        file.get(), IO::HDF5::Path::MeshXDMFTopology);
      ASSERT_FALSE(topology.empty());
      for (const auto node : topology)
        EXPECT_LT(node, count);

      const auto labels = IO::HDF5::readVectorDataset<IO::HDF5::U64>(
        file.get(), IO::HDF5::attributePath(d));
      ASSERT_FALSE(labels.empty());
      ASSERT_EQ(labels.size(), mesh.getCellCount());
      for (const Index cell : expected.inside)
        EXPECT_EQ(labels[cell], 1u);
      for (const Index cell : expected.outside)
        EXPECT_EQ(labels[cell], 2u);
      EXPECT_NE(std::find(labels.begin(), labels.end(), 1), labels.end());
      EXPECT_NE(std::find(labels.begin(), labels.end(), 2), labels.end());
      std::set<std::string> fields;
      for (const auto& [entry, attribute] : grid)
      {
        if (entry != "Attribute")
          continue;
        const auto name = attribute.get<std::string>("<xmlattr>.Name");
        fields.insert(name);
        if (name == "Attribute")
          continue;
        const auto reference = attribute.get<std::string>("DataItem");
        const auto separator = reference.rfind(':');
        ASSERT_NE(separator, std::string::npos);
        const auto path =
          std::filesystem::path(stem).parent_path() / reference.substr(0, separator);
        const auto fieldFile =
          IO::HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        ASSERT_TRUE(fieldFile);
        const auto values = IO::HDF5::readVectorDataset<IO::HDF5::F64>(
          fieldFile.get(), reference.substr(separator + 1));
        const bool nodal = name == "phi" || name == "displacement";
        EXPECT_EQ(
          attribute.get<std::string>("<xmlattr>.Center"), nodal ? "Node" : "Cell");
        EXPECT_EQ(values.size(),
          (nodal ? count : labels.size()) * (name == "displacement" ? d : 1));
        for (const auto value : values)
          EXPECT_TRUE(std::isfinite(value));
        if (name == "j" || name == "q_rel")
        {
          for (const auto value : values)
            EXPECT_GT(value, 0);
        }
        if (name == "cell_label")
        {
          ASSERT_EQ(values.size(), labels.size());
          for (size_t cell = 0; cell < labels.size(); ++cell)
            EXPECT_EQ(values[cell], labels[cell]);
        }
      }
      EXPECT_EQ(fields,
        (std::set<std::string>{
          "Attribute", "cell_label", "displacement", "phi", "j", "q_rel"}));
    }
    EXPECT_EQ(blocks, (std::set<std::string>{"background", "moved"}));
  }
}
