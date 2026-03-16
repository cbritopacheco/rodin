/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/IO.h>

#include <hdf5.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  // Helper to create a local mesh for each polytope type
  Geometry::Mesh<Context::Local> makeMesh(Polytope::Type type)
  {
    using LocalMesh = Geometry::Mesh<Context::Local>;
    switch (type)
    {
      case Polytope::Type::Segment:
        return LocalMesh::UniformGrid(type, { 11 });
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return LocalMesh::UniformGrid(type, { 4, 4 });
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
        return LocalMesh::UniformGrid(type, { 3, 3, 3 });
      default:
        ADD_FAILURE() << "Unsupported polytope type for makeMesh";
        return Geometry::Mesh<Context::Local>();
    }
  }

  std::string polytopeLabel(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Segment:       return "Segment1D";
      case Polytope::Type::Triangle:      return "Triangle2D";
      case Polytope::Type::Quadrilateral: return "Quad2D";
      case Polytope::Type::Tetrahedron:   return "Tet3D";
      case Polytope::Type::Hexahedron:    return "Hex3D";
      case Polytope::Type::Wedge:         return "Wedge3D";
      default:                            return "Unknown";
    }
  }

  // Parameterized test fixture
  class MPIMeshHDF5 : public ::testing::TestWithParam<Polytope::Type> {};

  // --- MPI Mesh HDF5 round-trip via MeshPrinter/MeshLoader -------------------

  TEST_P(MPIMeshHDF5, MeshRoundTripViaPrinterLoader)
  {
    const auto type = GetParam();
    const std::string meshFile =
        "/tmp/rodin_mpi_hdf5_rt_" + polytopeLabel(type) + ".h5";

    // Build a local mesh and wrap it in an MPI mesh (single-rank MPI)
    auto localMesh = makeMesh(type);
    ASSERT_GT(localMesh.getVertexCount(), 0u);
    ASSERT_GT(localMesh.getCellCount(), 0u);

    // Save the local mesh using MeshPrinter<HDF5, Context::Local> directly
    // (simulates what MPI MeshPrinter does by delegating to shard)
    MeshPrinter<FileFormat::HDF5, Context::Local> printer(localMesh);
    printer.print(meshFile);

    // Verify HDF5 dataset presence
    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Meta/SpaceDimension", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 1);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/XDMF", H5P_DEFAULT), 0);

    H5Fclose(h5);

    // Load back using MeshLoader<HDF5, Context::Local> and compare
    Mesh<Context::Local> loaded;
    MeshLoader<FileFormat::HDF5, Context::Local> loader(loaded);
    loader.load(meshFile);

    EXPECT_EQ(loaded.getSpaceDimension(), localMesh.getSpaceDimension());
    EXPECT_EQ(loaded.getDimension(), localMesh.getDimension());
    EXPECT_EQ(loaded.getVertexCount(), localMesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), localMesh.getCellCount());

    std::remove(meshFile.c_str());
  }

  // --- MPI Mesh persistence does not create XDMF datasets -------------------

  TEST_P(MPIMeshHDF5, MeshPersistenceNoXDMF)
  {
    const auto type = GetParam();
    const std::string meshFile =
        "/tmp/rodin_mpi_hdf5_noxdmf_" + polytopeLabel(type) + ".h5";

    auto localMesh = makeMesh(type);
    MeshPrinter<FileFormat::HDF5, Context::Local> printer(localMesh);
    printer.print(meshFile);

    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Canonical data present
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 1);

    // XDMF-specific must NOT be present
    EXPECT_EQ(H5Lexists(h5, "/Mesh/XDMF", H5P_DEFAULT), 0);

    H5Fclose(h5);
    std::remove(meshFile.c_str());
  }

  // --- MPI Mesh XDMF visualization -------------------------------------------

  TEST_P(MPIMeshHDF5, XDMFVisualizationTopology)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_mpi_xdmf_topo_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    auto localMesh = makeMesh(type);

    {
      XDMF xdmf(stem);
      // Passing a local mesh simulates what happens with an MPI mesh shard
      xdmf.setMesh(localMesh);
      xdmf.write(0.0);
      xdmf.close();
    }

    const auto meshH5 = testDir / "vis.mesh.h5";
    ASSERT_TRUE(boost::filesystem::exists(meshH5));

    hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Visualization-only data present
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/TopologySize", H5P_DEFAULT), 1);

    // Canonical persistence must NOT be present
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 0);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Transformations", H5P_DEFAULT), 0);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Meta", H5P_DEFAULT), 0);

    H5Fclose(h5);
    boost::filesystem::remove_all(testDir);
  }

  // --- XDMF write-and-close full workflow ------------------------------------

  TEST_P(MPIMeshHDF5, XDMFWriteAndClose)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_mpi_xdmf_wc_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "output";

    auto localMesh = makeMesh(type);
    P1 fes(localMesh);
    GridFunction gf(fes);
    gf.setName("temperature");
    gf = [](const Geometry::Point& p) { return p.x() + 1.0; };

    {
      XDMF xdmf(stem);
      xdmf.setMesh(localMesh);
      xdmf.add("temperature", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.write(1.0);
      xdmf.close();

      EXPECT_TRUE(xdmf.isClosed());
      EXPECT_EQ(xdmf.getSnapshotCount(), 2u);
    }

    // Verify XDMF XML content
    const auto xdmfFile = stem.string() + ".xdmf";
    std::ifstream ifs(xdmfFile);
    ASSERT_TRUE(ifs.good());
    std::ostringstream buffer;
    buffer << ifs.rdbuf();
    const auto text = buffer.str();
    EXPECT_NE(text.find("Xdmf"), std::string::npos);
    EXPECT_NE(text.find("Topology"), std::string::npos);
    EXPECT_NE(text.find("Geometry"), std::string::npos);
    EXPECT_NE(text.find("/Mesh/XDMF/Topology"), std::string::npos);
    EXPECT_NE(text.find("temperature"), std::string::npos);

    boost::filesystem::remove_all(testDir);
  }

  struct MPIPolytopeNameGenerator
  {
    std::string operator()(const ::testing::TestParamInfo<Polytope::Type>& info) const
    {
      return polytopeLabel(info.param);
    }
  };

  INSTANTIATE_TEST_SUITE_P(
      AllDimensions,
      MPIMeshHDF5,
      ::testing::Values(
          Polytope::Type::Segment,
          Polytope::Type::Triangle,
          Polytope::Type::Quadrilateral,
          Polytope::Type::Tetrahedron,
          Polytope::Type::Hexahedron,
          Polytope::Type::Wedge
      ),
      MPIPolytopeNameGenerator());
}
