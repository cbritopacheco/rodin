/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
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

#include <hdf5.h>

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_IO_HDF5, MeshRoundTrip)
  {
    const std::string meshFile = "/tmp/rodin_hdf5_mesh_rt.h5";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 2);

    mesh.save(meshFile, FileFormat::HDF5);

    // Verify HDF5 dataset presence
    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t vertices = H5Dopen2(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT);
    ASSERT_GE(vertices, 0);
    hid_t vspace = H5Dget_space(vertices);
    ASSERT_GE(vspace, 0);
    int rank = H5Sget_simple_extent_ndims(vspace);
    ASSERT_EQ(rank, 2);
    hsize_t dims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(vspace, dims, nullptr), 2);
    EXPECT_EQ(static_cast<size_t>(dims[0]), mesh.getVertexCount());
    EXPECT_EQ(static_cast<size_t>(dims[1]), mesh.getSpaceDimension());
    H5Sclose(vspace);
    H5Dclose(vertices);
    H5Fclose(h5);

    // Round-trip: load and compare
    Mesh loaded;
    loaded.load(meshFile, FileFormat::HDF5);
    EXPECT_EQ(loaded.getSpaceDimension(), mesh.getSpaceDimension());
    EXPECT_EQ(loaded.getDimension(), mesh.getDimension());
    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getPolytopeCount(0), mesh.getPolytopeCount(0));
    EXPECT_EQ(
        loaded.getPolytopeCount(mesh.getDimension()),
        mesh.getPolytopeCount(mesh.getDimension()));

    std::remove(meshFile.c_str());
  }

  TEST(Rodin_IO_HDF5, MeshPersistenceNoXDMF)
  {
    // The HDF5 persistence path (mesh.save) must NOT write /Mesh/XDMF datasets.
    // Those are only written by the XDMF visualization pipeline.
    const std::string meshFile = "/tmp/rodin_hdf5_no_xdmf.h5";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.save(meshFile, FileFormat::HDF5);

    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Canonical mesh data should be present
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Meta/SpaceDimension", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 1);

    // XDMF-specific datasets must NOT be present
    EXPECT_EQ(H5Lexists(h5, "/Mesh/XDMF", H5P_DEFAULT), 0);

    H5Fclose(h5);

    // Verify the canonical persistence data can be loaded back
    Mesh loaded;
    loaded.load(meshFile, FileFormat::HDF5);
    EXPECT_EQ(loaded.getSpaceDimension(), mesh.getSpaceDimension());
    EXPECT_EQ(loaded.getDimension(), mesh.getDimension());
    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());

    std::remove(meshFile.c_str());
  }

  TEST(Rodin_IO_HDF5, GridFunctionPersistenceRawDOFs)
  {
    // The HDF5 persistence path (gf.save) must write raw DOF vector data,
    // not evaluated vertex/cell values. This validates that the persistence
    // path is separate from the visualization path.
    const std::string gfFile = "/tmp/rodin_hdf5_gf_rawdofs.h5";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() * p.x() + p.y(); };

    gf.save(gfFile, FileFormat::HDF5);

    // Read back and compare with the raw DOF vector
    hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    hid_t dspace = H5Dget_space(dset);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);

    // The dataset size must match the raw DOF vector, not the vertex count
    EXPECT_EQ(static_cast<size_t>(count), static_cast<size_t>(gf.getData().size()));

    // Read the actual values and compare with DOFs
    std::vector<double> saved(static_cast<size_t>(count));
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, saved.data());

    const auto& dofs = gf.getData();
    for (size_t i = 0; i < saved.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(saved[i], static_cast<double>(dofs[static_cast<Eigen::Index>(i)]))
        << "DOF mismatch at index " << i;
    }

    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);

    // Also verify round-trip through GridFunction load
    GridFunction loaded(fes);
    loaded.load(gfFile, FileFormat::HDF5);
    EXPECT_EQ(loaded.getData().size(), gf.getData().size());
    for (Eigen::Index i = 0; i < gf.getData().size(); ++i)
    {
      EXPECT_DOUBLE_EQ(loaded.getData()[i], gf.getData()[i])
        << "Round-trip DOF mismatch at index " << i;
    }

    std::remove(gfFile.c_str());
  }

  TEST(Rodin_IO_HDF5, XDMFVisualizationTopology)
  {
    // The XDMF visualization path must write only minimal visualization data
    // to the mesh HDF5 file — not the full canonical Rodin persistence.
    const boost::filesystem::path testDir = "/tmp/rodin_xdmf_topo_test";
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    {
      XDMF xdmf(stem);
      xdmf.setMesh(mesh);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Find the mesh HDF5 file written by the XDMF pipeline
    const auto meshH5 = testDir / "vis.mesh.h5";
    ASSERT_TRUE(boost::filesystem::exists(meshH5));

    hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Visualization-only data should be present
    EXPECT_GE(H5Lexists(h5, "/Mesh", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/TopologySize", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Attributes", H5P_DEFAULT), 1);

    // Canonical persistence data must NOT be present — XDMF files are
    // visualization-only and do not contain full Rodin persistence.
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 0);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Transformations", H5P_DEFAULT), 0);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Meta", H5P_DEFAULT), 0);

    H5Fclose(h5);
    boost::filesystem::remove_all(testDir);
  }

  TEST(Rodin_IO_HDF5, XDMFVisualizationEvaluatedAttributes)
  {
    // The XDMF visualization path must write evaluated vertex data,
    // not raw DOFs. This validates the separation from the persistence path.
    const boost::filesystem::path testDir = "/tmp/rodin_xdmf_eval_test";
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    gf.setName("field");
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    {
      XDMF xdmf(stem);
      xdmf.setMesh(mesh);
      xdmf.add("field", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Find an attribute HDF5 file written by the XDMF pipeline
    boost::filesystem::path attrH5;
    for (auto& entry : boost::filesystem::directory_iterator(testDir))
    {
      const auto fn = entry.path().filename().string();
      if (fn.find("field") != std::string::npos && fn.find(".h5") != std::string::npos)
      {
        attrH5 = entry.path();
        break;
      }
    }
    ASSERT_FALSE(attrH5.empty()) << "No attribute HDF5 file found";

    hid_t h5 = H5Fopen(attrH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    hid_t dspace = H5Dget_space(dset);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);

    // For Node-centered data, the XDMF path evaluates at vertices,
    // so the dataset size should match the vertex count (not the DOF count).
    EXPECT_EQ(static_cast<size_t>(count), mesh.getVertexCount());

    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);
    boost::filesystem::remove_all(testDir);
  }

  TEST(Rodin_IO_HDF5, GridFunctionStandalone)
  {
    const std::string gfFile = "/tmp/rodin_hdf5_gf_sa.h5";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    gf.save(gfFile, FileFormat::HDF5);

    hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t metaSize = H5Dopen2(h5, "/GridFunction/Meta/Size", H5P_DEFAULT);
    ASSERT_GE(metaSize, 0);
    H5Dclose(metaSize);
    hid_t metaDimension = H5Dopen2(h5, "/GridFunction/Meta/Dimension", H5P_DEFAULT);
    ASSERT_GE(metaDimension, 0);
    H5Dclose(metaDimension);
    hid_t values = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(values, 0);
    hid_t dspace = H5Dget_space(values);
    ASSERT_GE(dspace, 0);
    int rank = H5Sget_simple_extent_ndims(dspace);
    ASSERT_EQ(rank, 1);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(dspace, &count, nullptr), 1);
    EXPECT_EQ(static_cast<size_t>(count), static_cast<size_t>(gf.getData().size()));
    H5Sclose(dspace);
    H5Dclose(values);
    H5Fclose(h5);

    std::remove(gfFile.c_str());
  }

  TEST(Rodin_IO_HDF5, XDMFWriteAndClose)
  {
    // Use a dedicated subdirectory for test output
    const boost::filesystem::path testDir = "/tmp/rodin_xdmf_test_dir";
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "output";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(2, 2);

    P1 fes(mesh);
    GridFunction gf(fes);
    gf.setName("temperature");
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    {
      XDMF xdmf(stem);
      xdmf.setMesh(mesh);
      xdmf.add("temperature", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.write(1.0);
      xdmf.close();

      EXPECT_TRUE(xdmf.isClosed());
      EXPECT_EQ(xdmf.getSnapshotCount(), 2u);
      EXPECT_EQ(xdmf.getGridCount(), 1u);
    }

    // Verify the XDMF XML file was written
    const auto xdmfFile = stem.string() + ".xdmf";
    std::ifstream ifs(xdmfFile);
    ASSERT_TRUE(ifs.good());
    std::ostringstream buffer;
    buffer << ifs.rdbuf();
    const auto text = buffer.str();
    EXPECT_NE(text.find("Xdmf"), std::string::npos);
    EXPECT_NE(text.find("Domain"), std::string::npos);
    EXPECT_NE(text.find("Topology"), std::string::npos);
    EXPECT_NE(text.find("Geometry"), std::string::npos);
    EXPECT_NE(text.find("/Mesh/XDMF/Topology"), std::string::npos);
    EXPECT_NE(text.find("/Mesh/Geometry/Vertices"), std::string::npos);
    EXPECT_NE(text.find("/GridFunction/Values/Data"), std::string::npos);
    EXPECT_NE(text.find("temperature"), std::string::npos);

    // Clean up the dedicated test directory
    boost::filesystem::remove_all(testDir);
  }

  // ===========================================================================
  // Parameterized tests: 1D, 2D, and 3D across all HDF5/XDMF paths
  // ===========================================================================

} // namespace Rodin::Tests::Unit

namespace
{
  using namespace Rodin;
  using namespace Rodin::IO;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  // Helper to create a mesh for each polytope type with appropriate grid shape
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

  // Helper to return a human-readable label for test names
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
  class HDF5MultiDim : public ::testing::TestWithParam<Polytope::Type> {};

  // --- Mesh round-trip persistence across dimensions -------------------------

  TEST_P(HDF5MultiDim, MeshRoundTrip)
  {
    const auto type = GetParam();
    const std::string meshFile =
        "/tmp/rodin_hdf5_rt_" + polytopeLabel(type) + ".h5";

    Mesh mesh = makeMesh(type);
    ASSERT_GT(mesh.getVertexCount(), 0u);
    ASSERT_GT(mesh.getCellCount(), 0u);

    mesh.save(meshFile, FileFormat::HDF5);

    // Verify vertex dataset presence and shape
    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t vertices = H5Dopen2(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT);
    ASSERT_GE(vertices, 0);
    hid_t vspace = H5Dget_space(vertices);
    ASSERT_GE(vspace, 0);
    hsize_t dims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(vspace, dims, nullptr), 2);
    EXPECT_EQ(static_cast<size_t>(dims[0]), mesh.getVertexCount());
    EXPECT_EQ(static_cast<size_t>(dims[1]), mesh.getSpaceDimension());
    H5Sclose(vspace);
    H5Dclose(vertices);
    H5Fclose(h5);

    // Load and compare
    Mesh loaded;
    loaded.load(meshFile, FileFormat::HDF5);
    EXPECT_EQ(loaded.getSpaceDimension(), mesh.getSpaceDimension());
    EXPECT_EQ(loaded.getDimension(), mesh.getDimension());
    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getPolytopeCount(0), mesh.getPolytopeCount(0));
    EXPECT_EQ(
        loaded.getPolytopeCount(mesh.getDimension()),
        mesh.getPolytopeCount(mesh.getDimension()));

    std::remove(meshFile.c_str());
  }

  // --- Mesh persistence must NOT contain XDMF datasets -----------------------

  TEST_P(HDF5MultiDim, MeshPersistenceNoXDMF)
  {
    const auto type = GetParam();
    const std::string meshFile =
        "/tmp/rodin_hdf5_noxdmf_" + polytopeLabel(type) + ".h5";

    Mesh mesh = makeMesh(type);
    mesh.save(meshFile, FileFormat::HDF5);

    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Canonical data present
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Meta/SpaceDimension", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 1);

    // XDMF-specific must NOT be present
    EXPECT_EQ(H5Lexists(h5, "/Mesh/XDMF", H5P_DEFAULT), 0);

    H5Fclose(h5);

    // Load-back verification
    Mesh loaded;
    loaded.load(meshFile, FileFormat::HDF5);
    EXPECT_EQ(loaded.getSpaceDimension(), mesh.getSpaceDimension());
    EXPECT_EQ(loaded.getDimension(), mesh.getDimension());
    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());

    std::remove(meshFile.c_str());
  }

  // --- GridFunction persistence writes raw DOFs across dimensions ------------

  TEST_P(HDF5MultiDim, GridFunctionPersistenceRawDOFs)
  {
    const auto type = GetParam();
    const std::string gfFile =
        "/tmp/rodin_hdf5_gf_dofs_" + polytopeLabel(type) + ".h5";

    Mesh mesh = makeMesh(type);
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + 1.0; };

    gf.save(gfFile, FileFormat::HDF5);

    // Verify raw DOF count matches
    hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    hid_t dspace = H5Dget_space(dset);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);
    EXPECT_EQ(static_cast<size_t>(count), static_cast<size_t>(gf.getData().size()));

    // Read values and compare DOFs
    std::vector<double> saved(static_cast<size_t>(count));
    H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, saved.data());
    const auto& dofs = gf.getData();
    for (size_t i = 0; i < saved.size(); ++i)
    {
      EXPECT_DOUBLE_EQ(saved[i], static_cast<double>(dofs[static_cast<Eigen::Index>(i)]))
          << "DOF mismatch at index " << i;
    }

    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);

    // Round-trip through GridFunction load
    GridFunction loaded(fes);
    loaded.load(gfFile, FileFormat::HDF5);
    EXPECT_EQ(loaded.getData().size(), gf.getData().size());
    for (Eigen::Index i = 0; i < gf.getData().size(); ++i)
    {
      EXPECT_DOUBLE_EQ(loaded.getData()[i], gf.getData()[i])
          << "Round-trip DOF mismatch at index " << i;
    }

    std::remove(gfFile.c_str());
  }

  // --- GridFunction standalone field file across dimensions -------------------

  TEST_P(HDF5MultiDim, GridFunctionStandalone)
  {
    const auto type = GetParam();
    const std::string gfFile =
        "/tmp/rodin_hdf5_gf_sa_" + polytopeLabel(type) + ".h5";

    Mesh mesh = makeMesh(type);
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + 1.0; };

    gf.save(gfFile, FileFormat::HDF5);

    hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Meta datasets present
    hid_t metaSize = H5Dopen2(h5, "/GridFunction/Meta/Size", H5P_DEFAULT);
    ASSERT_GE(metaSize, 0);
    H5Dclose(metaSize);
    hid_t metaDim = H5Dopen2(h5, "/GridFunction/Meta/Dimension", H5P_DEFAULT);
    ASSERT_GE(metaDim, 0);
    H5Dclose(metaDim);

    // Values dataset present with correct shape
    hid_t values = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(values, 0);
    hid_t dspace = H5Dget_space(values);
    int rank = H5Sget_simple_extent_ndims(dspace);
    ASSERT_EQ(rank, 1);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);
    EXPECT_EQ(static_cast<size_t>(count), static_cast<size_t>(gf.getData().size()));

    H5Sclose(dspace);
    H5Dclose(values);
    H5Fclose(h5);
    std::remove(gfFile.c_str());
  }

  // --- XDMF visualization topology across dimensions -------------------------

  TEST_P(HDF5MultiDim, XDMFVisualizationTopology)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_xdmf_topo_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    Mesh mesh = makeMesh(type);

    {
      XDMF xdmf(stem);
      xdmf.setMesh(mesh);
      xdmf.write(0.0);
      xdmf.close();
    }

    const auto meshH5 = testDir / "vis.mesh.h5";
    ASSERT_TRUE(boost::filesystem::exists(meshH5));

    hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    // Visualization-only data present
    EXPECT_GE(H5Lexists(h5, "/Mesh", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/TopologySize", H5P_DEFAULT), 1);
    EXPECT_GE(H5Lexists(h5, "/Mesh/Attributes", H5P_DEFAULT), 1);

    // Canonical persistence must NOT be present
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 0);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Transformations", H5P_DEFAULT), 0);
    EXPECT_EQ(H5Lexists(h5, "/Mesh/Meta", H5P_DEFAULT), 0);

    // Verify vertex dimensions
    hid_t vdset = H5Dopen2(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT);
    ASSERT_GE(vdset, 0);
    hid_t vspace = H5Dget_space(vdset);
    hsize_t vdims[2] = {0, 0};
    H5Sget_simple_extent_dims(vspace, vdims, nullptr);
    EXPECT_EQ(static_cast<size_t>(vdims[0]), mesh.getVertexCount());
    EXPECT_EQ(static_cast<size_t>(vdims[1]), mesh.getSpaceDimension());
    H5Sclose(vspace);
    H5Dclose(vdset);

    H5Fclose(h5);
    boost::filesystem::remove_all(testDir);
  }

  // --- XDMF visualization evaluated attributes across dimensions -------------

  TEST_P(HDF5MultiDim, XDMFVisualizationEvaluatedAttributes)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_xdmf_eval_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    Mesh mesh = makeMesh(type);
    P1 fes(mesh);
    GridFunction gf(fes);
    gf.setName("field");
    gf = [](const Geometry::Point& p) { return p.x() + 1.0; };

    {
      XDMF xdmf(stem);
      xdmf.setMesh(mesh);
      xdmf.add("field", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Find the attribute HDF5 file
    boost::filesystem::path attrH5;
    for (auto& entry : boost::filesystem::directory_iterator(testDir))
    {
      const auto fn = entry.path().filename().string();
      if (fn.find("field") != std::string::npos && fn.find(".h5") != std::string::npos)
      {
        attrH5 = entry.path();
        break;
      }
    }
    ASSERT_FALSE(attrH5.empty()) << "No attribute HDF5 file found";

    hid_t h5 = H5Fopen(attrH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    hid_t dspace = H5Dget_space(dset);
    hsize_t count = 0;
    H5Sget_simple_extent_dims(dspace, &count, nullptr);

    // For Node-centered data, size = vertex count (not DOF count)
    EXPECT_EQ(static_cast<size_t>(count), mesh.getVertexCount());

    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);
    boost::filesystem::remove_all(testDir);
  }

  // --- XDMF write-and-close across dimensions (full workflow) ----------------

  TEST_P(HDF5MultiDim, XDMFWriteAndClose)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_xdmf_wc_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "output";

    Mesh mesh = makeMesh(type);
    P1 fes(mesh);
    GridFunction gf(fes);
    gf.setName("temperature");
    gf = [](const Geometry::Point& p) { return p.x() + 1.0; };

    {
      XDMF xdmf(stem);
      xdmf.setMesh(mesh);
      xdmf.add("temperature", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.write(1.0);
      xdmf.close();

      EXPECT_TRUE(xdmf.isClosed());
      EXPECT_EQ(xdmf.getSnapshotCount(), 2u);
      EXPECT_EQ(xdmf.getGridCount(), 1u);
    }

    // Verify XDMF XML content
    const auto xdmfFile = stem.string() + ".xdmf";
    std::ifstream ifs(xdmfFile);
    ASSERT_TRUE(ifs.good());
    std::ostringstream buffer;
    buffer << ifs.rdbuf();
    const auto text = buffer.str();
    EXPECT_NE(text.find("Xdmf"), std::string::npos);
    EXPECT_NE(text.find("Domain"), std::string::npos);
    EXPECT_NE(text.find("Topology"), std::string::npos);
    EXPECT_NE(text.find("Geometry"), std::string::npos);
    EXPECT_NE(text.find("/Mesh/XDMF/Topology"), std::string::npos);
    EXPECT_NE(text.find("/Mesh/Geometry/Vertices"), std::string::npos);
    EXPECT_NE(text.find("/GridFunction/Values/Data"), std::string::npos);
    EXPECT_NE(text.find("temperature"), std::string::npos);

    boost::filesystem::remove_all(testDir);
  }

  // --- Custom name generator for readable test output ------------------------
  struct PolytopeNameGenerator
  {
    std::string operator()(const ::testing::TestParamInfo<Polytope::Type>& info) const
    {
      return polytopeLabel(info.param);
    }
  };

  // Instantiate parameterized tests for 1D, 2D, and 3D polytope types
  INSTANTIATE_TEST_SUITE_P(
      AllDimensions,
      HDF5MultiDim,
      ::testing::Values(
          Polytope::Type::Segment,       // 1D
          Polytope::Type::Triangle,      // 2D
          Polytope::Type::Quadrilateral, // 2D
          Polytope::Type::Tetrahedron,   // 3D
          Polytope::Type::Hexahedron,    // 3D
          Polytope::Type::Wedge          // 3D
      ),
      PolytopeNameGenerator());
}
