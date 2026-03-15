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

    hid_t xdmfTopology = H5Dopen2(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT);
    ASSERT_GE(xdmfTopology, 0);
    H5Dclose(xdmfTopology);
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
}
