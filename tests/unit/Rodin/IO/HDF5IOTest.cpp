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

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>

#if defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5
  #include <hdf5.h>
#endif

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_IO_HDF5, SaveMeshGridFunctionAndXDMF)
  {
#if !(defined(RODIN_IO_HAS_HDF5) && RODIN_IO_HAS_HDF5)
    GTEST_SKIP() << "Rodin was built without HDF5 support.";
#else
    const std::string meshFile = "/tmp/rodin_hdf5_mesh.h5";
    const std::string gfFile = "/tmp/rodin_hdf5_gf.h5";
    const std::string xdmfFile = "/tmp/rodin_hdf5.xdmf";

    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    mesh.save(meshFile, FileFormat::HDF5);
    gf.save(gfFile, FileFormat::HDF5);

    hid_t h5 = H5Fopen(meshFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t vertices = H5Dopen2(h5, "/Mesh/Vertices", H5P_DEFAULT);
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

    h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t values = H5Dopen2(h5, "/GridFunction/Values", H5P_DEFAULT);
    ASSERT_GE(values, 0);
    hid_t dspace = H5Dget_space(values);
    ASSERT_GE(dspace, 0);
    rank = H5Sget_simple_extent_ndims(dspace);
    ASSERT_EQ(rank, 1);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(dspace, &count, nullptr), 1);
    EXPECT_EQ(static_cast<size_t>(count), static_cast<size_t>(gf.getData().size()));
    H5Sclose(dspace);
    H5Dclose(values);
    H5Fclose(h5);

    XDMF xdmf(mesh, meshFile);
    xdmf.addGridFunction("u", gf, gfFile);
    xdmf.save(xdmfFile);

    std::ifstream ifs(xdmfFile);
    ASSERT_TRUE(ifs.good());
    std::ostringstream buffer;
    buffer << ifs.rdbuf();
    const auto text = buffer.str();
    EXPECT_NE(text.find(meshFile + ":/Mesh/XDMFTopology"), std::string::npos);
    EXPECT_NE(text.find(meshFile + ":/Mesh/Vertices"), std::string::npos);
    EXPECT_NE(text.find(gfFile + ":/GridFunction/Values"), std::string::npos);

    std::remove(meshFile.c_str());
    std::remove(gfFile.c_str());
    std::remove(xdmfFile.c_str());
#endif
  }
}
