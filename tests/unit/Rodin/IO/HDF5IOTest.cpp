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
#include <Rodin/Geometry/ParametricTransformation.h>
#include <Rodin/Geometry/Shard.h>
#include <Rodin/Geometry/Sharder.h>
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

  TEST(Rodin_IO_HDF5, XDMFQuadraticTriangleVisualizationTopology)
  {
    const boost::filesystem::path testDir = "/tmp/rodin_xdmf_p2_triangle_test";
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    Mesh mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    RealH1Element<2> fe(Polytope::Type::Triangle);
    PointCloud pm(2, fe.getCount());
    for (size_t i = 0; i < fe.getCount(); ++i)
    {
      const auto& rc = fe.getNode(i);
      pm(0, i) = rc[0];
      pm(1, i) = rc[1];
    }
    mesh.setPolytopeTransformation(
        {2, 0},
        new ParametricTransformation<RealH1Element<2>>(pm, fe));

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

    hid_t vdset = H5Dopen2(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT);
    ASSERT_GE(vdset, 0);
    hid_t vspace = H5Dget_space(vdset);
    hsize_t vdims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(vspace, vdims, nullptr), 2);
    EXPECT_EQ(static_cast<size_t>(vdims[0]), 6u);
    EXPECT_EQ(static_cast<size_t>(vdims[1]), 2u);
    H5Sclose(vspace);
    H5Dclose(vdset);

    hid_t tdset = H5Dopen2(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT);
    ASSERT_GE(tdset, 0);
    hid_t tspace = H5Dget_space(tdset);
    hsize_t tdims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(tspace, tdims, nullptr), 2);
    ASSERT_EQ(static_cast<size_t>(tdims[0]), 1u);
    ASSERT_EQ(static_cast<size_t>(tdims[1]), 6u);
    std::vector<unsigned long long> topology(6);
    ASSERT_GE(H5Dread(
          tdset,
          H5T_NATIVE_ULLONG,
          H5S_ALL,
          H5S_ALL,
          H5P_DEFAULT,
          topology.data()), 0);
    for (size_t i = 0; i < topology.size(); ++i)
      EXPECT_EQ(topology[i], static_cast<unsigned long long>(i));
    H5Sclose(tspace);
    H5Dclose(tdset);
    H5Fclose(h5);

    const auto xdmfFile = testDir / "vis.xdmf";
    ASSERT_TRUE(boost::filesystem::exists(xdmfFile));
    std::ifstream in(xdmfFile.string());
    ASSERT_TRUE(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    const auto xml = buffer.str();
    EXPECT_NE(xml.find("TopologyType=\"Triangle_6\""), std::string::npos);
    EXPECT_NE(xml.find("Dimensions=\"1 6\""), std::string::npos);
    EXPECT_NE(xml.find("Dimensions=\"6 2\""), std::string::npos);

    boost::filesystem::remove_all(testDir);
  }

  static Geometry::Mesh<Context::Local> makeQuadraticXDMFBaseMesh(Polytope::Type type)
  {
    using LocalMesh = Geometry::Mesh<Context::Local>;
    switch (type)
    {
      case Polytope::Type::Segment:
        return LocalMesh::UniformGrid(type, { 3 });
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return LocalMesh::UniformGrid(type, { 2, 2 });
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
        return LocalMesh::UniformGrid(type, { 2, 2, 2 });
      case Polytope::Type::Point:
        break;
    }

    Alert::Exception()
      << "Unsupported quadratic XDMF test polytope."
      << Alert::Raise;
    return {};
  }

  static std::string quadraticXDMFLabel(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Segment:       return "Segment1D";
      case Polytope::Type::Triangle:      return "Triangle2D";
      case Polytope::Type::Quadrilateral: return "Quad2D";
      case Polytope::Type::Tetrahedron:   return "Tet3D";
      case Polytope::Type::Pyramid:       return "Pyramid3D";
      case Polytope::Type::Wedge:         return "Wedge3D";
      case Polytope::Type::Hexahedron:    return "Hex3D";
      case Polytope::Type::Point:         return "Point0D";
    }
    return "Unknown";
  }

  static void setQuadraticCellTransformations(Geometry::Mesh<Context::Local>& mesh)
  {
    const size_t D = mesh.getDimension();
    for (Index i = 0; i < static_cast<Index>(mesh.getCellCount()); ++i)
    {
      const auto geometry = mesh.getConnectivity().getGeometry(D, i);
      RealH1Element<2> fe(geometry);
      PointCloud pm(mesh.getSpaceDimension(), fe.getCount());
      const auto& base = mesh.getPolytopeTransformation(D, i);
      for (size_t j = 0; j < fe.getCount(); ++j)
      {
        Math::SpatialPoint pc;
        base.transform(pc, fe.getNode(j));
        for (size_t d = 0; d < mesh.getSpaceDimension(); ++d)
          pm(d, j) = pc(static_cast<Eigen::Index>(d));
      }
      mesh.setPolytopeTransformation(
          {D, i},
          new ParametricTransformation<RealH1Element<2>>(pm, fe));
    }
  }

  struct QuadraticXDMFCase
  {
    Polytope::Type type;
    unsigned long long topologyId;
    size_t nodesPerCell;
  };

  class Rodin_IO_HDF5_QuadraticXDMF
    : public testing::TestWithParam<QuadraticXDMFCase>
  {};

  TEST_P(Rodin_IO_HDF5_QuadraticXDMF, SupportedQuadraticTopology)
  {
    const auto c = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_xdmf_p2_" + quadraticXDMFLabel(c.type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "vis";

    auto mesh = makeQuadraticXDMFBaseMesh(c.type);
    setQuadraticCellTransformations(mesh);

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

    hid_t vdset = H5Dopen2(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT);
    ASSERT_GE(vdset, 0);
    hid_t vspace = H5Dget_space(vdset);
    hsize_t vdims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(vspace, vdims, nullptr), 2);
    EXPECT_EQ(static_cast<size_t>(vdims[0]), mesh.getCellCount() * c.nodesPerCell);
    EXPECT_EQ(static_cast<size_t>(vdims[1]), mesh.getSpaceDimension());
    H5Sclose(vspace);
    H5Dclose(vdset);

    hid_t tdset = H5Dopen2(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT);
    ASSERT_GE(tdset, 0);
    hid_t tspace = H5Dget_space(tdset);
    hsize_t tdims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(tspace, tdims, nullptr), 2);
    EXPECT_EQ(static_cast<size_t>(tdims[0]), mesh.getCellCount());
    EXPECT_EQ(static_cast<size_t>(tdims[1]), c.nodesPerCell);

    std::vector<unsigned long long> topology(
        static_cast<size_t>(tdims[0] * tdims[1]));
    ASSERT_GE(H5Dread(
          tdset,
          H5T_NATIVE_ULLONG,
          H5S_ALL,
          H5S_ALL,
          H5P_DEFAULT,
          topology.data()), 0);
    ASSERT_FALSE(topology.empty());
    for (size_t i = 0; i < topology.size(); ++i)
      EXPECT_EQ(topology[i], static_cast<unsigned long long>(i));
    H5Sclose(tspace);
    H5Dclose(tdset);
    H5Fclose(h5);

    boost::filesystem::remove_all(testDir);
  }

  TEST_P(Rodin_IO_HDF5_QuadraticXDMF, ForcedMixedQuadraticTopologyUsesXDMFIds)
  {
    const auto c = GetParam();
    const boost::filesystem::path meshH5 =
        "/tmp/rodin_xdmf_p2_mixed_" + quadraticXDMFLabel(c.type) + ".h5";

    auto mesh = makeQuadraticXDMFBaseMesh(c.type);
    setQuadraticCellTransformations(mesh);

    HDF5::writeXDMFMesh(meshH5, mesh, false);

    hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    hid_t tdset = H5Dopen2(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT);
    ASSERT_GE(tdset, 0);
    hid_t tspace = H5Dget_space(tdset);
    ASSERT_EQ(H5Sget_simple_extent_ndims(tspace), 1);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(tspace, &count, nullptr), 1);

    std::vector<unsigned long long> topology(static_cast<size_t>(count));
    ASSERT_GE(H5Dread(
          tdset,
          H5T_NATIVE_ULLONG,
          H5S_ALL,
          H5S_ALL,
          H5P_DEFAULT,
          topology.data()), 0);

    ASSERT_EQ(
        topology.size(),
        mesh.getCellCount() * (1 + c.nodesPerCell));
    for (size_t cell = 0; cell < mesh.getCellCount(); ++cell)
    {
      const size_t offset = cell * (1 + c.nodesPerCell);
      EXPECT_EQ(topology[offset], c.topologyId);
      for (size_t i = 0; i < c.nodesPerCell; ++i)
      {
        EXPECT_EQ(
            topology[offset + 1 + i],
            static_cast<unsigned long long>(cell * c.nodesPerCell + i));
      }
    }

    H5Sclose(tspace);
    H5Dclose(tdset);
    H5Fclose(h5);
    boost::filesystem::remove(meshH5);
  }

  TEST(Rodin_IO_HDF5, XDMFMixedLinearAndCurvedCellsUsesMixedTopology)
  {
    const boost::filesystem::path meshH5 =
        "/tmp/rodin_xdmf_mixed_linear_curved.h5";

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

    RealH1Element<2> fe(Polytope::Type::Triangle);
    PointCloud pm(2, fe.getCount());
    const auto& base = mesh.getPolytopeTransformation(2, 0);
    for (size_t i = 0; i < fe.getCount(); ++i)
    {
      Math::SpatialPoint pc;
      base.transform(pc, fe.getNode(i));
      pm(0, i) = pc(0);
      pm(1, i) = pc(1);
    }
    mesh.setPolytopeTransformation(
        {2, 0},
        new ParametricTransformation<RealH1Element<2>>(pm, fe));

    HDF5::writeXDMFMesh(meshH5, mesh, true);

    hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t tdset = H5Dopen2(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT);
    ASSERT_GE(tdset, 0);
    hid_t tspace = H5Dget_space(tdset);
    ASSERT_EQ(H5Sget_simple_extent_ndims(tspace), 1);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(tspace, &count, nullptr), 1);
    ASSERT_EQ(static_cast<size_t>(count), HDF5::getXDMFMixedTopologySize(mesh));

    std::vector<unsigned long long> topology(static_cast<size_t>(count));
    ASSERT_GE(H5Dread(
          tdset,
          H5T_NATIVE_ULLONG,
          H5S_ALL,
          H5S_ALL,
          H5P_DEFAULT,
          topology.data()), 0);
    ASSERT_GE(topology.size(), 11u);
    EXPECT_EQ(topology[0], HDF5::getXDMFQuadraticMixedTopologyId(Polytope::Type::Triangle));
    EXPECT_EQ(topology[7], HDF5::getXDMFMixedTopologyId(Polytope::Type::Triangle));

    H5Sclose(tspace);
    H5Dclose(tdset);
    H5Fclose(h5);
    boost::filesystem::remove(meshH5);
  }

  TEST(Rodin_IO_HDF5, XDMFCurvedAttributesEvaluateOnVisualizationPoints)
  {
    const boost::filesystem::path nodeH5 =
        "/tmp/rodin_xdmf_p2_attr_node.h5";
    const boost::filesystem::path cellH5 =
        "/tmp/rodin_xdmf_p2_attr_cell.h5";

    Mesh mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({2, 0})
      .vertex({0, 3})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();
    setQuadraticCellTransformations(mesh);

    P1 fes(mesh);
    GridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return 2.0 * p.x() + 3.0 * p.y(); };

    HDF5::writeXDMFNodeAttribute(gf, mesh, nodeH5);
    HDF5::writeXDMFCellAttribute(gf, mesh, cellH5);

    hid_t h5 = H5Fopen(nodeH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    hid_t dspace = H5Dget_space(dset);
    hsize_t count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(dspace, &count, nullptr), 1);
    ASSERT_EQ(static_cast<size_t>(count), 6u);
    std::vector<double> nodeValues(static_cast<size_t>(count));
    ASSERT_GE(H5Dread(
          dset,
          H5T_NATIVE_DOUBLE,
          H5S_ALL,
          H5S_ALL,
          H5P_DEFAULT,
          nodeValues.data()), 0);
    const auto nodes = HDF5::getXDMFReferenceNodes(Polytope::Type::Triangle, 2);
    const auto& trans = mesh.getPolytopeTransformation(2, 0);
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      Math::SpatialPoint pc;
      trans.transform(pc, nodes[i]);
      EXPECT_NEAR(nodeValues[i], 2.0 * pc(0) + 3.0 * pc(1), 1e-14);
    }
    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);

    h5 = H5Fopen(cellH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);
    dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
    ASSERT_GE(dset, 0);
    dspace = H5Dget_space(dset);
    count = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(dspace, &count, nullptr), 1);
    ASSERT_EQ(static_cast<size_t>(count), 1u);
    double cellValue = 0.0;
    ASSERT_GE(H5Dread(
          dset,
          H5T_NATIVE_DOUBLE,
          H5S_ALL,
          H5S_ALL,
          H5P_DEFAULT,
          &cellValue), 0);
    EXPECT_NEAR(cellValue, 13.0 / 3.0, 1e-14);
    H5Sclose(dspace);
    H5Dclose(dset);
    H5Fclose(h5);

    boost::filesystem::remove(nodeH5);
    boost::filesystem::remove(cellH5);
  }

  TEST(Rodin_IO_HDF5, XDMFShardWritesOwnedCurvedCellsOnly)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(2, 2);
    mesh.getConnectivity().compute(2, 0);
    setQuadraticCellTransformations(mesh);

    BalancedCompactPartitioner partitioner(mesh);
    partitioner.partition(2);
    SharderBase<Context::Local> sharder(Context::Local{});
    sharder.shard(partitioner);
    const auto& shard = sharder.getShards().front();

    const boost::filesystem::path meshH5 =
        "/tmp/rodin_xdmf_shard_owned_curved.h5";
    HDF5::writeXDMFMesh(meshH5, shard, false);

    hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(h5, 0);

    hid_t vdset = H5Dopen2(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT);
    ASSERT_GE(vdset, 0);
    hid_t vspace = H5Dget_space(vdset);
    hsize_t vdims[2] = {0, 0};
    ASSERT_EQ(H5Sget_simple_extent_dims(vspace, vdims, nullptr), 2);
    EXPECT_EQ(
        static_cast<size_t>(vdims[0]),
        HDF5::getXDMFVisualizationVertexCount(shard));
    EXPECT_EQ(static_cast<size_t>(vdims[1]), shard.getSpaceDimension());
    H5Sclose(vspace);
    H5Dclose(vdset);

    hid_t tdset = H5Dopen2(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT);
    ASSERT_GE(tdset, 0);
    hid_t tspace = H5Dget_space(tdset);
    hsize_t tcount = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(tspace, &tcount, nullptr), 1);
    EXPECT_EQ(
        static_cast<size_t>(tcount),
        HDF5::getXDMFMixedTopologySize(shard));
    H5Sclose(tspace);
    H5Dclose(tdset);

    hid_t adset = H5Dopen2(h5, HDF5::attributePath(shard.getDimension()).c_str(), H5P_DEFAULT);
    ASSERT_GE(adset, 0);
    hid_t aspace = H5Dget_space(adset);
    hsize_t acount = 0;
    ASSERT_EQ(H5Sget_simple_extent_dims(aspace, &acount, nullptr), 1);
    EXPECT_EQ(
        static_cast<size_t>(acount),
        HDF5::getXDMFRenderedCellCount(shard));
    H5Sclose(aspace);
    H5Dclose(adset);

    H5Fclose(h5);
    boost::filesystem::remove(meshH5);
  }

  INSTANTIATE_TEST_SUITE_P(
      Rodin_IO_HDF5,
      Rodin_IO_HDF5_QuadraticXDMF,
      testing::Values(
        QuadraticXDMFCase{Polytope::Type::Segment,       34u, 3u},
        QuadraticXDMFCase{Polytope::Type::Triangle,      36u, 6u},
        QuadraticXDMFCase{Polytope::Type::Quadrilateral, 35u, 9u},
        QuadraticXDMFCase{Polytope::Type::Tetrahedron,   38u, 10u},
        QuadraticXDMFCase{Polytope::Type::Pyramid,       39u, 13u},
        QuadraticXDMFCase{Polytope::Type::Wedge,         41u, 18u},
        QuadraticXDMFCase{Polytope::Type::Hexahedron,    50u, 27u}));

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
      case Polytope::Type::Point:
      {
        Mesh<Context::Local>::Builder builder;
        return builder
          .initialize(1)
          .nodes(1)
          .vertex({ 0.0 })
          .finalize();
      }
      case Polytope::Type::Segment:
        return LocalMesh::UniformGrid(type, { 11 });
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return LocalMesh::UniformGrid(type, { 4, 4 });
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
        return LocalMesh::UniformGrid(type, { 3, 3, 3 });
      default:
        ADD_FAILURE() << "Unsupported polytope type for makeMesh";
        return Geometry::Mesh<Context::Local>();
    }
  }

  Geometry::Mesh<Context::Local> makeAttributeRegressionMesh(size_t dim)
  {
    Mesh<Context::Local>::Builder builder;
    builder.initialize(dim);

    switch (dim)
    {
      case 1:
      {
        builder
          .nodes(3)
          .vertex({ 0.0 })
          .vertex({ 1.0 })
          .vertex({ 2.0 })
          .polytope(Polytope::Type::Segment, { 0, 1 })
          .attribute({ 1, 0 }, 11)
          .polytope(Polytope::Type::Segment, { 1, 2 })
          .attribute({ 1, 1 }, 22);
        break;
      }
      case 2:
      {
        builder
          .nodes(5)
          .vertex({ 0.0, 0.0 })
          .vertex({ 1.0, 0.0 })
          .vertex({ 0.0, 1.0 })
          .vertex({ 2.0, 0.0 })
          .vertex({ 2.0, 1.0 })
          .polytope(Polytope::Type::Triangle, { 0, 1, 2 })
          .attribute({ 2, 0 }, 11)
          .polytope(Polytope::Type::Quadrilateral, { 1, 3, 4, 2 })
          .attribute({ 2, 1 }, 22);
        break;
      }
      case 3:
      {
        builder
          .nodes(8)
          .vertex({ 0.0, 0.0, 0.0 })
          .vertex({ 1.0, 0.0, 0.0 })
          .vertex({ 1.0, 1.0, 0.0 })
          .vertex({ 0.0, 1.0, 0.0 })
          .vertex({ 0.0, 0.0, 1.0 })
          .vertex({ 1.0, 0.0, 1.0 })
          .vertex({ 1.0, 1.0, 1.0 })
          .vertex({ 0.0, 1.0, 1.0 })
          .polytope(Polytope::Type::Tetrahedron, { 0, 1, 2, 4 })
          .attribute({ 3, 0 }, 11)
          .polytope(Polytope::Type::Hexahedron, { 0, 1, 2, 3, 4, 5, 6, 7 })
          .attribute({ 3, 1 }, 22);
        break;
      }
      default:
        ADD_FAILURE() << "Unsupported dimension for attribute regression mesh";
    }
    return builder.finalize();
  }

  void expectAttributeRegressionMesh(
      const Geometry::Mesh<Context::Local>& mesh,
      size_t dim)
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

  // Helper to return a human-readable label for test names
  std::string polytopeLabel(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Point:         return "Point0D";
      case Polytope::Type::Segment:       return "Segment1D";
      case Polytope::Type::Triangle:      return "Triangle2D";
      case Polytope::Type::Quadrilateral: return "Quad2D";
      case Polytope::Type::Tetrahedron:   return "Tet3D";
      case Polytope::Type::Pyramid:       return "Pyramid3D";
      case Polytope::Type::Hexahedron:    return "Hex3D";
      case Polytope::Type::Wedge:         return "Wedge3D";
      default:                            return "Unknown";
    }
  }
  // Parameterized test fixture
  class HDF5MultiDim : public ::testing::TestWithParam<Polytope::Type> {};
  class HDF5AttributeRegression : public ::testing::TestWithParam<size_t> {};

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

  TEST_P(HDF5AttributeRegression, MeshRoundTripPreservesDimensionAttributes)
  {
    const size_t dim = GetParam();
    const std::string meshFile =
        "/tmp/rodin_hdf5_attr_regression_" + std::to_string(dim) + "d.h5";

    Mesh mesh = makeAttributeRegressionMesh(dim);
    mesh.save(meshFile, FileFormat::HDF5);

    Mesh loaded;
    loaded.load(meshFile, FileFormat::HDF5);

    expectAttributeRegressionMesh(loaded, dim);

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
    if (type == Polytope::Type::Point)
    {
      P0 fes(mesh);
      GridFunction gf(fes);
      gf[0] = 1.0;

      gf.save(gfFile, FileFormat::HDF5);

      hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(h5, 0);
      hid_t dset = H5Dopen2(h5, "/GridFunction/Values/Data", H5P_DEFAULT);
      ASSERT_GE(dset, 0);
      hid_t dspace = H5Dget_space(dset);
      hsize_t count = 0;
      H5Sget_simple_extent_dims(dspace, &count, nullptr);
      EXPECT_EQ(static_cast<size_t>(count), static_cast<size_t>(gf.getData().size()));
      H5Sclose(dspace);
      H5Dclose(dset);
      H5Fclose(h5);

      GridFunction loaded(fes);
      loaded.load(gfFile, FileFormat::HDF5);
      ASSERT_EQ(loaded.getData().size(), gf.getData().size());
      EXPECT_DOUBLE_EQ(loaded.getData()[0], gf.getData()[0]);

      std::remove(gfFile.c_str());
      return;
    }

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
    if (type == Polytope::Type::Point)
    {
      P0 fes(mesh);
      GridFunction gf(fes);
      gf[0] = 1.0;

      gf.save(gfFile, FileFormat::HDF5);

      hid_t h5 = H5Fopen(gfFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(h5, 0);
      hid_t metaSize = H5Dopen2(h5, "/GridFunction/Meta/Size", H5P_DEFAULT);
      ASSERT_GE(metaSize, 0);
      H5Dclose(metaSize);
      hid_t metaDim = H5Dopen2(h5, "/GridFunction/Meta/Dimension", H5P_DEFAULT);
      ASSERT_GE(metaDim, 0);
      H5Dclose(metaDim);
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
      return;
    }

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
    EXPECT_EQ(static_cast<size_t>(vdims[0]), HDF5::getXDMFVisualizationVertexCount(mesh));
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
    if (type == Polytope::Type::Point)
    {
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
      EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT), 1);
      H5Fclose(h5);
      boost::filesystem::remove_all(testDir);
      return;
    }

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

    // For Node-centered XDMF data, size follows visualization nodes. Curved
    // visualization samples order > 1 transformations at quadratic XDMF nodes.
    EXPECT_EQ(static_cast<size_t>(count), HDF5::getXDMFVisualizationVertexCount(mesh));

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
    if (type == Polytope::Type::Point)
    {
      {
        XDMF xdmf(stem);
        xdmf.setMesh(mesh);
        xdmf.write(0.0);
        xdmf.write(1.0);
        xdmf.close();

        EXPECT_TRUE(xdmf.isClosed());
        EXPECT_EQ(xdmf.getSnapshotCount(), 2u);
        EXPECT_EQ(xdmf.getGridCount(), 1u);
      }

      const auto xdmfFile = stem.string() + ".xdmf";
      std::ifstream ifs(xdmfFile);
      ASSERT_TRUE(ifs.good());
      std::ostringstream buffer;
      buffer << ifs.rdbuf();
      const auto text = buffer.str();
      EXPECT_NE(text.find("Xdmf"), std::string::npos);
      EXPECT_NE(text.find("Topology"), std::string::npos);
      EXPECT_NE(text.find("/Mesh/XDMF/Topology"), std::string::npos);

      boost::filesystem::remove_all(testDir);
      return;
    }

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
          Polytope::Type::Point,         // 0D
          Polytope::Type::Triangle,      // 2D
          Polytope::Type::Quadrilateral, // 2D
          Polytope::Type::Tetrahedron,   // 3D
          Polytope::Type::Pyramid,       // 3D
          Polytope::Type::Hexahedron,    // 3D
          Polytope::Type::Wedge          // 3D
      ),
      PolytopeNameGenerator());

  INSTANTIATE_TEST_SUITE_P(
      Dimensions,
      HDF5AttributeRegression,
      ::testing::Values(1, 2, 3),
      [](const ::testing::TestParamInfo<size_t>& info)
      {
        return "Dim" + std::to_string(info.param) + "D";
      });
}
