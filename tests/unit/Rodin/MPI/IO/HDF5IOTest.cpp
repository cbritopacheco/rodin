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

  // --- Shard metadata HDF5 dataset structure verification --------------------

  /**
   * Verifies that the `/Shard` group hierarchy written by the MPI printer
   * has the expected layout.  Because MPI is not available in the unit-test
   * environment we create the datasets directly using the HDF5 path helpers
   * and then verify them with raw HDF5 API calls.
   */
  TEST(ShardMetadata, GroupAndDatasetLayout)
  {
    const std::string testFile = "/tmp/rodin_shard_meta_layout.h5";

    // -- Write phase: create a file with the shard layout -------------------
    {
      const auto file = HDF5::File(
          H5Fcreate(testFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      ASSERT_TRUE(file);

      // Create the group hierarchy for dimension 0
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::Shard,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::ShardFlags,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::ShardPolytopeMap,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::ShardOwner,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::ShardHalo,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }

      // /Shard/Flags/0 — ownership flags (U8)
      {
        std::vector<HDF5::U8> flags = { 1, 1, 2, 0 };
        HDF5::writeVectorDataset(file.get(), HDF5::shardFlagsPath(0), flags);
      }

      // /Shard/PolytopeMap/0
      {
        const auto g = HDF5::Group(
            H5Gcreate2(file.get(), HDF5::shardPolytopeMapGroupPath(0).c_str(),
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> left = { 10, 20, 30, 40 };
        HDF5::writeVectorDataset(file.get(), HDF5::shardPolytopeMapLeftPath(0), left);
      }
      {
        const auto g = HDF5::Group(
            H5Gcreate2(file.get(), HDF5::shardPolytopeMapRightGroupPath(0).c_str(),
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> keys = { 10, 20, 30, 40 };
        std::vector<HDF5::U64> vals = { 0, 1, 2, 3 };
        HDF5::writeVectorDataset(file.get(), HDF5::shardPolytopeMapRightGroupPath(0) + "/Keys", keys);
        HDF5::writeVectorDataset(file.get(), HDF5::shardPolytopeMapRightGroupPath(0) + "/Values", vals);
      }

      // /Shard/Owner/0
      {
        const auto g = HDF5::Group(
            H5Gcreate2(file.get(), HDF5::shardOwnerGroupPath(0).c_str(),
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> keys = { 2 };
        std::vector<HDF5::U64> vals = { 3 };
        HDF5::writeVectorDataset(file.get(), HDF5::shardOwnerGroupPath(0) + "/Keys", keys);
        HDF5::writeVectorDataset(file.get(), HDF5::shardOwnerGroupPath(0) + "/Values", vals);
      }

      // /Shard/Halo/0  (CSR: keys=[0,1], offsets=[0,1,3], indices=[1,2,3])
      {
        const auto g = HDF5::Group(
            H5Gcreate2(file.get(), HDF5::shardHaloGroupPath(0).c_str(),
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> keys    = { 0, 1 };
        std::vector<HDF5::U64> offsets  = { 0, 1, 3 };
        std::vector<HDF5::U64> indices  = { 1, 2, 3 };
        HDF5::writeVectorDataset(file.get(), HDF5::shardHaloGroupPath(0) + "/Keys", keys);
        HDF5::writeVectorDataset(file.get(), HDF5::shardHaloGroupPath(0) + "/Offsets", offsets);
        HDF5::writeVectorDataset(file.get(), HDF5::shardHaloGroupPath(0) + "/Indices", indices);
      }
    }

    // -- Read-back phase: verify all groups/datasets exist -------------------
    {
      hid_t h5 = H5Fopen(testFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(h5, 0);

      EXPECT_TRUE(HDF5::exists(h5, HDF5::Path::Shard));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::Path::ShardFlags));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardFlagsPath(0)));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardPolytopeMapLeftPath(0)));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardPolytopeMapRightGroupPath(0) + "/Keys"));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardPolytopeMapRightGroupPath(0) + "/Values"));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardOwnerGroupPath(0) + "/Keys"));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardOwnerGroupPath(0) + "/Values"));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardHaloGroupPath(0) + "/Keys"));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardHaloGroupPath(0) + "/Offsets"));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardHaloGroupPath(0) + "/Indices"));

      // Verify data round-trip for flags
      {
        auto flags = HDF5::readVectorDataset<HDF5::U8>(h5, HDF5::shardFlagsPath(0));
        ASSERT_EQ(flags.size(), 4u);
        EXPECT_EQ(flags[0], 1);
        EXPECT_EQ(flags[1], 1);
        EXPECT_EQ(flags[2], 2);
        EXPECT_EQ(flags[3], 0);
      }

      // Verify data round-trip for polytope map left
      {
        auto left = HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardPolytopeMapLeftPath(0));
        ASSERT_EQ(left.size(), 4u);
        EXPECT_EQ(left[0], 10u);
        EXPECT_EQ(left[3], 40u);
      }

      // Verify halo CSR structure
      {
        auto keys    = HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardHaloGroupPath(0) + "/Keys");
        auto offsets = HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardHaloGroupPath(0) + "/Offsets");
        auto indices = HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardHaloGroupPath(0) + "/Indices");
        ASSERT_EQ(keys.size(), 2u);
        ASSERT_EQ(offsets.size(), 3u);
        ASSERT_EQ(indices.size(), 3u);
        EXPECT_EQ(offsets[0], 0u);
        EXPECT_EQ(offsets[2], 3u);
      }

      H5Fclose(h5);
    }

    std::remove(testFile.c_str());
  }

  // --- Shard path helpers ---------------------------------------------------

  TEST(ShardMetadata, PathHelpers)
  {
    EXPECT_EQ(HDF5::shardFlagsPath(0), "/Shard/Flags/0");
    EXPECT_EQ(HDF5::shardFlagsPath(2), "/Shard/Flags/2");
    EXPECT_EQ(HDF5::shardPolytopeMapGroupPath(1), "/Shard/PolytopeMap/1");
    EXPECT_EQ(HDF5::shardPolytopeMapLeftPath(3), "/Shard/PolytopeMap/3/Left");
    EXPECT_EQ(HDF5::shardPolytopeMapRightGroupPath(0), "/Shard/PolytopeMap/0/Right");
    EXPECT_EQ(HDF5::shardOwnerGroupPath(2), "/Shard/Owner/2");
    EXPECT_EQ(HDF5::shardHaloGroupPath(1), "/Shard/Halo/1");
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

  // ===========================================================================
  // Distributed XDMF file-per-rank tests (simulated with local meshes)
  // ===========================================================================

  /**
   * Parameterized test: simulates a 2-rank MPI XDMF export using the
   * distributed XDMF constructor. Each "rank" writes its own rank-specific
   * HDF5 file, and only the root rank writes the master XDMF XML.
   */
  class MPIDistributedXDMF : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(MPIDistributedXDMF, FilePerRankVisualization)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_mpi_dist_xdmf_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "sim";

    auto localMesh = makeMesh(type);
    const size_t numRanks = 2;

    // Simulate both ranks
    // Non-root ranks must run first so their mesh files exist before root's
    // close()/flush() reads all rank files.
    for (size_t r = 1; r < numRanks; ++r)
    {
      XDMF xdmf(stem, r, numRanks, /*rootRank=*/0);
      xdmf.setMesh(localMesh);
      xdmf.write(0.0);
      xdmf.close();
    }
    {
      XDMF xdmf(stem, /*rank=*/0, numRanks, /*rootRank=*/0);
      xdmf.setMesh(localMesh);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Verify per-rank HDF5 files
    for (size_t r = 0; r < numRanks; ++r)
    {
      const auto meshH5 = testDir / ("sim.r" + std::to_string(r) + ".mesh.h5");
      ASSERT_TRUE(boost::filesystem::exists(meshH5)) << "Missing " << meshH5;

      hid_t h5 = H5Fopen(meshH5.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(h5, 0);

      // Visualization data
      EXPECT_GE(H5Lexists(h5, "/Mesh/Geometry/Vertices", H5P_DEFAULT), 1);
      EXPECT_GE(H5Lexists(h5, "/Mesh/XDMF/Topology", H5P_DEFAULT), 1);

      // No persistence data
      EXPECT_EQ(H5Lexists(h5, "/Mesh/Connectivity", H5P_DEFAULT), 0);

      H5Fclose(h5);
    }

    // Verify master XDMF XML
    const auto xdmfFile = testDir / "sim.xdmf";
    ASSERT_TRUE(boost::filesystem::exists(xdmfFile));

    std::ifstream ifs(xdmfFile.string());
    std::ostringstream buf;
    buf << ifs.rdbuf();
    const auto xml = buf.str();

    EXPECT_NE(xml.find("CollectionType=\"Spatial\""), std::string::npos);
    EXPECT_NE(xml.find("sim.r0.mesh.h5"), std::string::npos);
    EXPECT_NE(xml.find("sim.r1.mesh.h5"), std::string::npos);

    boost::filesystem::remove_all(testDir);
  }

  TEST_P(MPIDistributedXDMF, FilePerRankWithAttributes)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
        "/tmp/rodin_mpi_dist_xdmf_attr_" + polytopeLabel(type);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "field";

    auto localMesh = makeMesh(type);
    P1 fes(localMesh);
    GridFunction gf(fes);
    gf.setName("u");
    gf = [](const Geometry::Point& p) { return p.x(); };

    const size_t numRanks = 2;

    // Write non-root ranks first so root's close()/flush() can read their files.
    for (size_t r = 1; r < numRanks; ++r)
    {
      XDMF xdmf(stem, r, numRanks);
      xdmf.setMesh(localMesh);
      xdmf.add("u", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.close();
    }
    {
      XDMF xdmf(stem, /*rank=*/0, numRanks);
      xdmf.setMesh(localMesh);
      xdmf.add("u", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Verify rank-specific attribute files
    for (size_t r = 0; r < numRanks; ++r)
    {
      const auto attrH5 = testDir / ("field.u.r" + std::to_string(r) + ".000000.h5");
      EXPECT_TRUE(boost::filesystem::exists(attrH5)) << "Missing " << attrH5;
    }

    // XDMF references attribute files
    const auto xdmfFile = testDir / "field.xdmf";
    ASSERT_TRUE(boost::filesystem::exists(xdmfFile));
    std::ifstream ifs(xdmfFile.string());
    std::ostringstream buf;
    buf << ifs.rdbuf();
    const auto xml = buf.str();
    EXPECT_NE(xml.find("field.u.r0.000000.h5"), std::string::npos);
    EXPECT_NE(xml.find("field.u.r1.000000.h5"), std::string::npos);

    boost::filesystem::remove_all(testDir);
  }

  INSTANTIATE_TEST_SUITE_P(
      AllDimensions,
      MPIDistributedXDMF,
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
