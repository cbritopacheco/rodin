/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <array>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <mpi.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO.h>
#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/IO.h>
#include <Rodin/MPI/Variational/P1.h>

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
        return LocalMesh::UniformGrid(type, {11});
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return LocalMesh::UniformGrid(type, {4, 4});
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Wedge:
        return LocalMesh::UniformGrid(type, {3, 3, 3});
      default:
        ADD_FAILURE() << "Unsupported polytope type for makeMesh";
        return Geometry::Mesh<Context::Local>();
    }
  }

  std::string polytopeLabel(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Segment:
        return "Segment1D";
      case Polytope::Type::Triangle:
        return "Triangle2D";
      case Polytope::Type::Quadrilateral:
        return "Quad2D";
      case Polytope::Type::Tetrahedron:
        return "Tet3D";
      case Polytope::Type::Hexahedron:
        return "Hex3D";
      case Polytope::Type::Pyramid:
        return "Pyramid3D";
      case Polytope::Type::Wedge:
        return "Wedge3D";
      default:
        return "Unknown";
    }
  }

  int worldRank()
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }

  boost::filesystem::path rankTempPath(const std::string& stem)
  {
    return boost::filesystem::temp_directory_path() /
      (stem + "_r" + std::to_string(worldRank()));
  }

  /// @brief MPI HDF5 stream overloads reject operations requiring a file path.
  TEST(ShardMetadata, RejectsStreamPersistence)
  {
    boost::mpi::environment environment;
    boost::mpi::communicator communicator(MPI_COMM_SELF, boost::mpi::comm_attach);
    const Context::MPI context(environment, communicator);
    Geometry::Mesh<Context::MPI> mesh(context);
    MeshLoader<FileFormat::HDF5, Context::MPI> loader(mesh);
    MeshPrinter<FileFormat::HDF5, Context::MPI> printer(mesh);
    std::stringstream stream;
    EXPECT_THROW(loader.load(stream), Alert::Exception);
    EXPECT_THROW(printer.print(stream), Alert::Exception);
  }

  /// @brief Malformed shard state, maps, ownership and halo CSR data are rejected on load.
  TEST(ShardMetadata, RejectsMalformedMetadata)
  {
    boost::mpi::environment environment;
    boost::mpi::communicator communicator(MPI_COMM_SELF, boost::mpi::comm_attach);
    const Context::MPI context(environment, communicator);
    auto mesh =
      Geometry::Mesh<Context::MPI>::UniformGrid(context, Polytope::Type::Segment, {3});
    for (int corruption = 0; corruption < 5; ++corruption)
    {
      SCOPED_TRACE(corruption);
      const auto path = rankTempPath("rodin-malformed-shard.h5");
      mesh.save(path, FileFormat::HDF5);
      Geometry::Mesh<Context::MPI> restored(context);
      MeshLoader<FileFormat::HDF5, Context::MPI> validLoader(restored);
      validLoader.load(path);
      ASSERT_EQ(restored.getDimension(), mesh.getDimension());
      ASSERT_EQ(restored.getVertexCount(), mesh.getVertexCount());
      ASSERT_EQ(restored.getCellCount(), mesh.getCellCount());
      {
        const auto file = HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
        ASSERT_TRUE(file);
        const auto replace = [&]<class T>(
                               const std::string& dataset, const std::vector<T>& values) {
          ASSERT_GE(H5Ldelete(file.get(), dataset.c_str(), H5P_DEFAULT), 0);
          HDF5::writeVectorDataset(file.get(), dataset, values);
        };
        if (corruption == 0)
        {
          replace(HDF5::shardStatePath(0), std::vector<HDF5::U8>{});
        }
        else if (corruption == 1)
        {
          replace(HDF5::shardPolytopeMapLeftPath(0), std::vector<HDF5::U64>{});
        }
        else if (corruption == 2)
        {
          replace(HDF5::shardOwnerGroupPath(0) + "/Values",
            std::vector<HDF5::U64>(mesh.getVertexCount() + 1, 0));
        }
        else if (corruption == 3)
        {
          replace(HDF5::shardHaloGroupPath(0) + "/Offsets", std::vector<HDF5::U64>{});
        }
        else
        {
          replace(HDF5::shardHaloGroupPath(0) + "/Keys", std::vector<HDF5::U64>{0});
          replace(HDF5::shardHaloGroupPath(0) + "/Offsets", std::vector<HDF5::U64>{0, 1});
          replace(HDF5::shardHaloGroupPath(0) + "/Indices", std::vector<HDF5::U64>{});
        }
      }
      Geometry::Mesh<Context::MPI> loaded(context);
      MeshLoader<FileFormat::HDF5, Context::MPI> loader(loaded);
      EXPECT_THROW(loader.load(path), Alert::Exception);
      EXPECT_EQ(loaded.getDimension(), 0u);
      boost::filesystem::remove(path);
    }
  }

  /// @brief Convenience exporters reject MPI-backed fields before filesystem effects.
  TEST(ShardMetadata, ConvenienceAttributeExportRequiresExplicitShard)
  {
    boost::mpi::environment environment;
    boost::mpi::communicator communicator;
    const Context::MPI context(environment, communicator);
    auto mesh =
      Geometry::Mesh<Context::MPI>::UniformGrid(context, Polytope::Type::Segment, {17});
    const auto path = rankTempPath("rodin-mpi-attribute-" +
      boost::filesystem::unique_path("%%%%-%%%%-%%%%.h5").string());
    ASSERT_FALSE(boost::filesystem::exists(path));
    const auto verify = [&]<class Range>() {
      using FES = P1<Range, Geometry::Mesh<Context::MPI>>;
      auto space = [&] {
        if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
          return FES(mesh, 2, 3);
        else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
          return FES(mesh, 2);
        else
          return FES(mesh);
      }();
      GridFunction field(space);
      field.getData().setConstant(-17);
      const auto before = field.getData().eval();
      for (const bool node : {false, true})
      {
        SCOPED_TRACE(node ? "node" : "cell");
        try
        {
          if (node)
            HDF5::writeXDMFNodeAttribute(field, path);
          else
            HDF5::writeXDMFCellAttribute(field, path);
          FAIL() << "A convenience exporter accepted an MPI-backed field";
        }
        catch (const Alert::Exception& error)
        {
          const std::string message = error.what();
          EXPECT_NE(
            message.find("mesh is not a Mesh<Context::Local>"), std::string::npos);
          EXPECT_NE(message.find("shard"), std::string::npos);
        }
        catch (const std::exception& error)
        {
          ADD_FAILURE() << "Expected a native exporter diagnostic: " << error.what();
        }
        EXPECT_EQ(field.getData(), before);
        EXPECT_FALSE(boost::filesystem::exists(path));
        boost::filesystem::remove(path);
      }
    };
    verify.template operator()<Real>();
    verify.template operator()<Math::SpatialVector<Real>>();
    verify.template operator()<Math::SpatialMatrix<Real>>();
  }

  // Parameterized test fixture
  class MPIMeshHDF5 : public ::testing::TestWithParam<Polytope::Type>
  {};

  // --- MPI Mesh HDF5 round-trip via MeshPrinter/MeshLoader -------------------

  /// @brief Verifies mesh round trip via printer loader for MPI mesh HDF 5 by checking exact expected values, MPI behavior.
  TEST_P(MPIMeshHDF5, MeshRoundTripViaPrinterLoader)
  {
    const auto type = GetParam();
    const auto meshFile =
      rankTempPath("rodin_mpi_hdf5_rt_" + polytopeLabel(type) + ".h5");

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

    boost::filesystem::remove(meshFile);
  }

  // --- MPI Mesh persistence does not create XDMF datasets -------------------

  /// @brief Verifies mesh persistence no XDMF for MPI mesh HDF 5 by checking exact expected values, MPI behavior.
  TEST_P(MPIMeshHDF5, MeshPersistenceNoXDMF)
  {
    const auto type = GetParam();
    const auto meshFile =
      rankTempPath("rodin_mpi_hdf5_noxdmf_" + polytopeLabel(type) + ".h5");

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
    boost::filesystem::remove(meshFile);
  }

  // --- MPI Mesh XDMF visualization -------------------------------------------

  /// @brief Verifies XDMF visualization topology for MPI mesh HDF 5 by checking exact expected values, true predicates, MPI behavior.
  TEST_P(MPIMeshHDF5, XDMFVisualizationTopology)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
      rankTempPath("rodin_mpi_xdmf_topo_" + polytopeLabel(type));
    boost::filesystem::remove_all(testDir);
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

  /// @brief Verifies XDMF write and close for MPI mesh HDF 5 by checking exact expected values, true predicates, MPI behavior.
  TEST_P(MPIMeshHDF5, XDMFWriteAndClose)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
      rankTempPath("rodin_mpi_xdmf_wc_" + polytopeLabel(type));
    boost::filesystem::remove_all(testDir);
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
    const auto testFile = rankTempPath("rodin_shard_meta_layout.h5");

    // -- Write phase: create a file with the shard layout -------------------
    {
      const auto file =
        HDF5::File(H5Fcreate(testFile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      ASSERT_TRUE(file);

      // Create the group hierarchy for dimension 0
      {
        const auto g = HDF5::Group(H5Gcreate2(
          file.get(), HDF5::Path::Shard, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(
          file.get(), HDF5::Path::ShardState, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(), HDF5::Path::ShardPolytopeMap,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(
          file.get(), HDF5::Path::ShardOwner, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        const auto g = HDF5::Group(H5Gcreate2(
          file.get(), HDF5::Path::ShardHalo, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }

      // /Shard/Flags/0 — ownership flags (U8)
      {
        std::vector<HDF5::U8> flags = {1, 1, 2, 0};
        HDF5::writeVectorDataset(file.get(), HDF5::shardStatePath(0), flags);
      }

      // /Shard/PolytopeMap/0
      {
        const auto g =
          HDF5::Group(H5Gcreate2(file.get(), HDF5::shardPolytopeMapGroupPath(0).c_str(),
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> left = {10, 20, 30, 40};
        HDF5::writeVectorDataset(file.get(), HDF5::shardPolytopeMapLeftPath(0), left);
      }
      {
        const auto g = HDF5::Group(
          H5Gcreate2(file.get(), HDF5::shardPolytopeMapRightGroupPath(0).c_str(),
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> keys = {10, 20, 30, 40};
        std::vector<HDF5::U64> vals = {0, 1, 2, 3};
        HDF5::writeVectorDataset(
          file.get(), HDF5::shardPolytopeMapRightGroupPath(0) + "/Keys", keys);
        HDF5::writeVectorDataset(
          file.get(), HDF5::shardPolytopeMapRightGroupPath(0) + "/Values", vals);
      }

      // /Shard/Owner/0
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(),
          HDF5::shardOwnerGroupPath(0).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> keys = {2};
        std::vector<HDF5::U64> vals = {3};
        HDF5::writeVectorDataset(
          file.get(), HDF5::shardOwnerGroupPath(0) + "/Keys", keys);
        HDF5::writeVectorDataset(
          file.get(), HDF5::shardOwnerGroupPath(0) + "/Values", vals);
      }

      // /Shard/Halo/0  (CSR: keys=[0,1], offsets=[0,1,3], indices=[1,2,3])
      {
        const auto g = HDF5::Group(H5Gcreate2(file.get(),
          HDF5::shardHaloGroupPath(0).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        ASSERT_TRUE(g);
      }
      {
        std::vector<HDF5::U64> keys = {0, 1};
        std::vector<HDF5::U64> offsets = {0, 1, 3};
        std::vector<HDF5::U64> indices = {1, 2, 3};
        HDF5::writeVectorDataset(file.get(), HDF5::shardHaloGroupPath(0) + "/Keys", keys);
        HDF5::writeVectorDataset(
          file.get(), HDF5::shardHaloGroupPath(0) + "/Offsets", offsets);
        HDF5::writeVectorDataset(
          file.get(), HDF5::shardHaloGroupPath(0) + "/Indices", indices);
      }
    }

    // -- Read-back phase: verify all groups/datasets exist -------------------
    {
      hid_t h5 = H5Fopen(testFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      ASSERT_GE(h5, 0);

      EXPECT_TRUE(HDF5::exists(h5, HDF5::Path::Shard));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::Path::ShardState));
      EXPECT_TRUE(HDF5::exists(h5, HDF5::shardStatePath(0)));
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
        auto flags = HDF5::readVectorDataset<HDF5::U8>(h5, HDF5::shardStatePath(0));
        ASSERT_EQ(flags.size(), 4u);
        EXPECT_EQ(flags[0], 1);
        EXPECT_EQ(flags[1], 1);
        EXPECT_EQ(flags[2], 2);
        EXPECT_EQ(flags[3], 0);
      }

      // Verify data round-trip for polytope map left
      {
        auto left =
          HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardPolytopeMapLeftPath(0));
        ASSERT_EQ(left.size(), 4u);
        EXPECT_EQ(left[0], 10u);
        EXPECT_EQ(left[3], 40u);
      }

      // Verify halo CSR structure
      {
        auto keys =
          HDF5::readVectorDataset<HDF5::U64>(h5, HDF5::shardHaloGroupPath(0) + "/Keys");
        auto offsets = HDF5::readVectorDataset<HDF5::U64>(
          h5, HDF5::shardHaloGroupPath(0) + "/Offsets");
        auto indices = HDF5::readVectorDataset<HDF5::U64>(
          h5, HDF5::shardHaloGroupPath(0) + "/Indices");
        ASSERT_EQ(keys.size(), 2u);
        ASSERT_EQ(offsets.size(), 3u);
        ASSERT_EQ(indices.size(), 3u);
        EXPECT_EQ(offsets[0], 0u);
        EXPECT_EQ(offsets[2], 3u);
      }

      H5Fclose(h5);
    }

    boost::filesystem::remove(testFile);
  }

  // --- Shard path helpers ---------------------------------------------------

  /// @brief Verifies path helpers for shard metadata by checking exact expected values.
  TEST(ShardMetadata, PathHelpers)
  {
    EXPECT_EQ(HDF5::shardStatePath(0), "/Shard/State/0");
    EXPECT_EQ(HDF5::shardStatePath(2), "/Shard/State/2");
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

  INSTANTIATE_TEST_SUITE_P(AllDimensions, MPIMeshHDF5,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge),
    MPIPolytopeNameGenerator());

  // ===========================================================================
  // Distributed XDMF tests using MPI_COMM_SELF (single-rank distributed mode)
  // ===========================================================================

  /**
   * Tests the distributed XDMF code path using MPI_COMM_SELF
   * (1 rank, exercises the collective logic without multi-process setup).
   */
  class MPIDistributedXDMF : public ::testing::TestWithParam<Polytope::Type>
  {};

  /// @brief Verifies file per rank visualization for MPI distributed XDMF by checking exact expected values, true predicates, MPI behavior.
  TEST_P(MPIDistributedXDMF, FilePerRankVisualization)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
      rankTempPath("rodin_mpi_dist_xdmf_" + polytopeLabel(type));
    boost::filesystem::remove_all(testDir);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "sim";

    auto localMesh = makeMesh(type);

    {
      XDMF xdmf(MPI_COMM_SELF, stem);
      xdmf.setMesh(localMesh);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Verify rank 0 HDF5 file
    const auto meshH5 = testDir / "sim.r0.mesh.h5";
    ASSERT_TRUE(boost::filesystem::exists(meshH5)) << "Missing " << meshH5;

    {
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

    boost::filesystem::remove_all(testDir);
  }

  /// @brief Verifies file per rank with attributes for MPI distributed XDMF by checking true predicates, MPI behavior.
  TEST_P(MPIDistributedXDMF, FilePerRankWithAttributes)
  {
    const auto type = GetParam();
    const boost::filesystem::path testDir =
      rankTempPath("rodin_mpi_dist_xdmf_attr_" + polytopeLabel(type));
    boost::filesystem::remove_all(testDir);
    boost::filesystem::create_directories(testDir);
    const boost::filesystem::path stem = testDir / "field";

    auto localMesh = makeMesh(type);
    P1 fes(localMesh);
    GridFunction gf(fes);
    gf.setName("u");
    gf = [](const Geometry::Point& p) { return p.x(); };

    {
      XDMF xdmf(MPI_COMM_SELF, stem);
      xdmf.setMesh(localMesh);
      xdmf.add("u", gf, XDMF::Center::Node);
      xdmf.write(0.0);
      xdmf.close();
    }

    // Verify rank-specific attribute files
    const auto attrH5 = testDir / "field.u.r0.000000.h5";
    EXPECT_TRUE(boost::filesystem::exists(attrH5)) << "Missing " << attrH5;

    // XDMF references attribute files
    const auto xdmfFile = testDir / "field.xdmf";
    ASSERT_TRUE(boost::filesystem::exists(xdmfFile));
    std::ifstream ifs(xdmfFile.string());
    std::ostringstream buf;
    buf << ifs.rdbuf();
    const auto xml = buf.str();
    EXPECT_NE(xml.find("field.u.r0.000000.h5"), std::string::npos);

    boost::filesystem::remove_all(testDir);
  }

  /// @brief Every rank retains its shapes, static values and transient snapshots.
  TEST(MPIDistributedXDMFSeries, MixedShapesWithNonzeroPublishingRoot)
  {
    const int rank = worldRank();
    int ranks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    const int root = ranks - 1;
    std::array<char, 256> directory{};
    if (rank == root)
    {
      const auto path = (boost::filesystem::temp_directory_path() /
        boost::filesystem::unique_path("rodin-mpi-series-%%%%-%%%%"))
                          .string();
      ASSERT_LT(path.size(), directory.size());
      std::copy(path.begin(), path.end(), directory.begin());
      boost::filesystem::create_directory(path);
    }
    MPI_Bcast(directory.data(), directory.size(), MPI_CHAR, root, MPI_COMM_WORLD);
    const boost::filesystem::path testDir(directory.data());
    auto mesh = LocalMesh::UniformGrid(
      Polytope::Type::Triangle, {static_cast<size_t>(rank + 2), 2});
    P1 tensorSpace(mesh, 3, 3);
    P0 matrixSpace(mesh, 2, 3);
    P1 vectorSpace(mesh, 2);
    GridFunction tensor(tensorSpace);
    GridFunction matrix(matrixSpace);
    GridFunction vector(vectorSpace);
    Math::SpatialMatrix<Real> tensorValue(3, 3);
    Math::SpatialMatrix<Real> matrixValue(2, 3);
    for (size_t i = 0; i < 3; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
        tensorValue(i, j) = 100 * rank + 10 * i + j + 1;
    }
    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
        matrixValue(i, j) = 100 * rank + 10 * i + j + 2;
    }
    tensor = tensorValue;
    matrix = matrixValue;
    vector = Math::SpatialVector<Real>{Real(rank + 3), Real(-rank - 4)};
    XDMF output(MPI_COMM_WORLD, testDir / "series", root);
    XDMF::FilePatterns patterns;
    patterns.transientMesh = "mesh-r{rank}-{index}.h5";
    patterns.staticAttribute = "{name}-r{rank}-fixed.h5";
    patterns.attribute = "{name}-r{rank}-{index}.h5";
    output.setFilePatterns(patterns)
      .setPadding(1)
      .setMesh(mesh, XDMF::MeshPolicy::Transient)
      .add("tensor", tensor, XDMF::Center::Node, XDMF::AttributePolicy::Static)
      .add("matrix", matrix, XDMF::Center::Cell)
      .add("vector", vector, XDMF::Center::Node);
    output.write(Real(0.25));
    tensor = Real(1000) * tensorValue;
    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
        matrixValue(i, j) += Real(1000);
    }
    matrix = matrixValue;
    vector = Math::SpatialVector<Real>{Real(rank + 13), Real(-rank + 6)};
    output.write(Real(0.75)).close();
    EXPECT_EQ(output.getSnapshotCount(), 2u);
    EXPECT_EQ(output.getRootRank(), static_cast<size_t>(root));
    EXPECT_EQ(output.getNumRanks(), static_cast<size_t>(ranks));
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != root)
      return;
    std::ifstream input((testDir / "series.xdmf").string());
    ASSERT_TRUE(input.good());
    const std::string xml((std::istreambuf_iterator<char>(input)), {});
    EXPECT_NE(xml.find("Time Value=\"0.25\""), std::string::npos);
    EXPECT_NE(xml.find("Time Value=\"0.75\""), std::string::npos);
    for (int piece = 0; piece < ranks; ++piece)
    {
      SCOPED_TRACE(piece);
      const size_t vertices = 2 * (piece + 2);
      const size_t cells = 2 * (piece + 1);
      for (size_t step = 0; step < 2; ++step)
      {
        SCOPED_TRACE(step);
        const auto meshName =
          "mesh-r" + std::to_string(piece) + "-" + std::to_string(step) + ".h5";
        EXPECT_TRUE(boost::filesystem::exists(testDir / meshName));
        EXPECT_NE(xml.find(meshName + ":/Mesh/Geometry/Vertices"), std::string::npos);
        for (const std::string name : {"tensor", "matrix", "vector"})
        {
          SCOPED_TRACE(name);
          const bool isTensor = name == "tensor";
          const bool isMatrix = name == "matrix";
          const size_t count = isMatrix ? cells : vertices;
          const size_t rows = isTensor ? 3 : 2;
          const size_t columns = name == "vector" ? 1 : 3;
          const auto fileName = name + "-r" + std::to_string(piece) +
            (isTensor ? "-fixed.h5" : "-" + std::to_string(step) + ".h5");
          const std::string dimensions = std::to_string(count) + " " +
            std::to_string(rows) + (columns == 1 ? "" : " " + std::to_string(columns));
          const std::string dataItem = "Dimensions=\"" + dimensions + "\">" + fileName +
            ":/GridFunction/Values/Data</DataItem>";
          EXPECT_NE(xml.find(dataItem), std::string::npos);
          const auto file = HDF5::File(
            H5Fopen((testDir / fileName).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
          ASSERT_TRUE(file);
          const auto dataset =
            HDF5::DataSet(H5Dopen2(file.get(), "/GridFunction/Values/Data", H5P_DEFAULT));
          ASSERT_TRUE(dataset);
          const auto space = HDF5::Space(H5Dget_space(dataset.get()));
          ASSERT_TRUE(space);
          std::array<hsize_t, 3> shape{};
          ASSERT_EQ(H5Sget_simple_extent_dims(space.get(), shape.data(), nullptr),
            columns == 1 ? 2 : 3);
          EXPECT_EQ(shape[0], count);
          EXPECT_EQ(shape[1], rows);
          if (columns != 1)
          {
            EXPECT_EQ(shape[2], columns);
          }
          std::vector<double> values(count * rows * columns);
          ASSERT_GE(H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, values.data()),
            0);
          for (size_t entity = 0; entity < count; ++entity)
          {
            for (size_t i = 0; i < rows; ++i)
            {
              for (size_t j = 0; j < columns; ++j)
              {
                const Real expected = columns == 1
                  ? (i == 0 ? piece + 3 : -piece - 4) + Real(10 * step)
                  : 100 * piece + 10 * i + j + (isTensor ? 1 : 2 + 1000 * step);
                EXPECT_DOUBLE_EQ(values[(entity * rows + i) * columns + j], expected);
              }
            }
          }
        }
      }
    }
    EXPECT_NE(xml.find("AttributeType=\"Tensor\" Center=\"Node\""), std::string::npos);
    EXPECT_NE(xml.find("AttributeType=\"Matrix\" Center=\"Cell\""), std::string::npos);
    EXPECT_NE(xml.find("AttributeType=\"Vector\" Center=\"Node\""), std::string::npos);
    boost::filesystem::remove_all(testDir);
  }

  INSTANTIATE_TEST_SUITE_P(AllDimensions, MPIDistributedXDMF,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron, Polytope::Type::Pyramid, Polytope::Type::Wedge),
    MPIPolytopeNameGenerator());
}

int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

/// @brief Default distributed static filenames keep each rank's values in a separate file.
TEST(MPIDistributedXDMFSeries, DefaultStaticAttributesAreRankLocal)
{
  const int rank = worldRank();
  int ranks = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  std::array<char, 256> directory{};
  if (rank == 0)
  {
    const auto path = (boost::filesystem::temp_directory_path() /
      boost::filesystem::unique_path("rodin-mpi-static-%%%%-%%%%"))
                        .string();
    std::copy(path.begin(), path.end(), directory.begin());
    boost::filesystem::create_directory(path);
  }
  MPI_Bcast(directory.data(), directory.size(), MPI_CHAR, 0, MPI_COMM_WORLD);
  const boost::filesystem::path testDir(directory.data());
  const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  P1 space(mesh);
  GridFunction field(space);
  field = Real(rank + 1);
  XDMF output(MPI_COMM_WORLD, testDir / "series");
  output.setMesh(mesh).add(
    "value", field, XDMF::Center::Node, XDMF::AttributePolicy::Static);
  output.write().close();
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
  {
    for (int r = 0; r < ranks; ++r)
    {
      const auto path = testDir / ("series.value.r" + std::to_string(r) + ".h5");
      ASSERT_TRUE(boost::filesystem::exists(path));
      auto file = HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      const auto values = HDF5::readVectorDataset<HDF5::F64>(
        file.get(), HDF5::Path::GridFunctionValuesData);
      ASSERT_FALSE(values.empty());
      for (const auto value : values)
        EXPECT_EQ(value, Real(r + 1));
    }
    boost::filesystem::remove_all(testDir);
  }
}
