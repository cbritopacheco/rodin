/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <string>
#include <vector>

#include "Rodin/Geometry.h"
#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/IO/HDF5.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::IO;

namespace
{
  Optional<std::string> rejectedGroup;
  Optional<size_t> rejectedOpen;
  size_t openCalls = 0;
  size_t failures = 0;
}

extern "C"
{
  hid_t __real_H5Gcreate2(
    hid_t file, const char* name, hid_t link, hid_t creation, hid_t access);
  hid_t __wrap_H5Gcreate2(
    hid_t file, const char* name, hid_t link, hid_t creation, hid_t access)
  {
    if (rejectedGroup && *rejectedGroup == name)
    {
      ++failures;
      return -1;
    }
    return __real_H5Gcreate2(file, name, link, creation, access);
  }

  hid_t __real_H5Fopen(const char* name, unsigned flags, hid_t access);
  hid_t __wrap_H5Fopen(const char* name, unsigned flags, hid_t access)
  {
    if (rejectedOpen && ++openCalls == *rejectedOpen)
    {
      ++failures;
      return -1;
    }
    return __real_H5Fopen(name, flags, access);
  }
}

namespace Rodin::Tests::Unit
{
  class MPIHDF5FailureTest : public ::testing::Test
  {
    protected:
      boost::filesystem::path filename;

      void SetUp() override
      {
        rejectedGroup.reset();
        rejectedOpen.reset();
        failures = 0;
        filename = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-mpi-hdf5-failure-%%%%-%%%%.h5");
      }

      void TearDown() override
      {
        rejectedGroup.reset();
        rejectedOpen.reset();
        EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
        boost::filesystem::remove(filename);
      }

      template <class Action>
      void expectFailure(Action&& action, const std::string& diagnostic)
      {
        failures = 0;
        try
        {
          action();
          ADD_FAILURE() << "The selected backend failure was accepted";
        }
        catch (const Alert::Exception& error)
        {
          EXPECT_NE(std::string(error.what()).find(diagnostic), std::string::npos);
        }
        rejectedGroup.reset();
        rejectedOpen.reset();
        EXPECT_EQ(failures, 1u);
        EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
      }

      void expectShard(const Geometry::Mesh<Context::MPI>& restored,
        const Geometry::Mesh<Context::MPI>& expected)
      {
        EXPECT_EQ(restored.getDimension(), expected.getDimension());
        EXPECT_EQ(restored.getShard().getVertices().getMatrix(),
          expected.getShard().getVertices().getMatrix());
        const auto& actual = restored.getShard();
        const auto& original = expected.getShard();
        for (size_t d = 0; d <= expected.getDimension(); ++d)
        {
          SCOPED_TRACE(d);
          EXPECT_EQ(actual.getState(d), original.getState(d));
          EXPECT_EQ(actual.getPolytopeMap(d).left, original.getPolytopeMap(d).left);
          EXPECT_EQ(actual.getPolytopeMap(d).right, original.getPolytopeMap(d).right);
          EXPECT_EQ(actual.getOwner(d), original.getOwner(d));
          EXPECT_EQ(actual.getHalo(d), original.getHalo(d));
        }
      }
  };

  /// @brief Every shard metadata group failure releases files before collective recovery.
  TEST_F(MPIHDF5FailureTest, ShardGroupFailuresPermitCollectiveRetry)
  {
    boost::mpi::environment environment;
    boost::mpi::communicator communicator;
    const Context::MPI context(environment, communicator);
    const auto mesh = Geometry::Mesh<Context::MPI>::UniformGrid(
      context, Polytope::Type::Triangle, {4, 4});
    MeshPrinter<FileFormat::HDF5, Context::MPI> printer(mesh);
    const std::vector<std::string> groups{HDF5::Path::Shard, HDF5::Path::ShardState,
      HDF5::Path::ShardPolytopeMap, HDF5::Path::ShardOwner, HDF5::Path::ShardHalo,
      HDF5::shardPolytopeMapGroupPath(0), HDF5::shardPolytopeMapRightGroupPath(0),
      HDF5::shardOwnerGroupPath(0), HDF5::shardHaloGroupPath(0)};
    for (const auto& group : groups)
    {
      SCOPED_TRACE(group);
      rejectedGroup = group;
      expectFailure([&] { printer.print(filename); }, "Failed to create " + group);
      printer.print(filename);
      Geometry::Mesh<Context::MPI> restored(context);
      MeshLoader<FileFormat::HDF5, Context::MPI>(restored).load(filename);
      expectShard(restored, mesh);
    }
  }

  /// @brief Failed metadata reopens preserve mesh state and allow a later collective load.
  TEST_F(MPIHDF5FailureTest, MetadataReopenFailuresPermitCollectiveRetry)
  {
    boost::mpi::environment environment;
    boost::mpi::communicator communicator;
    const Context::MPI context(environment, communicator);
    const auto mesh = Geometry::Mesh<Context::MPI>::UniformGrid(
      context, Polytope::Type::Triangle, {4, 4});
    MeshPrinter<FileFormat::HDF5, Context::MPI> printer(mesh);
    openCalls = 0;
    rejectedOpen = 1;
    expectFailure([&] { printer.print(filename); }, "Failed to reopen HDF5 file");
    printer.print(filename);
    Geometry::Mesh<Context::MPI> restored(context);
    MeshLoader<FileFormat::HDF5, Context::MPI> loader(restored);
    openCalls = 0;
    rejectedOpen = 2;
    expectFailure(
      [&] { loader.load(filename); }, "Failed to open HDF5 file for shard metadata");
    EXPECT_EQ(restored.getCellCount(), 0u);
    loader.load(filename);
    expectShard(restored, mesh);
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
