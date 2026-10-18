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
#include "Rodin/Variational.h"
#include "Rodin/IO/HDF5.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace
{
  enum class Operation
  {
    CreateSpace,
    GetSpace,
    Rank,
    Dimensions,
    Write,
    CreateGroup
  };

  struct Failure
  {
    Failure()
      : hits(0)
    {}
    Optional<Operation> operation;
    std::string group;
    size_t hits;

    void arm(Operation selected, std::string path = {})
    {
      operation = selected;
      group = std::move(path);
      hits = 0;
    }

    bool reject(Operation selected, const char* path = nullptr)
    {
      if (operation != selected ||
        (selected == Operation::CreateGroup && group != path))
        return false;
      ++hits;
      return true;
    }
  };

  Failure failure;
}

extern "C"
{
  hid_t __real_H5Screate_simple(int rank, const hsize_t* dims, const hsize_t* maxdims);
  hid_t __wrap_H5Screate_simple(int rank, const hsize_t* dims, const hsize_t* maxdims)
  {
    return failure.reject(Operation::CreateSpace)
      ? -1
      : __real_H5Screate_simple(rank, dims, maxdims);
  }

  hid_t __real_H5Dget_space(hid_t dataset);
  hid_t __wrap_H5Dget_space(hid_t dataset)
  {
    return failure.reject(Operation::GetSpace) ? -1 : __real_H5Dget_space(dataset);
  }

  int __real_H5Sget_simple_extent_ndims(hid_t space);
  int __wrap_H5Sget_simple_extent_ndims(hid_t space)
  {
    return failure.reject(Operation::Rank) ? -1
                                           : __real_H5Sget_simple_extent_ndims(space);
  }

  int __real_H5Sget_simple_extent_dims(hid_t space, hsize_t* dims, hsize_t* maxdims);
  int __wrap_H5Sget_simple_extent_dims(hid_t space, hsize_t* dims, hsize_t* maxdims)
  {
    return failure.reject(Operation::Dimensions)
      ? -1
      : __real_H5Sget_simple_extent_dims(space, dims, maxdims);
  }

  herr_t __real_H5Dwrite(hid_t dataset, hid_t type, hid_t memory, hid_t file,
    hid_t properties, const void* values);
  herr_t __wrap_H5Dwrite(hid_t dataset, hid_t type, hid_t memory, hid_t file,
    hid_t properties, const void* values)
  {
    return failure.reject(Operation::Write)
      ? -1
      : __real_H5Dwrite(dataset, type, memory, file, properties, values);
  }

  hid_t __real_H5Gcreate2(
    hid_t file, const char* name, hid_t link, hid_t creation, hid_t access);
  hid_t __wrap_H5Gcreate2(
    hid_t file, const char* name, hid_t link, hid_t creation, hid_t access)
  {
    return failure.reject(Operation::CreateGroup, name)
      ? -1
      : __real_H5Gcreate2(file, name, link, creation, access);
  }
}

namespace Rodin::Tests::Unit
{
  class HDF5FailureTest : public ::testing::Test
  {
    protected:
      boost::filesystem::path filename;

      void SetUp() override
      {
        failure.operation.reset();
        filename = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-hdf5-failure-%%%%-%%%%.h5");
      }

      void TearDown() override
      {
        failure.operation.reset();
        EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
        boost::filesystem::remove(filename);
      }

      template <class Action>
      void expectFailure(Action&& action, const std::string& diagnostic)
      {
        try
        {
          action();
          FAIL() << "The backend failure was accepted";
        }
        catch (const Alert::Exception& error)
        {
          EXPECT_NE(std::string(error.what()).find(diagnostic), std::string::npos);
        }
        EXPECT_EQ(failure.hits, 1u);
        failure.operation.reset();
      }
  };

  /// @brief Failed dataspace creation and writes release datasets and allow replacement.
  TEST_F(HDF5FailureTest, DatasetOutputFailuresPermitRetry)
  {
    auto file =
      HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_TRUE(file);
    const std::vector<HDF5::F64> values{1, -2, 3, -4};
    for (const auto operation : {Operation::CreateSpace, Operation::Write})
    {
      for (bool matrix : {false, true})
      {
        SCOPED_TRACE(matrix);
        const auto write = [&] {
          if (matrix)
            HDF5::writeMatrixDataset(file.get(), "/Data", values, 2, 2);
          else
            HDF5::writeVectorDataset(file.get(), "/Data", values);
        };
        failure.arm(operation);
        expectFailure(write,
          operation == Operation::CreateSpace ? "Failed to create HDF5 dataspace"
                                              : "Failed to write HDF5 dataset");
        EXPECT_EQ(H5Fget_obj_count(file.get(), H5F_OBJ_ALL), 1);
        if (H5Lexists(file.get(), "/Data", H5P_DEFAULT) > 0)
        {
          ASSERT_GE(H5Ldelete(file.get(), "/Data", H5P_DEFAULT), 0);
        }
        write();
        EXPECT_EQ(HDF5::readVectorDataset<HDF5::F64>(file.get(), "/Data"), values);
        ASSERT_GE(H5Ldelete(file.get(), "/Data", H5P_DEFAULT), 0);
      }
    }
  }

  /// @brief Failed reader queries close their dataset and leave the data readable.
  TEST_F(HDF5FailureTest, DatasetQueriesPermitRetry)
  {
    auto file =
      HDF5::File(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
    ASSERT_TRUE(file);
    const std::vector<HDF5::F64> values{1, -2, 3, -4};
    HDF5::writeMatrixDataset(file.get(), "/Data", values, 2, 2);
    for (size_t query = 0; query < 3; ++query)
    {
      for (const auto operation :
        {Operation::GetSpace, Operation::Rank, Operation::Dimensions})
      {
        if (query == 0 && operation != Operation::GetSpace)
          continue;
        SCOPED_TRACE(query);
        failure.arm(operation);
        expectFailure(
          [&] {
            if (query == 0)
              (void)HDF5::readVectorDataset<HDF5::F64>(file.get(), "/Data");
            else if (query == 1)
              (void)HDF5::readMatrixShape(file.get(), "/Data");
            else
              (void)HDF5::readDatasetShape(file.get(), "/Data");
          },
          query == 1 && operation == Operation::Rank ? "Expected rank-2" : "Failed to");
        EXPECT_EQ(H5Fget_obj_count(file.get(), H5F_OBJ_ALL), 1);
        EXPECT_EQ(HDF5::readVectorDataset<HDF5::F64>(file.get(), "/Data"), values);
        EXPECT_EQ(
          HDF5::readDatasetShape(file.get(), "/Data"), (std::vector<hsize_t>{2, 2}));
      }
    }
  }

  /// @brief Mesh group creation errors carry member context and leave no open handles.
  TEST_F(HDF5FailureTest, MeshGroupsPermitRetry)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    const std::vector<std::string> groups{HDF5::Path::Mesh, HDF5::Path::MeshMeta,
      HDF5::Path::MeshGeometry, HDF5::Path::MeshConnectivity,
      HDF5::Path::MeshConnectivityMeta, HDF5::Path::MeshConnectivityCounts,
      HDF5::Path::MeshConnectivityEntities, HDF5::Path::MeshConnectivityState,
      HDF5::Path::MeshConnectivityIncidence, HDF5::Path::MeshAttributes,
      HDF5::Path::MeshTransformations, HDF5::entityGroupPath(1),
      HDF5::incidenceGroupPath(2, 1), HDF5::transformationGroupPath(2)};
    MeshPrinter<FileFormat::HDF5, Context::Local> printer(mesh);
    for (const auto& group : groups)
    {
      SCOPED_TRACE(group);
      failure.arm(Operation::CreateGroup, group);
      expectFailure([&] { printer.print(filename); }, "In member function");
      EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
      printer.print(filename);
      LocalMesh restored;
      MeshLoader<FileFormat::HDF5, Context::Local>(restored).load(filename);
      EXPECT_EQ(restored.getCellCount(), mesh.getCellCount());
      EXPECT_EQ(restored.getVertices().getMatrix(), mesh.getVertices().getMatrix());
    }
  }

  /// @brief Field group creation errors release their file before a successful retry.
  TEST_F(HDF5FailureTest, FieldGroupsPermitRetry)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P0 fes(mesh);
    GridFunction field(fes);
    field = RealFunction(3);
    for (const auto* group : {HDF5::Path::GridFunction, HDF5::Path::GridFunctionMeta,
           HDF5::Path::GridFunctionValues})
    {
      SCOPED_TRACE(group);
      failure.arm(Operation::CreateGroup, group);
      expectFailure(
        [&] { field.save(filename, FileFormat::HDF5); }, "In member function");
      EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
      field.save(filename, FileFormat::HDF5);
      GridFunction restored(fes);
      restored.load(filename, FileFormat::HDF5);
      EXPECT_EQ(restored.getData(), field.getData());
    }
  }

  /// @brief Visualization group failures release files before a successful replacement.
  TEST_F(HDF5FailureTest, VisualizationGroupsPermitRetry)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    P1 nodeSpace(mesh);
    P0 cellSpace(mesh);
    GridFunction nodes(nodeSpace);
    GridFunction cells(cellSpace);
    nodes = RealFunction(3);
    cells = RealFunction(4);
    for (bool cellCentered : {false, true})
    {
      for (const auto* group : {HDF5::Path::GridFunction, HDF5::Path::GridFunctionMeta,
             HDF5::Path::GridFunctionValues})
      {
        SCOPED_TRACE(group);
        SCOPED_TRACE(cellCentered);
        const auto write = [&] {
          if (cellCentered)
            HDF5::writeXDMFCellAttribute(cells, mesh, filename);
          else
            HDF5::writeXDMFNodeAttribute(nodes, mesh, filename);
        };
        failure.arm(Operation::CreateGroup, group);
        expectFailure(write, std::string("Failed to create ") + group);
        EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
        write();
        auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        ASSERT_TRUE(file);
        EXPECT_EQ(HDF5::readVectorDataset<HDF5::F64>(
                    file.get(), HDF5::Path::GridFunctionValuesData),
          std::vector<HDF5::F64>(
            cellCentered ? mesh.getCellCount() : mesh.getVertexCount(),
            cellCentered ? 4 : 3));
      }
    }
    failure.arm(Operation::CreateGroup, HDF5::Path::Mesh);
    expectFailure(
      [&] { HDF5::writeXDMFMesh(filename, mesh); }, "Failed to create /Mesh group");
    EXPECT_EQ(H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL), 0);
    HDF5::writeXDMFMesh(filename, mesh);
    auto file = HDF5::File(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    ASSERT_TRUE(file);
    EXPECT_EQ(HDF5::readMatrixShape(file.get(), HDF5::Path::MeshGeometryVertices),
      (std::pair<hsize_t, hsize_t>{mesh.getVertexCount(), mesh.getSpaceDimension()}));
  }

}
