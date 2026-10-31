/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <hdf5.h>

#include <Rodin/IO/XDMF.h>

#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /**
   * @brief Reads a whole dataset as doubles.
   * @param path HDF5 file.
   * @param dataset Dataset path inside the file.
   * @return Values of the dataset in storage order.
   */
  std::vector<double> readDataset(
    const boost::filesystem::path& path, const char* dataset)
  {
    const auto file = IO::HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    EXPECT_TRUE(file);
    const auto data = IO::HDF5::DataSet(H5Dopen2(file.get(), dataset, H5P_DEFAULT));
    EXPECT_TRUE(data);
    const auto space = IO::HDF5::Space(H5Dget_space(data.get()));
    std::vector<double> values(H5Sget_simple_extent_npoints(space.get()));
    EXPECT_GE(H5Dread(data.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                values.data()),
      0);
    return values;
  }

  /// @brief Isolates the files written by each test.
  class Kokkos_IO_XDMF : public ::testing::Test
  {
    protected:
      /// @brief Creates an independent output directory.
      void SetUp() override
      {
        m_directory = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-kokkos-xdmf-%%%%-%%%%");
        boost::filesystem::create_directory(m_directory);
      }

      /// @brief Removes the output directory.
      void TearDown() override
      {
        boost::filesystem::remove_all(m_directory);
      }

      /**
       * @brief Writes a mesh and a field with fixed file names.
       * @param stem Output stem inside the directory.
       * @param mesh Mesh of the grid.
       * @param gf Field written on the vertices.
       * @param policy Mesh export policy.
       */
      template <class MeshType, class GF>
      void write(const std::string& stem, const MeshType& mesh, const GF& gf,
        IO::XDMF::MeshPolicy policy = IO::XDMF::MeshPolicy::Static)
      {
        IO::XDMF writer(m_directory / stem);
        IO::XDMF::FilePatterns patterns;
        patterns.staticMesh = stem + ".mesh.h5";
        patterns.transientMesh = stem + ".mesh.{index}.h5";
        patterns.attribute = stem + ".{index}.h5";
        writer.setFilePatterns(patterns).setPadding(1);
        writer.grid().setMesh(mesh, policy).add("u", gf);
        writer.write(0.0).write(1.0).close();
      }

      /**
       * @brief Expects two written datasets to be equal.
       * @param lhs First HDF5 file.
       * @param rhs Second HDF5 file.
       * @param dataset Dataset path inside both files.
       */
      void expectSameDataset(
        const std::string& lhs, const std::string& rhs, const char* dataset) const
      {
        const auto expected = readDataset(m_directory / lhs, dataset);
        const auto actual = readDataset(m_directory / rhs, dataset);
        ASSERT_EQ(actual.size(), expected.size()) << dataset;
        for (size_t i = 0; i < expected.size(); ++i)
          EXPECT_DOUBLE_EQ(actual[i], expected[i]) << dataset << " entry " << i;
      }

      /**
       * @brief Device output matches the host output of the same field.
       * @param source Host mesh.
       * @param vdim Number of field components.
       */
      void checkSameOutput(const HostMesh& source, size_t vdim)
      {
        const DeviceMesh mesh(TestContext{}, source);
        P1 hfes(source, vdim);
        P1 dfes(mesh, vdim);
        GridFunction hg(hfes);
        fill(hg, source, vdim);
        GridFunction dg(dfes);
        dg.getData() = hg.getData();

        write("host", source, hg);
        write("device", mesh, dg);
        for (const char* dataset : {"/Mesh/Geometry/Vertices", "/Mesh/XDMF/Topology"})
          expectSameDataset("host.mesh.h5", "device.mesh.h5", dataset);
        for (const char* attribute : {"0", "1"})
        {
          expectSameDataset(std::string("host.") + attribute + ".h5",
            std::string("device.") + attribute + ".h5", "/GridFunction/Values/Data");
        }
      }

      boost::filesystem::path m_directory; ///< Output directory.
  };

  TEST_F(Kokkos_IO_XDMF, ScalarFieldOnTriangles)
  {
    checkSameOutput(makeTriangleMesh(), 1);
  }

  TEST_F(Kokkos_IO_XDMF, VectorFieldOnMixedMesh)
  {
    checkSameOutput(makeMixedMesh(), 2);
  }

  TEST_F(Kokkos_IO_XDMF, ScalarFieldOnTetrahedra)
  {
    checkSameOutput(makeTetrahedronMesh(), 1);
  }

  /// @brief A transient device mesh is copied again at every snapshot.
  TEST_F(Kokkos_IO_XDMF, TransientMeshFollowsDeviceCoordinates)
  {
    const HostMesh source = makeTriangleMesh();
    DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    GridFunction hg(hfes);
    fill(hg, source, 1);
    GridFunction dg(dfes);
    dg.getData() = hg.getData();

    IO::XDMF writer(m_directory / "moving");
    IO::XDMF::FilePatterns patterns;
    patterns.transientMesh = "moving.mesh.{index}.h5";
    patterns.attribute = "moving.{index}.h5";
    writer.setFilePatterns(patterns).setPadding(1);
    writer.grid().setMesh(mesh, IO::XDMF::MeshPolicy::Transient).add("u", dg);
    writer.write(0.0);
    mesh.scale(2);
    writer.write(1.0).close();

    const auto before =
      readDataset(m_directory / "moving.mesh.0.h5", "/Mesh/Geometry/Vertices");
    const auto after =
      readDataset(m_directory / "moving.mesh.1.h5", "/Mesh/Geometry/Vertices");
    ASSERT_EQ(after.size(), before.size());
    for (size_t i = 0; i < before.size(); ++i)
      EXPECT_DOUBLE_EQ(after[i], 2 * before[i]);
  }
}
