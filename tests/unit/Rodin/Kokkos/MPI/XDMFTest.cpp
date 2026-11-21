/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <hdf5.h>

#include <boost/mpi/collectives/broadcast.hpp>

#include <Rodin/IO/XDMF.h>

#include "Common.h"

namespace Rodin::Tests::Unit
{
  /**
   * @brief Reads a whole dataset as doubles.
   * @param path HDF5 file.
   * @param dataset Dataset path inside the file.
   * @return Values of the dataset in storage order.
   */
  std::vector<double> readDistributedDataset(
    const boost::filesystem::path& path, const char* dataset)
  {
    const auto file = IO::HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
    EXPECT_TRUE(file) << path;
    const auto data = IO::HDF5::DataSet(H5Dopen2(file.get(), dataset, H5P_DEFAULT));
    EXPECT_TRUE(data) << path << dataset;
    const auto space = IO::HDF5::Space(H5Dget_space(data.get()));
    std::vector<double> values(H5Sget_simple_extent_npoints(space.get()));
    EXPECT_GE(H5Dread(data.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                values.data()),
      0);
    return values;
  }

  /**
   * @brief Writes one field of a distributed mesh, one file set per rank.
   * @param directory Output directory shared by the ranks.
   * @param stem Output stem inside the directory.
   * @param mesh Distributed mesh of the grid.
   * @param gf Field written on the vertices.
   */
  template <class MeshType, class GF>
  void writeDistributedField(const boost::filesystem::path& directory,
    const std::string& stem, const MeshType& mesh, const GF& gf)
  {
    IO::XDMF writer(static_cast<MPI_Comm>(*g_world), directory / stem);
    IO::XDMF::FilePatterns patterns;
    patterns.staticMesh = stem + ".mesh.{rank}.h5";
    patterns.transientMesh = stem + ".mesh.{rank}.{index}.h5";
    patterns.attribute = stem + ".{rank}.{index}.h5";
    writer.setFilePatterns(patterns).setPadding(1);
    writer.grid().setMesh(mesh).add("u", gf);
    writer.write(0.0).close();
  }

  /**
   * @brief The XDMF output of a distributed device field matches the host
   * output of the same problem.
   */
  void checkDistributedXDMF()
  {
    using namespace Rodin::Variational;
    Context::MPI ctx(*g_env, *g_world);
    const HostMesh host = distributeFromRoot(ctx, Polytope::Type::Triangle, 12);
    const DeviceMesh mesh(DeviceContext(ctx), host);
    P1<Real, HostMesh> hfes(host);
    P1<Real, DeviceMesh> dfes(mesh);
    PETSc::Variational::TrialFunction hu(hfes);
    PETSc::Variational::TestFunction hv(hfes);
    PETSc::Variational::TrialFunction du(dfes);
    PETSc::Variational::TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction df(
      KOKKOS_LAMBDA(const ::Rodin::Kokkos::View<::Rodin::Geometry::Point>& p)->Real {
        return 1 + p.x() * p.y();
      });
    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) + Integral(hu, hv) - Integral(hf, hv);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv);
    solve(hp);
    solve(dp);

    long id = static_cast<long>(getpid());
    boost::mpi::broadcast(*g_world, id, 0);
    const boost::filesystem::path directory = boost::filesystem::temp_directory_path() /
      ("rodin_kokkos_mpi_xdmf_" + std::to_string(id));
    if (g_world->rank() == 0)
      boost::filesystem::create_directory(directory);
    g_world->barrier();

    writeDistributedField(directory, "host", host, hu.getSolution());
    writeDistributedField(directory, "device", mesh, du.getSolution());
    g_world->barrier();

    const std::string rank = std::to_string(g_world->rank());
    for (const char* dataset : {"/Mesh/Geometry/Vertices", "/Mesh/XDMF/Topology"})
    {
      const auto expected =
        readDistributedDataset(directory / ("host.mesh." + rank + ".h5"), dataset);
      const auto actual =
        readDistributedDataset(directory / ("device.mesh." + rank + ".h5"), dataset);
      ASSERT_EQ(actual.size(), expected.size()) << dataset;
      for (size_t i = 0; i < expected.size(); ++i)
        EXPECT_DOUBLE_EQ(actual[i], expected[i]) << dataset << " entry " << i;
    }
    const char* values = "/GridFunction/Values/Data";
    const auto expected =
      readDistributedDataset(directory / ("host." + rank + ".0.h5"), values);
    const auto actual =
      readDistributedDataset(directory / ("device." + rank + ".0.h5"), values);
    ASSERT_EQ(actual.size(), expected.size());
    ASSERT_FALSE(expected.empty());
    double scale = 0;
    for (const double v : expected)
      scale = std::max(scale, std::abs(v));
    ASSERT_GT(scale, 0);
    for (size_t i = 0; i < expected.size(); ++i)
      EXPECT_NEAR(actual[i], expected[i], 1e-9 * scale) << "vertex value " << i;

    g_world->barrier();
    if (g_world->rank() == 0)
      boost::filesystem::remove_all(directory);
  }

  /// @brief A distributed device field is written through its host shard.
  TEST(Kokkos_MPI_IO, XDMFMatchesHost)
  {
    checkDistributedXDMF();
  }
}
