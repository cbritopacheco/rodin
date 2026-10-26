/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
/// @file @brief Native HDF5 persistence preserves independently owned geometry charts.
#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include <fstream>
#include <vector>
#include "Rodin/Geometry/IdentityTransformation.h"

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/IO/HDF5.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::IO;

namespace
{
  class HDF5Geometry : public testing::Test
  {
    protected:
      void SetUp() override
      {
        directory = boost::filesystem::temp_directory_path() /
          boost::filesystem::unique_path("rodin-geometry-%%%%-%%%%");
        boost::filesystem::create_directory(directory);
      }

      void TearDown() override
      {
        boost::filesystem::remove_all(directory);
      }
      boost::filesystem::path directory;
  };

  template <size_t Degree>
  void roundtrip(const boost::filesystem::path& file, bool embedded)
  {
    auto mesh = [&] {
      if (embedded)
      {
        return LocalMesh::Builder()
          .initialize(3)
          .nodes(3)
          .vertex({0, 0, 0})
          .vertex({1.1, 0, 0})
          .vertex({0, 1, 0})
          .polytope(Polytope::Type::Triangle, {0, 1, 2})
          .finalize();
      }
      return LocalMesh::Builder()
        .initialize(2)
        .nodes(3)
        .vertex({0, 0})
        .vertex({1.1, 0})
        .vertex({0, 1})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .finalize();
    }();
    RealH1Element<Degree> element(Polytope::Type::Triangle);
    PointCloud points(mesh.getSpaceDimension(), element.getCount());
    for (size_t local = 0; local < element.getCount(); ++local)
    {
      const auto r = element.getNode(local);
      points(0, local) = r[0] + 0.1 * r[0] * r[0];
      points(1, local) = r[1] + 0.2 * r[0] * r[1];
      if (embedded)
        points(2, local) = 0.3 * r[0] * r[1];
    }
    mesh.setPolytopeTransformation({2, 0},
      new ParametricTransformation<RealH1Element<Degree>>(
        std::move(points), std::move(element)));
    MeshPrinter<FileFormat::HDF5, Context::Local>(mesh).print(file);
    LocalMesh loaded;
    MeshLoader<FileFormat::HDF5, Context::Local>(loaded).load(file);
    const auto& source = mesh.getCell()->getTransformation();
    const auto& restored = loaded.getCell()->getTransformation();
    EXPECT_EQ(restored.getFactorOrder(), Degree);
    for (const auto& reference :
      {Math::SpatialPoint{0.2, 0.3}, Math::SpatialPoint{0.6, 0.1}})
    {
      Math::SpatialPoint original, actual;
      Math::SpatialMatrix<Real> originalJacobian, actualJacobian;
      source.transform(original, reference);
      restored.transform(actual, reference);
      source.jacobian(originalJacobian, reference);
      restored.jacobian(actualJacobian, reference);
      EXPECT_LT((actual - original).norm(), 1e-13);
      EXPECT_LT((actualJacobian - originalJacobian).norm(), 1e-12);
    }
  }
}

TEST_F(HDF5Geometry, CurvedAndEmbeddedFamilies)
{
  const auto file = directory / "curved.h5";
  const auto check = [&]<size_t Degree>() {
    roundtrip<Degree>(file, false);
    roundtrip<Degree>(file, true);
  };
  check.template operator()<1>();
  check.template operator()<2>();
  check.template operator()<3>();
  check.template operator()<4>();
  check.template operator()<5>();
  check.template operator()<6>();
}

TEST_F(HDF5Geometry, ExplicitP1ChartAndMalformedFamily)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  RealP1Element element(Polytope::Type::Segment);
  PointCloud points(1, 2);
  points(0, 0) = 2;
  points(0, 1) = 4;
  mesh.setPolytopeTransformation({1, 0},
    new ParametricTransformation<RealP1Element>(std::move(points), std::move(element)));
  const auto file = directory / "p1.h5";
  MeshPrinter<FileFormat::HDF5, Context::Local>(mesh).print(file);
  LocalMesh restored;
  MeshLoader<FileFormat::HDF5, Context::Local> loader(restored);
  loader.load(file);
  Math::SpatialPoint physical;
  restored.getCell()->getTransformation().transform(physical, Math::SpatialPoint{0.5});
  EXPECT_EQ(physical[0], 3);
  const auto before = restored.getVertexCoordinates(0);
  {
    const auto handle = HDF5::File(H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    ASSERT_TRUE(handle);
    const auto dataset =
      HDF5::DataSet(H5Dopen2(handle.get(), "/Mesh/Transformations/1/Kind", H5P_DEFAULT));
    const HDF5::I32 invalid = -1;
    ASSERT_GE(
      H5Dwrite(dataset.get(), H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &invalid),
      0);
  }
  EXPECT_THROW(loader.load(file), Alert::Exception);
  EXPECT_EQ(restored.getVertexCoordinates(0).getData(), before.getData());
}

/// @brief Unregistered chart families fail before truncating an existing checkpoint.
TEST_F(HDF5Geometry, UnsupportedChartPreservesExistingFile)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
  mesh.setPolytopeTransformation({1, 0}, new IdentityTransformation(1));
  const auto file = directory / "existing.h5";
  {
    std::ofstream output(file.string());
    output << "retained checkpoint";
  }
  EXPECT_THROW(
    (MeshPrinter<FileFormat::HDF5, Context::Local>(mesh).print(file)), Alert::Exception);
  std::ifstream input(file.string());
  const std::string contents((std::istreambuf_iterator<char>(input)), {});
  EXPECT_EQ(contents, "retained checkpoint");
}

/// @brief XDMF mixed topology encodes the node count required by Polyvertex cells.
TEST_F(HDF5Geometry, XDMFPointTopology)
{
  const auto mesh =
    LocalMesh::Build().initialize(1).nodes(2).vertex({0}).vertex({1}).finalize();
  const auto path = directory / "points.h5";
  HDF5::writeXDMFMesh(path.string(), mesh);
  auto file = HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  const auto topology =
    HDF5::readVectorDataset<HDF5::U64>(file.get(), HDF5::Path::MeshXDMFTopology);
  const std::vector<HDF5::U64> expected{1, 1, 0, 1, 1, 1};
  EXPECT_EQ(topology, expected);
}
