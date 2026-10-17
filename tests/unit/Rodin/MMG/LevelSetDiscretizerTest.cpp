/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Unsupported MMG modes preserve their input and release native resources. */
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  enum class UnsupportedMode
  {
    PlanarSurface,
    EmbeddedSurface,
    EmbeddedMaterial
  };

  class LevelSetOptions : public ::testing::TestWithParam<UnsupportedMode>
  {};

  TEST_P(LevelSetOptions, RejectsUnsupportedModeWithoutPoisoningInput)
  {
    const bool embedded = GetParam() != UnsupportedMode::PlanarSurface;
    constexpr Attribute Material = 17;
    constexpr Attribute Contour = 71;
    constexpr Real Isovalue = 0.37;
    auto planar = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    planar.scale(Real(0.5));
    auto mesh = [&] {
      if (!embedded)
        return planar;
      auto builder = LocalMesh::Builder();
      builder.initialize(3).nodes(planar.getVertexCount());
      for (auto vertex = planar.getVertex(); vertex; ++vertex)
      {
        const auto& coordinates = vertex->getCoordinates();
        builder.vertex({coordinates.x(), coordinates.y(), 0});
      }
      for (auto cell = planar.getCell(); cell; ++cell)
      {
        const auto& vertices = cell->getVertices();
        builder.polytope(
          Polytope::Type::Triangle, {vertices[0], vertices[1], vertices[2]});
      }
      return builder.finalize();
    }();
    mesh.getConnectivity().compute(1, 2);
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute({2, cell->getIndex()}, Material);
    P1 space(mesh);
    MMG::RealGridFunction levelSet(space);
    levelSet = [](const Point& point) { return point.x(); };
    const auto before = levelSet.getData().eval();
    std::vector<Math::SpatialPoint> coordinates;
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      coordinates.push_back(vertex->getCoordinates());

    MMG::LevelSetDiscretizer discretizer;
    discretizer.setLevelSet(Isovalue).setBoundaryReference(Contour).setHMax(0.3);
    const bool material = GetParam() == UnsupportedMode::EmbeddedMaterial;
    if (material)
      discretizer.noSplit(Material);
    else
      discretizer.surface(true);
    try
    {
      discretizer.discretize(levelSet);
      FAIL() << "Unsupported MMG mode was accepted";
    }
    catch (const Alert::Exception& error)
    {
      const std::string message = error.what();
      EXPECT_NE(message.find("member function"), std::string::npos);
      EXPECT_NE(
        message.find(embedded ? "discretizeMMGS" : "discretizeMMG2D"), std::string::npos);
      EXPECT_NE(message.find(material ? "Material splitting is not supported"
                                      : "Meshing the surface"),
        std::string::npos);
    }
    EXPECT_EQ(levelSet.getData(), before);
    ASSERT_EQ(mesh.getVertexCount(), coordinates.size());
    EXPECT_EQ(mesh.getCellCount(), 8u);
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      for (size_t axis = 0; axis < mesh.getSpaceDimension(); ++axis)
      {
        EXPECT_DOUBLE_EQ(
          vertex->getCoordinates()[axis], coordinates[vertex->getIndex()][axis]);
      }
    }
    for (auto cell = mesh.getCell(); cell; ++cell)
      EXPECT_EQ(cell->getAttribute(), Material);

    if (!material)
    {
      // The same object can be retried after correcting its unsupported option.
      discretizer.surface(false);
      const auto result = discretizer.discretize(levelSet);
      ASSERT_FALSE(result.isEmpty());
      Real area = 0;
      Real contourLength = 0;
      for (auto cell = result.getCell(); cell; ++cell)
        area += cell->getMeasure();
      for (auto edge = result.getFace(); edge; ++edge)
      {
        if (edge->getAttribute() != Contour)
          continue;
        contourLength += edge->getMeasure();
        for (const Index vertex : edge->getVertices())
          EXPECT_NEAR(result.getVertexCoordinates(vertex).x(), Isovalue, 1e-10);
      }
      EXPECT_NEAR(area, 1, 1e-10);
      EXPECT_NEAR(contourLength, 1, 1e-10);
      EXPECT_EQ(levelSet.getData(), before);
    }
  }

  INSTANTIATE_TEST_SUITE_P(Unsupported, LevelSetOptions,
    ::testing::Values(UnsupportedMode::PlanarSurface, UnsupportedMode::EmbeddedSurface,
      UnsupportedMode::EmbeddedMaterial));
}
