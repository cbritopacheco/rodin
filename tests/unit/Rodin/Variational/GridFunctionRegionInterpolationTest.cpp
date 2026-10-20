/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE, or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <type_traits>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

template <class Scalar>
class GridFunctionRegionInterpolation : public ::testing::Test
{};

using RegionInterpolationScalars = ::testing::Types<Real, Complex>;
TYPED_TEST_SUITE(GridFunctionRegionInterpolation, RegionInterpolationScalars);

/// @brief Facet interpolation reproduces an affine trace, preserves unselected vertices, and is idempotent in P1/P2 scalar/vector spaces.
TYPED_TEST(GridFunctionRegionInterpolation, SelectiveFacesAndInterfaces)
{
  using Scalar = TypeParam;
  constexpr Real tolerance = 1e-12;
  const Scalar factor = [] {
    if constexpr (std::is_same_v<Scalar, Complex>)
      return Scalar(1, -2);
    else
      return Scalar(1);
  }();
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(1, 2);
  const auto check = [&](const auto& space) {
    using Range = typename FormLanguage::Traits<std::decay_t<decltype(space)>>::RangeType;
    const auto source = [&](const Point& point) -> Range {
      if constexpr (std::is_same_v<Range, Scalar>)
        return factor * (2 + 3 * point.x() + 5 * point.y());
      else
      {
        return Range{factor * (2 + 3 * point.x() + 5 * point.y()),
          factor * (-1 + 7 * point.x() - 2 * point.y())};
      }
    };
    const Range initial = [&]() -> Range {
      if constexpr (std::is_same_v<Range, Scalar>)
        return Scalar(-11) * factor;
      else
        return Range{Scalar(-11) * factor, Scalar(13) * factor};
    }();
    const auto expectValue = [&](const Range& value, const Range& expected) {
      if constexpr (std::is_same_v<Range, Scalar>)
        EXPECT_LE(std::abs(value - expected), tolerance);
      else
        EXPECT_LE((value - expected).norm(), tolerance);
    };
    for (Region region : {Region::Faces, Region::Interface})
    {
      SCOPED_TRACE(static_cast<int>(region));
      auto selected = region == Region::Faces ? mesh.getBoundary() : mesh.getInterface();
      ASSERT_TRUE(selected);
      const Index faceIndex = selected->getIndex();
      const auto& vertices = selected->getVertices();
      GridFunction field(space);
      field = initial;
      const auto untouched = field.getData().eval();
      field.project(region, source, [](const Polytope&) { return false; });
      EXPECT_EQ(field.getData(), untouched);
      const auto predicate = [&](const Polytope& polytope) {
        return polytope.getDimension() == 1 && polytope.getIndex() == faceIndex;
      };
      for (size_t repetition = 0; repetition < 2; ++repetition)
      {
        field.project(region, source, predicate);
        for (Real coordinate : {Real(0), Real(0.25), Real(0.5), Real(0.75), Real(1)})
        {
          const Point point(*selected, Math::SpatialPoint{coordinate});
          expectValue(field.getValue(point), source(point));
        }
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          const Point point(*mesh.getPolytope(0, vertex), Math::SpatialPoint{});
          const bool onSelectedFace =
            std::find(vertices.begin(), vertices.end(), vertex) != vertices.end();
          expectValue(field.getValue(point), onSelectedFace ? source(point) : initial);
        }
      }
    }
  };
  check(P1<Scalar>(mesh));
  check(P1<Math::SpatialVector<Scalar>>(mesh, 2));
  check(H1<2, Scalar>(std::integral_constant<size_t, 2>{}, mesh));
  check(H1<2, Math::SpatialVector<Scalar>>(std::integral_constant<size_t, 2>{}, mesh, 2));
}

/// @brief A selected submesh evaluates retained root points and rejects other root cells without modifying its field.
TYPED_TEST(GridFunctionRegionInterpolation, SubMeshPointOwnership)
{
  using Scalar = TypeParam;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  ASSERT_EQ(mesh.getCellCount(), 2);
  mesh.getConnectivity().compute(1, 2);
  mesh.setAttribute({2, 0}, Attribute(7));
  mesh.setAttribute({2, 1}, Attribute(8));
  auto submesh = mesh.keep(Attribute(7));
  submesh.getConnectivity().compute(1, 0);
  const auto coordinates = Polytope::Traits(Polytope::Type::Triangle).getCentroid();
  const Point retained(*mesh.getCell(0), coordinates);
  const Point rejected(*mesh.getCell(1), coordinates);
  const auto check = [&](const auto& space) {
    GridFunction field(space);
    field = Scalar(3);
    const auto original = field.getData().eval();
    EXPECT_LE(std::abs(field.getValue(retained) - Scalar(3)), Real(1e-12));
    EXPECT_THROW(field.getValue(rejected), Alert::Exception);
    EXPECT_EQ(field.getData(), original);
    EXPECT_LE(std::abs(field.getValue(retained) - Scalar(3)), Real(1e-12));
  };
  check(P1<Scalar>(submesh));
  check(H1<2, Scalar>(std::integral_constant<size_t, 2>{}, submesh));
}
