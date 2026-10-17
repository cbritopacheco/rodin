/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @file @brief Physical-node oracles for high-order quadrilateral face orientations.
#include <array>
#include <tuple>
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  class FaceOrientation
    : public testing::TestWithParam<std::tuple<Polytope::Type, size_t>>
  {};

  /// @brief Every square symmetry preserves physical identity between face and volume nodes.
  TEST_P(FaceOrientation, PhysicalNodesShareGlobalDOFs)
  {
    const auto [geometry, orientation] = GetParam();
    const Polytope::Traits traits(geometry);
    LocalMesh::Builder builder;
    builder.initialize(3).nodes(traits.getVertexCount());
    IndexArray vertices(traits.getVertexCount());
    for (Index i = 0; i < traits.getVertexCount(); ++i)
    {
      builder.vertex(traits.getVertex(i));
      vertices[i] = i;
    }
    const std::array<Index, 4> corners = geometry == Polytope::Type::Wedge
      ? std::array<Index, 4>{0, 1, 4, 3}
      : std::array<Index, 4>{0, 1, 2, 3};
    IndexArray ordered(4);
    for (size_t i = 0; i < 4; ++i)
      ordered[i] = corners[(orientation / 2 + (orientation % 2 ? 4 - i : i)) % 4];
    Index faceIndex;
    builder.polytope(Polytope::Type::Quadrilateral, ordered, faceIndex);
    builder.polytope(geometry, vertices);
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    const auto face = mesh.getPolytope(2, faceIndex);
    const auto cell = mesh.getCell();
    ASSERT_TRUE(face);
    ASSERT_TRUE(cell);
    ASSERT_EQ(face->getGeometry(), Polytope::Type::Quadrilateral);
    ASSERT_EQ(face->getVertices().size(), 4u);
    for (size_t i = 0; i < 4; ++i)
      ASSERT_EQ(face->getVertices()[i], ordered[i]);
    ASSERT_EQ(cell->getGeometry(), geometry);
    // Reference GLL nodes are computed to roughly 1e-14; mapping and matching
    // allow 1e-12, far below the separation of distinct degree-three/four nodes.
    constexpr Real coordinateTolerance = 1e-12;
    const auto check = [&]<size_t Degree>() {
      SCOPED_TRACE(Degree);
      H1<Degree, Real, LocalMesh> space(std::integral_constant<size_t, Degree>{}, mesh);
      const auto& faceElement = space.getFiniteElement(2, faceIndex);
      const auto& cellElement = space.getFiniteElement(3, cell->getIndex());
      const auto& faceDOFs = space.getDOFs(2, faceIndex);
      const auto& cellDOFs = space.getDOFs(3, cell->getIndex());
      ASSERT_EQ(faceDOFs.size(), faceElement.getCount());
      ASSERT_EQ(cellDOFs.size(), cellElement.getCount());
      for (size_t localFace = 0; localFace < faceElement.getCount(); ++localFace)
      {
        const Geometry::Point facePoint(*face, faceElement.getNode(localFace));
        const auto expected = facePoint.getPhysicalCoordinates();
        size_t matches = 0;
        for (size_t localCell = 0; localCell < cellElement.getCount(); ++localCell)
        {
          const Geometry::Point cellPoint(*cell, cellElement.getNode(localCell));
          if ((cellPoint.getPhysicalCoordinates() - expected).norm() <
            coordinateTolerance)
          {
            ++matches;
            EXPECT_EQ(faceDOFs[localFace], cellDOFs[localCell]);
          }
        }
        EXPECT_EQ(matches, 1u) << "Face node " << localFace;
      }
    };
    check.template operator()<3>();
    check.template operator()<4>();
  }

  INSTANTIATE_TEST_SUITE_P(AllSquareSymmetries, FaceOrientation,
    testing::Combine(testing::Values(Polytope::Type::Hexahedron, Polytope::Type::Pyramid,
                       Polytope::Type::Wedge),
      testing::Range(size_t(0), size_t(8))));
}
