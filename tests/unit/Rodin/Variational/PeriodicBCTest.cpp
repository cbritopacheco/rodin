/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Partition-of-unity and lifecycle contracts for periodic constraints. */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Variational/PeriodicBC.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Scalar>
  void checkPeriodicConstraints()
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    P1<Scalar> space(mesh);
    TrialFunction u(space);
    const IndexMap<IndexSet> adjacency{{0, {1, 2, 3}}, {1, {2}}};
    PeriodicBC boundary(u, adjacency);
    EXPECT_EQ(&boundary.getOperand(), &u);
    EXPECT_EQ(boundary.getAdjacency(), adjacency);
    EXPECT_FALSE(boundary.isComponent());
    EXPECT_TRUE(boundary.getDOFs().empty());
    const auto verify = [&adjacency](const auto& constraints) {
      ASSERT_EQ(constraints.getDOFs().size(), adjacency.size());
      for (const auto& [parent, children] : adjacency)
      {
        const auto& [indices, coefficients] = constraints.getDOFs().at(parent);
        ASSERT_EQ(static_cast<size_t>(indices.size()), children.size());
        ASSERT_EQ(coefficients.size(), indices.size());
        size_t position = 0;
        Scalar constantValue = 0;
        for (const auto child : children)
        {
          EXPECT_EQ(indices[position], child);
          EXPECT_NEAR(
            std::abs(coefficients[position] - Scalar(1) / Scalar(children.size())), 0,
            1e-15);
          constantValue += coefficients[position] * Scalar(7);
          ++position;
        }
        EXPECT_NEAR(std::abs(constantValue - Scalar(7)), 0, 1e-14);
      }
    };
    boundary.assemble();
    verify(boundary);
    boundary.assemble();
    verify(boundary);
    std::unique_ptr<decltype(boundary)> clone(boundary.copy());
    verify(*clone);
    clone->assemble();
    verify(*clone);
    auto moved = std::move(*clone);
    verify(moved);
    moved.assemble();
    verify(moved);
    PeriodicBC temporary(u, IndexMap<IndexSet>(adjacency));
    temporary.assemble();
    verify(temporary);
    PeriodicBC empty(u, IndexMap<IndexSet>{});
    empty.assemble();
    EXPECT_TRUE(empty.getDOFs().empty());
  }
}

TEST(Rodin_Variational_PeriodicBC, RealConstraintsPreserveConstants)
{
  checkPeriodicConstraints<Real>();
}

TEST(Rodin_Variational_PeriodicBC, ComplexConstraintsPreserveConstants)
{
  checkPeriodicConstraints<Complex>();
}
