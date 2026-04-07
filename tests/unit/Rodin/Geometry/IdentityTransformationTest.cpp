/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- IdentityTransformation basic properties ----

  TEST(Geometry_IdentityTransformation, Construction)
  {
    IdentityTransformation t(2);
    EXPECT_EQ(t.getReferenceDimension(), 2);
    EXPECT_EQ(t.getPhysicalDimension(), 2);
    EXPECT_EQ(t.getOrder(), 1);
  }

  TEST(Geometry_IdentityTransformation, Construction3D)
  {
    IdentityTransformation t(3);
    EXPECT_EQ(t.getReferenceDimension(), 3);
    EXPECT_EQ(t.getPhysicalDimension(), 3);
  }

  TEST(Geometry_IdentityTransformation, Transform2D)
  {
    IdentityTransformation t(2);
    Math::SpatialPoint rc{0.3, 0.7};
    Math::SpatialPoint pc(2);
    t.transform(pc, rc);
    EXPECT_NEAR(pc(0), 0.3, 1e-14);
    EXPECT_NEAR(pc(1), 0.7, 1e-14);
  }

  TEST(Geometry_IdentityTransformation, Transform3D)
  {
    IdentityTransformation t(3);
    Math::SpatialPoint rc{0.1, 0.2, 0.3};
    Math::SpatialPoint pc(3);
    t.transform(pc, rc);
    EXPECT_NEAR(pc(0), 0.1, 1e-14);
    EXPECT_NEAR(pc(1), 0.2, 1e-14);
    EXPECT_NEAR(pc(2), 0.3, 1e-14);
  }

  TEST(Geometry_IdentityTransformation, Jacobian2D)
  {
    IdentityTransformation t(2);
    Math::SpatialPoint rc{0.5, 0.5};
    Math::SpatialMatrix<Real> jac(2, 2);
    t.jacobian(jac, rc);
    // Should be identity matrix
    EXPECT_NEAR(jac(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(jac(0, 1), 0.0, 1e-14);
    EXPECT_NEAR(jac(1, 0), 0.0, 1e-14);
    EXPECT_NEAR(jac(1, 1), 1.0, 1e-14);
  }

  TEST(Geometry_IdentityTransformation, Jacobian3D)
  {
    IdentityTransformation t(3);
    Math::SpatialPoint rc{0.1, 0.2, 0.3};
    Math::SpatialMatrix<Real> jac(3, 3);
    t.jacobian(jac, rc);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        EXPECT_NEAR(jac(i, j), (i == j) ? 1.0 : 0.0, 1e-14);
  }

  TEST(Geometry_IdentityTransformation, Copy)
  {
    IdentityTransformation t(2);
    auto* c = t.copy();
    EXPECT_EQ(c->getReferenceDimension(), 2);
    EXPECT_EQ(c->getPhysicalDimension(), 2);
    EXPECT_EQ(c->getOrder(), 1);

    Math::SpatialPoint rc{0.4, 0.6};
    Math::SpatialPoint pc(2);
    c->transform(pc, rc);
    EXPECT_NEAR(pc(0), 0.4, 1e-14);
    EXPECT_NEAR(pc(1), 0.6, 1e-14);
    delete c;
  }

  TEST(Geometry_IdentityTransformation, MoveConstruction)
  {
    IdentityTransformation t(3);
    IdentityTransformation moved(std::move(t));
    EXPECT_EQ(moved.getReferenceDimension(), 3);
    EXPECT_EQ(moved.getPhysicalDimension(), 3);
  }

  TEST(Geometry_IdentityTransformation, CopyConstruction)
  {
    IdentityTransformation t(2);
    IdentityTransformation copy(t);
    EXPECT_EQ(copy.getReferenceDimension(), 2);
    EXPECT_EQ(copy.getPhysicalDimension(), 2);
  }
}
