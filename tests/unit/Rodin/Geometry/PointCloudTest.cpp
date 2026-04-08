/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry/PointCloud.h>
#include <Rodin/Math/SpatialVector.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- Construction ----

  TEST(Geometry_PointCloud, DefaultConstruction)
  {
    PointCloud pc;
    EXPECT_EQ(pc.getDimension(), 0);
    EXPECT_EQ(pc.getCount(), 0);
    EXPECT_EQ(pc.cols(), 0);
  }

  TEST(Geometry_PointCloud, SizedConstruction2D)
  {
    PointCloud pc(2, 10);
    EXPECT_EQ(pc.getDimension(), 2);
    EXPECT_EQ(pc.rows(), 2);
    EXPECT_EQ(pc.getCount(), 10);
    EXPECT_EQ(pc.cols(), 10);
  }

  TEST(Geometry_PointCloud, SizedConstruction3D)
  {
    PointCloud pc(3, 5);
    EXPECT_EQ(pc.getDimension(), 3);
    EXPECT_EQ(pc.getCount(), 5);
  }

  // ---- Push back ----

  TEST(Geometry_PointCloud, PushBack1D)
  {
    PointCloud pc;
    pc.setDimension(1);
    pc.push_back(std::array<Real, 1>{1.5});
    EXPECT_EQ(pc.getCount(), 1);
    EXPECT_NEAR(pc(0, 0), 1.5, 1e-14);
  }

  TEST(Geometry_PointCloud, PushBack2D)
  {
    PointCloud pc;
    pc.setDimension(2);
    pc.push_back(std::array<Real, 2>{1.0, 2.0});
    EXPECT_EQ(pc.getCount(), 1);
    EXPECT_NEAR(pc(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(pc(1, 0), 2.0, 1e-14);
  }

  TEST(Geometry_PointCloud, PushBack3D)
  {
    PointCloud pc;
    pc.setDimension(3);
    pc.push_back(std::array<Real, 3>{1.0, 2.0, 3.0});
    EXPECT_EQ(pc.getCount(), 1);
    EXPECT_NEAR(pc(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(pc(1, 0), 2.0, 1e-14);
    EXPECT_NEAR(pc(2, 0), 3.0, 1e-14);
  }

  TEST(Geometry_PointCloud, PushBackSpatialPoint)
  {
    PointCloud pc;
    pc.setDimension(2);
    Math::SpatialPoint sp{3.0, 4.0};
    pc.push_back(sp);
    EXPECT_EQ(pc.getCount(), 1);
    EXPECT_NEAR(pc(0, 0), 3.0, 1e-14);
    EXPECT_NEAR(pc(1, 0), 4.0, 1e-14);
  }

  // ---- Element access ----

  TEST(Geometry_PointCloud, ElementAccessReadWrite)
  {
    PointCloud pc(2, 3);
    pc(0, 0) = 1.0; pc(1, 0) = 2.0;
    pc(0, 1) = 3.0; pc(1, 1) = 4.0;
    pc(0, 2) = 5.0; pc(1, 2) = 6.0;

    EXPECT_NEAR(pc(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(pc(1, 0), 2.0, 1e-14);
    EXPECT_NEAR(pc(0, 1), 3.0, 1e-14);
    EXPECT_NEAR(pc(1, 1), 4.0, 1e-14);
    EXPECT_NEAR(pc(0, 2), 5.0, 1e-14);
    EXPECT_NEAR(pc(1, 2), 6.0, 1e-14);
  }

  // ---- Copy and move ----

  TEST(Geometry_PointCloud, CopyConstruction)
  {
    PointCloud pc(2, 2);
    pc(0, 0) = 1.0; pc(1, 0) = 2.0;
    pc(0, 1) = 3.0; pc(1, 1) = 4.0;

    PointCloud copy(pc);
    EXPECT_EQ(copy.getDimension(), 2);
    EXPECT_EQ(copy.getCount(), 2);
    EXPECT_NEAR(copy(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(copy(1, 1), 4.0, 1e-14);
  }

  TEST(Geometry_PointCloud, MoveConstruction)
  {
    PointCloud pc(2, 2);
    pc(0, 0) = 1.0; pc(1, 0) = 2.0;

    PointCloud moved(std::move(pc));
    EXPECT_EQ(moved.getDimension(), 2);
    EXPECT_EQ(moved.getCount(), 2);
    EXPECT_NEAR(moved(0, 0), 1.0, 1e-14);
  }

  // ---- Resize and reserve ----

  TEST(Geometry_PointCloud, Resize)
  {
    PointCloud pc;
    pc.resize(3, 10);
    EXPECT_EQ(pc.getDimension(), 3);
    EXPECT_EQ(pc.getCount(), 10);
  }

  TEST(Geometry_PointCloud, Clear)
  {
    PointCloud pc(2, 5);
    EXPECT_EQ(pc.getCount(), 5);
    pc.clear();
    EXPECT_EQ(pc.getCount(), 0);
  }

  // ---- Matrix views ----

  TEST(Geometry_PointCloud, GetMatrix)
  {
    PointCloud pc(2, 3);
    pc(0, 0) = 1.0; pc(1, 0) = 2.0;
    pc(0, 1) = 3.0; pc(1, 1) = 4.0;
    pc(0, 2) = 5.0; pc(1, 2) = 6.0;

    auto mat = pc.getMatrix();
    EXPECT_EQ(mat.rows(), 2);
    EXPECT_EQ(mat.cols(), 3);
    EXPECT_NEAR(mat(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(mat(1, 2), 6.0, 1e-14);
  }

  TEST(Geometry_PointCloud, GetPackedMatrix)
  {
    PointCloud pc(2, 2);
    pc(0, 0) = 1.0; pc(1, 0) = 2.0;
    pc(0, 1) = 3.0; pc(1, 1) = 4.0;

    auto pmat = pc.getPackedMatrix();
    EXPECT_EQ(pmat.rows(), 3);
    EXPECT_EQ(pmat.cols(), 2);
    EXPECT_NEAR(pmat(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(pmat(1, 0), 2.0, 1e-14);
  }

  // ---- Multiple push_backs ----

  TEST(Geometry_PointCloud, MultiplePushBacks)
  {
    PointCloud pc;
    pc.setDimension(2);
    for (int i = 0; i < 100; ++i)
      pc.push_back(std::array<Real, 2>{static_cast<Real>(i), static_cast<Real>(i * 2)});

    EXPECT_EQ(pc.getCount(), 100);
    EXPECT_NEAR(pc(0, 50), 50.0, 1e-14);
    EXPECT_NEAR(pc(1, 50), 100.0, 1e-14);
  }
}
