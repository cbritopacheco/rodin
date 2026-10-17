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

  /// @brief Verifies default construction for geometry point cloud by checking exact expected values.
  TEST(Geometry_PointCloud, DefaultConstruction)
  {
    PointCloud pc;
    EXPECT_EQ(pc.getDimension(), 0);
    EXPECT_EQ(pc.getCount(), 0);
    EXPECT_EQ(pc.cols(), 0);
  }

  /// @brief Verifies sized construction 2 D for geometry point cloud by checking exact expected values.
  TEST(Geometry_PointCloud, SizedConstruction2D)
  {
    PointCloud pc(2, 10);
    EXPECT_EQ(pc.getDimension(), 2);
    EXPECT_EQ(pc.rows(), 2);
    EXPECT_EQ(pc.getCount(), 10);
    EXPECT_EQ(pc.cols(), 10);
  }

  /// @brief Verifies sized construction 3 D for geometry point cloud by checking exact expected values.
  TEST(Geometry_PointCloud, SizedConstruction3D)
  {
    PointCloud pc(3, 5);
    EXPECT_EQ(pc.getDimension(), 3);
    EXPECT_EQ(pc.getCount(), 5);
  }

  // ---- Push back ----

  /// @brief Verifies push back 1 D for geometry point cloud by checking tolerance-based numerical results, exact expected values.
  TEST(Geometry_PointCloud, PushBack1D)
  {
    PointCloud pc;
    pc.setDimension(1);
    pc.push_back(std::array<Real, 1>{1.5});
    EXPECT_EQ(pc.getCount(), 1);
    EXPECT_NEAR(pc(0, 0), 1.5, 1e-14);
  }

  /// @brief Verifies push back 2 D for geometry point cloud by checking tolerance-based numerical results, exact expected values.
  TEST(Geometry_PointCloud, PushBack2D)
  {
    PointCloud pc;
    pc.setDimension(2);
    pc.push_back(std::array<Real, 2>{1.0, 2.0});
    EXPECT_EQ(pc.getCount(), 1);
    EXPECT_NEAR(pc(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(pc(1, 0), 2.0, 1e-14);
  }

  /// @brief Verifies push back 3 D for geometry point cloud by checking tolerance-based numerical results, exact expected values.
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

  /// @brief Verifies push back spatial point for geometry point cloud by checking tolerance-based numerical results, exact expected values.
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

  /// @brief Verifies element access read write for geometry point cloud by checking tolerance-based numerical results.
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

  /// @brief Verifies copy construction for geometry point cloud by checking tolerance-based numerical results, exact expected values, copy semantics.
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

  /// @brief Verifies move construction for geometry point cloud by checking tolerance-based numerical results, exact expected values, move semantics.
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

  /// @brief Verifies resize for geometry point cloud by checking exact expected values.
  TEST(Geometry_PointCloud, Resize)
  {
    PointCloud pc;
    pc.resize(3, 10);
    EXPECT_EQ(pc.getDimension(), 3);
    EXPECT_EQ(pc.getCount(), 10);
  }

  /// @brief Verifies clear for geometry point cloud by checking exact expected values.
  TEST(Geometry_PointCloud, Clear)
  {
    PointCloud pc(2, 5);
    EXPECT_EQ(pc.getCount(), 5);
    pc.clear();
    EXPECT_EQ(pc.getCount(), 0);
  }

  // ---- Matrix views ----

  /// @brief Verifies get matrix for geometry point cloud by checking tolerance-based numerical results, exact expected values.
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

  /// @brief Verifies get packed matrix for geometry point cloud by checking tolerance-based numerical results, exact expected values.
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

  /// @brief Verifies multiple push backs for geometry point cloud by checking tolerance-based numerical results, exact expected values.
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

  /// @brief Empty and zero-dimensional coordinate sets preserve shape and inactive packed rows.
  TEST(Geometry_PointCloud, EmptyCoordinateOperationsPreserveShapeAndInactiveRows)
  {
    for (auto dimension :
      {std::uint8_t{0}, std::uint8_t{1}, std::uint8_t{2}, std::uint8_t{3}})
    {
      for (size_t count : {size_t{0}, size_t{3}})
      {
        SCOPED_TRACE(static_cast<unsigned>(dimension));
        SCOPED_TRACE(count);
        PointCloud cloud(dimension, count);
        cloud.getPackedMatrix().setConstant(9);
        cloud.setZero();
        for (size_t point = 0; point < count; ++point)
        {
          EXPECT_EQ(cloud[point].size(), dimension);
          for (size_t row = 0; row < 3; ++row)
            EXPECT_EQ(cloud.getPackedMatrix()(row, point), row < dimension ? 0 : 9);
        }
        cloud.setConstant(5);
        for (size_t point = 0; point < count; ++point)
        {
          for (size_t row = 0; row < 3; ++row)
            EXPECT_EQ(cloud.getPackedMatrix()(row, point), row < dimension ? 5 : 9);
          for (size_t row = 0; row < dimension; ++row)
            cloud(row, point) = static_cast<Real>(1 + 10 * point + row);
        }
        EXPECT_EQ(&(cloud *= Real(2)), &cloud);
        const auto left = Real(3) * cloud;
        const auto right = cloud * Real(3);
        for (const auto* result : {&left, &right})
        {
          EXPECT_EQ(result->getDimension(), dimension);
          EXPECT_EQ(result->getCount(), count);
          EXPECT_EQ(result->getMatrix().rows(), dimension);
          EXPECT_EQ(result->getMatrix().cols(), static_cast<Eigen::Index>(count));
          for (size_t point = 0; point < count; ++point)
          {
            for (size_t row = 0; row < dimension; ++row)
            {
              EXPECT_EQ(
                (*result)(row, point), Real(6) * static_cast<Real>(1 + 10 * point + row));
            }
            if (dimension != 0)
            {
              for (size_t row = dimension; row < 3; ++row)
                EXPECT_EQ(result->getPackedMatrix()(row, point), 9);
            }
          }
        }
        for (size_t point = 0; point < count; ++point)
        {
          for (size_t row = 0; row < 3; ++row)
          {
            EXPECT_EQ(cloud.getPackedMatrix()(row, point),
              row < dimension ? Real(2) * static_cast<Real>(1 + 10 * point + row)
                              : Real(9));
          }
        }
      }
    }
  }

}
