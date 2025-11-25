/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Variational/H1/Dubiner.h"

using namespace Rodin;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  //==========================================================================
  // DubinerTriangle Collapsed Coordinate Tests
  //==========================================================================

  TEST(DubinerTriangle, CollapsedCoords_Vertex_Origin)
  {
    // (x,y) = (0,0) -> (r,s) = (-1,-1)
    Real r, s;
    DubinerTriangle<2>::getCollapsed(r, s, 0.0, 0.0);
    EXPECT_NEAR(r, -1.0, 1e-14);
    EXPECT_NEAR(s, -1.0, 1e-14);
  }

  TEST(DubinerTriangle, CollapsedCoords_Vertex_X1)
  {
    // (x,y) = (1,0) -> (r,s) = (1,-1)
    Real r, s;
    DubinerTriangle<2>::getCollapsed(r, s, 1.0, 0.0);
    EXPECT_NEAR(r, 1.0, 1e-14);
    EXPECT_NEAR(s, -1.0, 1e-14);
  }

  TEST(DubinerTriangle, CollapsedCoords_Vertex_Y1)
  {
    // (x,y) = (0,1) -> (r,s) = (-1, 1) (collapsed at top)
    Real r, s;
    DubinerTriangle<2>::getCollapsed(r, s, 0.0, 1.0);
    EXPECT_NEAR(r, -1.0, 1e-14);  // Collapsed
    EXPECT_NEAR(s, 1.0, 1e-14);
  }

  TEST(DubinerTriangle, CollapsedCoords_Center)
  {
    // (x,y) = (1/3, 1/3) -> midpoint of triangle
    Real r, s;
    DubinerTriangle<2>::getCollapsed(r, s, 1.0/3.0, 1.0/3.0);

    // s = 2*y - 1 = 2/3 - 1 = -1/3
    // r = 2*x/(1-y) - 1 = 2*(1/3)/(2/3) - 1 = 1 - 1 = 0
    EXPECT_NEAR(s, -1.0/3.0, 1e-14);
    EXPECT_NEAR(r, 0.0, 1e-14);
  }

  //==========================================================================
  // DubinerTriangle Basis Function Tests
  //==========================================================================

  TEST(DubinerTriangle, Basis_P0Q0)
  {
    // ψ_{0,0} = P_0^{0,0}(a) * P_0^{1,0}(b) * 1 = 1 * 1 * 1 = 1
    Real basis;
    DubinerTriangle<2>::getBasis<0, 0>(basis, 0.0, 0.0);
    EXPECT_NEAR(basis, 1.0, 1e-14);

    DubinerTriangle<2>::getBasis<0, 0>(basis, 0.5, -0.5);
    EXPECT_NEAR(basis, 1.0, 1e-14);
  }

  TEST(DubinerTriangle, Basis_P1Q0)
  {
    // ψ_{1,0} involves P_1^{0,0}(a) = a
    Real basis;
    Real r = 0.5, s = -0.5;

    DubinerTriangle<2>::getBasis<1, 0>(basis, r, s);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  TEST(DubinerTriangle, Basis_P0Q1)
  {
    // ψ_{0,1} involves P_0^{0,0}(a) * P_1^{1,0}(b)
    Real basis;
    Real r = 0.0, s = 0.0;

    DubinerTriangle<2>::getBasis<0, 1>(basis, r, s);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  //==========================================================================
  // DubinerTriangle Gradient Tests
  //==========================================================================

  TEST(DubinerTriangle, Gradient_P0Q0)
  {
    // ψ_{0,0} = 1, so ∇ψ = (0, 0)
    Real dpsi_dr, dpsi_ds;
    DubinerTriangle<2>::getGradient<0, 0>(dpsi_dr, dpsi_ds, 0.0, 0.0);
    EXPECT_NEAR(dpsi_dr, 0.0, 1e-14);
    EXPECT_NEAR(dpsi_ds, 0.0, 1e-14);
  }

  TEST(DubinerTriangle, Gradient_P1Q0_NoNaN)
  {
    Real dpsi_dr, dpsi_ds;
    DubinerTriangle<3>::getGradient<1, 0>(dpsi_dr, dpsi_ds, 0.5, -0.5);
    EXPECT_FALSE(std::isnan(dpsi_dr));
    EXPECT_FALSE(std::isnan(dpsi_ds));
  }

  TEST(DubinerTriangle, Gradient_AtSingularity)
  {
    // At s = 1 (top of triangle), there's a singularity
    Real dpsi_dr, dpsi_ds;
    DubinerTriangle<2>::getGradient<1, 0>(dpsi_dr, dpsi_ds, -1.0, 1.0);
    // Should not crash or give NaN
    EXPECT_FALSE(std::isnan(dpsi_dr));
    EXPECT_FALSE(std::isnan(dpsi_ds));
  }

  //==========================================================================
  // VandermondeTriangle Tests
  //==========================================================================

  TEST(VandermondeTriangle, MatrixSize_K1)
  {
    const auto& V = VandermondeTriangle<1>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTriangle<1>::Count);
    EXPECT_EQ(V.cols(), FeketeTriangle<1>::Count);
  }

  TEST(VandermondeTriangle, MatrixSize_K2)
  {
    const auto& V = VandermondeTriangle<2>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTriangle<2>::Count);
    EXPECT_EQ(V.cols(), FeketeTriangle<2>::Count);
  }

  TEST(VandermondeTriangle, MatrixSize_K5)
  {
    const auto& V = VandermondeTriangle<5>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTriangle<5>::Count);
    EXPECT_EQ(V.cols(), FeketeTriangle<5>::Count);
  }

  TEST(VandermondeTriangle, InverseSize_K2)
  {
    const auto& Vinv = VandermondeTriangle<2>::getInverse();
    EXPECT_EQ(Vinv.rows(), FeketeTriangle<2>::Count);
    EXPECT_EQ(Vinv.cols(), FeketeTriangle<2>::Count);
  }

  TEST(VandermondeTriangle, InverseIsInverse_K2)
  {
    const auto& V = VandermondeTriangle<2>::getMatrix();
    const auto& Vinv = VandermondeTriangle<2>::getInverse();

    auto I = V * Vinv;

    for (size_t i = 0; i < V.rows(); ++i)
    {
      for (size_t j = 0; j < V.cols(); ++j)
      {
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(I(i, j), expected, 1e-10);
      }
    }
  }

  TEST(VandermondeTriangle, InverseIsInverse_K5)
  {
    const auto& V = VandermondeTriangle<5>::getMatrix();
    const auto& Vinv = VandermondeTriangle<5>::getInverse();

    auto I = V * Vinv;

    for (size_t i = 0; i < V.rows(); ++i)
    {
      for (size_t j = 0; j < V.cols(); ++j)
      {
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(I(i, j), expected, 1e-9);
      }
    }
  }

  //==========================================================================
  // DubinerTetrahedron Collapsed Coordinate Tests
  //==========================================================================

  TEST(DubinerTetrahedron, CollapsedCoords_Vertex_Origin)
  {
    // (x,y,z) = (0,0,0) -> (a,b,c) = (-1,-1,-1)
    Real a, b, c;
    DubinerTetrahedron<2>::getCollapsed(a, b, c, 0.0, 0.0, 0.0);
    EXPECT_NEAR(a, -1.0, 1e-14);
    EXPECT_NEAR(b, -1.0, 1e-14);
    EXPECT_NEAR(c, -1.0, 1e-14);
  }

  TEST(DubinerTetrahedron, CollapsedCoords_Vertex_X1)
  {
    // (x,y,z) = (1,0,0) -> (a,b,c) = (1,-1,-1)
    Real a, b, c;
    DubinerTetrahedron<2>::getCollapsed(a, b, c, 1.0, 0.0, 0.0);
    EXPECT_NEAR(a, 1.0, 1e-14);
    EXPECT_NEAR(b, -1.0, 1e-14);
    EXPECT_NEAR(c, -1.0, 1e-14);
  }

  TEST(DubinerTetrahedron, CollapsedCoords_Vertex_Y1)
  {
    // (x,y,z) = (0,1,0) -> (a,b,c)
    Real a, b, c;
    DubinerTetrahedron<2>::getCollapsed(a, b, c, 0.0, 1.0, 0.0);
    // c = 2*0 - 1 = -1
    // b = 2*1/(1-0) - 1 = 1
    // a = collapsed since 1-y-z = 0
    EXPECT_NEAR(c, -1.0, 1e-14);
    EXPECT_NEAR(b, 1.0, 1e-14);
    EXPECT_NEAR(a, -1.0, 1e-14);  // collapsed
  }

  TEST(DubinerTetrahedron, CollapsedCoords_Vertex_Z1)
  {
    // (x,y,z) = (0,0,1) -> (a,b,c)
    Real a, b, c;
    DubinerTetrahedron<2>::getCollapsed(a, b, c, 0.0, 0.0, 1.0);
    // c = 2*1 - 1 = 1
    // b = collapsed since z=1
    // a = collapsed
    EXPECT_NEAR(c, 1.0, 1e-14);
    EXPECT_NEAR(b, -1.0, 1e-14);  // collapsed
    EXPECT_NEAR(a, -1.0, 1e-14);  // collapsed
  }

  //==========================================================================
  // DubinerTetrahedron Basis Function Tests
  //==========================================================================

  TEST(DubinerTetrahedron, Basis_P0Q0R0)
  {
    // ψ_{0,0,0} = 1
    Real basis;
    DubinerTetrahedron<2>::getBasis<0, 0, 0>(basis, 0.0, 0.0, 0.0);
    EXPECT_NEAR(basis, 1.0, 1e-14);

    DubinerTetrahedron<2>::getBasis<0, 0, 0>(basis, 0.5, 0.3, -0.2);
    EXPECT_NEAR(basis, 1.0, 1e-14);
  }

  TEST(DubinerTetrahedron, Basis_P1Q0R0_NoNaN)
  {
    Real basis;
    DubinerTetrahedron<2>::getBasis<1, 0, 0>(basis, 0.5, -0.5, -0.5);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  //==========================================================================
  // DubinerTetrahedron Gradient Tests
  //==========================================================================

  TEST(DubinerTetrahedron, Gradient_P0Q0R0)
  {
    // ψ_{0,0,0} = 1, so ∇ψ = (0, 0, 0)
    Real dpsi_da, dpsi_db, dpsi_dc;
    DubinerTetrahedron<2>::getGradient<0, 0, 0>(dpsi_da, dpsi_db, dpsi_dc, 0.0, 0.0, 0.0);
    EXPECT_NEAR(dpsi_da, 0.0, 1e-14);
    EXPECT_NEAR(dpsi_db, 0.0, 1e-14);
    EXPECT_NEAR(dpsi_dc, 0.0, 1e-14);
  }

  TEST(DubinerTetrahedron, Gradient_P1Q0R0_NoNaN)
  {
    Real dpsi_da, dpsi_db, dpsi_dc;
    DubinerTetrahedron<3>::getGradient<1, 0, 0>(dpsi_da, dpsi_db, dpsi_dc, 0.5, -0.5, -0.5);
    EXPECT_FALSE(std::isnan(dpsi_da));
    EXPECT_FALSE(std::isnan(dpsi_db));
    EXPECT_FALSE(std::isnan(dpsi_dc));
  }

  //==========================================================================
  // VandermondeTetrahedron Tests
  //==========================================================================

  TEST(VandermondeTetrahedron, MatrixSize_K1)
  {
    const auto& V = VandermondeTetrahedron<1>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTetrahedron<1>::Count);
    EXPECT_EQ(V.cols(), FeketeTetrahedron<1>::Count);
  }

  TEST(VandermondeTetrahedron, MatrixSize_K2)
  {
    const auto& V = VandermondeTetrahedron<2>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTetrahedron<2>::Count);
    EXPECT_EQ(V.cols(), FeketeTetrahedron<2>::Count);
  }

  TEST(VandermondeTetrahedron, InverseSize_K2)
  {
    const auto& Vinv = VandermondeTetrahedron<2>::getInverse();
    EXPECT_EQ(Vinv.rows(), FeketeTetrahedron<2>::Count);
    EXPECT_EQ(Vinv.cols(), FeketeTetrahedron<2>::Count);
  }

  TEST(VandermondeTetrahedron, InverseIsInverse_K2)
  {
    const auto& V = VandermondeTetrahedron<2>::getMatrix();
    const auto& Vinv = VandermondeTetrahedron<2>::getInverse();

    auto I = V * Vinv;

    for (size_t i = 0; i < V.rows(); ++i)
    {
      for (size_t j = 0; j < V.cols(); ++j)
      {
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(I(i, j), expected, 1e-10);
      }
    }
  }

  //==========================================================================
  // Higher Order Tests (K = 5, 6)
  //==========================================================================

  TEST(DubinerTriangle, Basis_P2Q1_NoNaN)
  {
    Real basis;
    Real r = 0.3, s = -0.3;

    DubinerTriangle<5>::getBasis<2, 1>(basis, r, s);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  TEST(DubinerTriangle, Basis_P3Q2_NoNaN)
  {
    Real basis;
    Real r = 0.0, s = 0.5;

    DubinerTriangle<6>::getBasis<3, 2>(basis, r, s);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  TEST(DubinerTriangle, Gradient_P2Q1_NoNaN)
  {
    Real dpsi_dr, dpsi_ds;
    DubinerTriangle<5>::getGradient<2, 1>(dpsi_dr, dpsi_ds, 0.3, -0.5);
    EXPECT_FALSE(std::isnan(dpsi_dr));
    EXPECT_FALSE(std::isnan(dpsi_ds));
    EXPECT_FALSE(std::isinf(dpsi_dr));
    EXPECT_FALSE(std::isinf(dpsi_ds));
  }

  TEST(VandermondeTriangle, InverseIsInverse_K6)
  {
    const auto& V = VandermondeTriangle<6>::getMatrix();
    const auto& Vinv = VandermondeTriangle<6>::getInverse();

    auto I = V * Vinv;

    for (size_t i = 0; i < V.rows(); ++i)
    {
      for (size_t j = 0; j < V.cols(); ++j)
      {
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(I(i, j), expected, 1e-8);
      }
    }
  }

  TEST(VandermondeTriangle, MatrixSize_K6)
  {
    const auto& V = VandermondeTriangle<6>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTriangle<6>::Count);
    EXPECT_EQ(V.cols(), FeketeTriangle<6>::Count);
  }

  TEST(DubinerTetrahedron, Basis_P1Q1R0_NoNaN)
  {
    Real basis;
    DubinerTetrahedron<5>::getBasis<1, 1, 0>(basis, 0.3, -0.3, -0.5);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  TEST(DubinerTetrahedron, Basis_P2Q1R1_NoNaN)
  {
    Real basis;
    DubinerTetrahedron<6>::getBasis<2, 1, 1>(basis, 0.0, 0.2, -0.5);
    EXPECT_FALSE(std::isnan(basis));
    EXPECT_FALSE(std::isinf(basis));
  }

  TEST(DubinerTetrahedron, Gradient_P1Q1R0_NoNaN)
  {
    Real dpsi_da, dpsi_db, dpsi_dc;
    DubinerTetrahedron<5>::getGradient<1, 1, 0>(dpsi_da, dpsi_db, dpsi_dc, 0.3, -0.3, -0.5);
    EXPECT_FALSE(std::isnan(dpsi_da));
    EXPECT_FALSE(std::isnan(dpsi_db));
    EXPECT_FALSE(std::isnan(dpsi_dc));
  }

  TEST(VandermondeTetrahedron, InverseIsInverse_K5)
  {
    const auto& V = VandermondeTetrahedron<5>::getMatrix();
    const auto& Vinv = VandermondeTetrahedron<5>::getInverse();

    auto I = V * Vinv;

    for (size_t i = 0; i < V.rows(); ++i)
    {
      for (size_t j = 0; j < V.cols(); ++j)
      {
        Real expected = (i == j) ? 1.0 : 0.0;
        EXPECT_NEAR(I(i, j), expected, 1e-8);
      }
    }
  }

  TEST(VandermondeTetrahedron, MatrixSize_K5)
  {
    const auto& V = VandermondeTetrahedron<5>::getMatrix();
    EXPECT_EQ(V.rows(), FeketeTetrahedron<5>::Count);
    EXPECT_EQ(V.cols(), FeketeTetrahedron<5>::Count);
  }
}
