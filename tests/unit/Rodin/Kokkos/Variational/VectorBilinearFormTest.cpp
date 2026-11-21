/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  static_assert(std::is_same_v<decltype(P1(std::declval<const DeviceMesh&>(), size_t(2))),
    P1<Math::SpatialVector<Real>, DeviceMesh>>);

  /// @brief Matrix coefficients contracted with Jacobians match the host.
  void checkMatrixCoefficient(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    Math::SpatialMatrix<Real> a(d, d);
    for (size_t i = 0; i < d; ++i)
    {
      for (size_t j = 0; j < d; ++j)
        a(i, j) = 1 + Real(i) - 2 * Real(j) + Real(i == j);
    }
    MatrixFunction ha(a);
    MatrixFunction da(a);
    expectSameForm(hu, hv, Integral(Dot(ha, Jacobian(hu)), Dot(ha, Jacobian(hv))), du, dv,
      Integral(Dot(da, Jacobian(du)), Dot(da, Jacobian(dv))));

    GridFunction hw(hfes);
    GridFunction dw(dfes);
    fill(hw, source, d);
    fill(dw, source, d);
    expectSameForm(hu, hv,
      Integral(Dot(Jacobian(hw), Jacobian(hu)), Dot(ha, Jacobian(hv))), du, dv,
      Integral(Dot(Jacobian(dw), Jacobian(du)), Dot(da, Jacobian(dv))));
  }

  /// @brief Velocity and pressure coupling forms match the host.
  void checkMixedForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hvh(source, d);
    P1 hph(source);
    P1 dvh(mesh, d);
    P1 dph(mesh);
    TrialFunction hu(hvh);
    TestFunction hv(hvh);
    TrialFunction hp(hph);
    TestFunction hq(hph);
    TrialFunction du(dvh);
    TestFunction dv(dvh);
    TrialFunction dp(dph);
    TestFunction dq(dph);

    expectSameForm(hu, hq, Integral(Div(hu), hq), du, dq, Integral(Div(du), dq));
    expectSameForm(hp, hv, Integral(hp, Div(hv)), dp, dv, Integral(dp, Div(dv)));
  }

  /// @brief Block diagonal integrals can list the diagonals of their blocks only.
  void checkBlockDiagonal(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 fes(mesh, d);
    TrialFunction u(fes);
    TestFunction v(fes);
    auto mass = Integral(u, v);
    auto jacobian = Integral(Jacobian(u), Jacobian(v));
    auto divergence = Integral(Div(u), Div(v));
    EXPECT_TRUE(mass.isBlockDiagonal());
    EXPECT_TRUE(jacobian.isBlockDiagonal());
    EXPECT_FALSE(divergence.isBlockDiagonal());
    const size_t full = mass.prepare(false);
    EXPECT_GT(full, 0);
    EXPECT_EQ(mass.prepare(true) * d, full);
    EXPECT_EQ(jacobian.prepare(true) * d, full);
    EXPECT_EQ(divergence.prepare(true), full);
  }

  TEST(Kokkos_Variational_VectorBilinearForm, BlockDiagonalOnTriangles)
  {
    checkBlockDiagonal(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, BlockDiagonalOnTetrahedra)
  {
    checkBlockDiagonal(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, MatrixCoefficientOnTriangles)
  {
    checkMatrixCoefficient(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, MatrixCoefficientOnMixedMesh)
  {
    checkMatrixCoefficient(makeMixedMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, MatrixCoefficientOnTetrahedra)
  {
    checkMatrixCoefficient(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, MixedFormsOnTriangles)
  {
    checkMixedForms(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, MixedFormsOnMixedMesh)
  {
    checkMixedForms(makeMixedMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, MixedFormsOnTetrahedra)
  {
    checkMixedForms(makeTetrahedronMesh());
  }
}
