/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /// @brief Scalar grid function coefficients match the host.
  void checkScalarGridFunction(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    GridFunction hg(hfes);
    fill(hg, source, 1);
    GridFunction dg(dfes);
    dg.getData() = hg.getData();

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });

    expectSameForm(hu, hv, Integral(hg * hu, hv), du, dv, Integral(dg * du, dv));
    expectSameForm(hu, hv, Integral(hg * Grad(hu), Grad(hv)), du, dv,
      Integral(dg * Grad(du), Grad(dv)));
    expectSameForm(
      hu, hv, Integral((hg + hf) * hu, hv), du, dv, Integral((dg + df) * du, dv));
    expectSameForm(
      hu, hv, Integral((hg * hf) * hu, hv), du, dv, Integral((dg * df) * du, dv));
  }

  /// @brief Vector grid function coefficients of the Oseen terms match the host.
  void checkVectorGridFunction(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hvh(source, d);
    P1 hsh(source);
    P1 dvh(mesh, d);
    P1 dsh(mesh);
    TrialFunction hu(hvh);
    TestFunction hv(hvh);
    TrialFunction hp(hsh);
    TrialFunction du(dvh);
    TestFunction dv(dvh);
    TrialFunction dp(dsh);

    GridFunction hw(hvh);
    fill(hw, source, d);
    GridFunction dw(dvh);
    dw.getData() = hw.getData();
    GridFunction hg(hsh);
    fill(hg, source, 1);
    GridFunction dg(dsh);
    dg.getData() = hg.getData();

    // Convection (grad u) w and its Newton counterpart (grad w) u.
    expectSameForm(
      hu, hv, Integral(Jacobian(hu) * hw, hv), du, dv, Integral(Jacobian(du) * dw, dv));
    expectSameForm(
      hu, hv, Integral(Jacobian(hw) * hu, hv), du, dv, Integral(Jacobian(dw) * du, dv));
    expectSameForm(
      hu, hv, Integral(Div(hw) * hu, hv), du, dv, Integral(Div(dw) * du, dv));
    expectSameForm(
      hu, hv, Integral(hw * Div(hu), hv), du, dv, Integral(dw * Div(du), dv));
    expectSameForm(
      hp, hv, Integral(hp * Grad(hg), hv), dp, dv, Integral(dp * Grad(dg), dv));
  }

  /// @brief Component-wise vector function coefficients match the host.
  void checkVectorFunction(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source, 2);
    P1 dfes(mesh, 2);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    VectorFunction hw{
      [](const Geometry::Point& p) { return 1 + p.x() * p.y(); }, Real(0.5)};
    VectorFunction dw(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
        return 1 + p.x() * p.y();
      },
      Real(0.5));
    expectSameForm(
      hu, hv, Integral(Jacobian(hu) * hw, hv), du, dv, Integral(Jacobian(du) * dw, dv));

    VectorFunction hc{2, -1};
    VectorFunction dc{2, -1};
    expectSameForm(
      hu, hv, Integral(Jacobian(hu) * hc, hv), du, dv, Integral(Jacobian(du) * dc, dv));
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, ScalarGridFunctionOnTriangles)
  {
    checkScalarGridFunction(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, ScalarGridFunctionOnMixedMesh)
  {
    checkScalarGridFunction(makeMixedMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, ScalarGridFunctionOnTetrahedra)
  {
    checkScalarGridFunction(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, VectorGridFunctionOnTriangles)
  {
    checkVectorGridFunction(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, VectorGridFunctionOnMixedMesh)
  {
    checkVectorGridFunction(makeMixedMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, VectorGridFunctionOnTetrahedra)
  {
    checkVectorGridFunction(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, VectorFunctionOnTriangles)
  {
    checkVectorFunction(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, VectorFunctionOnMixedMesh)
  {
    checkVectorFunction(makeMixedMesh());
  }
}
