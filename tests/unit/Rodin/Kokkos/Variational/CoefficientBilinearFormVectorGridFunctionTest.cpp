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
}
