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

  /**
   * @brief Coefficients changed between assemblies reach the device.
   *
   * Two grid functions of the same space enter one form, their device
   * copies being reused while their coefficients do not change. The
   * coefficients are then changed, swapped and restored between
   * assemblies of the same form.
   */
  void checkGridFunctionChange(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    GridFunction hf(hfes);
    GridFunction hg(hfes);
    fill(hf, source, 1);
    hg.getData() = hf.getData().array() * Real(-0.5) + Real(2);
    GridFunction df(dfes);
    GridFunction dg(dfes);
    df.getData() = hf.getData();
    dg.getData() = hg.getData();

    BilinearForm ha(hu, hv);
    ha = Integral(hf * hu, hv) + Integral((hf * hg) * Grad(hu), Grad(hv));
    BilinearForm da(du, dv);
    da = Integral(df * du, dv) + Integral((df * dg) * Grad(du), Grad(dv));
    const auto check = [&]() {
      df.getData() = hf.getData();
      dg.getData() = hg.getData();
      ha.assemble();
      da.assemble();
      expectSameOperator(ha.getOperator(), da.getOperator());
    };
    check();
    check();
    const Math::Vector<Real> first = hf.getData();
    hf.getData() = hf.getData().array() + Real(1);
    check();
    hg.getData() = first;
    check();
    hf.getData() = first;
    hg.getData() = first.array() * Real(-0.5) + Real(2);
    check();
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

  TEST(Kokkos_Variational_CoefficientBilinearForm, GridFunctionChangeOnTriangles)
  {
    checkGridFunctionChange(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_CoefficientBilinearForm, GridFunctionChangeOnMixedMesh)
  {
    checkGridFunctionChange(makeMixedMesh());
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
