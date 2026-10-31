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

  /// @brief Scalar linear forms match the host.
  void checkScalarLinearForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source);
    P1 dfes(mesh);
    TestFunction hv(hfes);
    TestFunction dv(dfes);
    P1 hvh(source, d);
    P1 dvh(mesh, d);

    GridFunction hg(hfes);
    fill(hg, source, 1);
    GridFunction dg(dfes);
    dg.getData() = hg.getData();
    GridFunction hw(hvh);
    fill(hw, source, d);
    GridFunction dw(dvh);
    dw.getData() = hw.getData();

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });

    expectSameLinearForm(hv, Integral(hv), dv, Integral(dv));
    expectSameLinearForm(hv, Integral(hf, hv), dv, Integral(df, dv));
    expectSameLinearForm(hv, Integral(hg, hv), dv, Integral(dg, dv));
    expectSameLinearForm(hv, Integral(hf * hg, hv), dv, Integral(df * dg, dv));
    expectSameLinearForm(hv, Integral(hw, Grad(hv)), dv, Integral(dw, Grad(dv)));
  }

  /// @brief Vector linear forms of the Oseen right-hand side match the host.
  void checkVectorLinearForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hvh(source, d);
    P1 hsh(source);
    P1 dvh(mesh, d);
    P1 dsh(mesh);
    TestFunction hv(hvh);
    TestFunction dv(dvh);

    GridFunction hw(hvh);
    fill(hw, source, d);
    GridFunction dw(dvh);
    dw.getData() = hw.getData();
    GridFunction hg(hsh);
    fill(hg, source, 1);
    GridFunction dg(dsh);
    dg.getData() = hg.getData();

    expectSameLinearForm(hv, Integral(hw, hv), dv, Integral(dw, dv));
    expectSameLinearForm(
      hv, Integral(Jacobian(hw) * hw, hv), dv, Integral(Jacobian(dw) * dw, dv));
    expectSameLinearForm(hv, Integral(hg, Div(hv)), dv, Integral(dg, Div(dv)));
    expectSameLinearForm(
      hv, Integral(Jacobian(hw), Jacobian(hv)), dv, Integral(Jacobian(dw), Jacobian(dv)));
  }

  /// @brief Component-wise vector function sources match the host.
  void checkVectorFunctionSource(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source, 2);
    P1 dfes(mesh, 2);
    TestFunction hv(hfes);
    TestFunction dv(dfes);

    VectorFunction hf{
      [](const Geometry::Point& p) { return 1 + p.x() * p.y(); }, Real(0.5)};
    VectorFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
        return 1 + p.x() * p.y();
      },
      Real(0.5));
    expectSameLinearForm(hv, Integral(hf, hv), dv, Integral(df, dv));
  }

/// @brief Subtracted and restricted integrators match the host.
  void checkCombinedLinearForm(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TestFunction hv(hfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.x(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.x();
    });

    LinearForm hl(hv);
    hl = Integral(hv);
    hl -= Integral(hf, hv).over(2);
    hl.assemble();
    LinearForm dl(dv);
    dl = Integral(dv);
    dl -= Integral(df, dv).over(2);
    dl.assemble();
    expectSameVector(hl.getVector(), dl.getVector());
  }

  TEST(Kokkos_Variational_LinearForm, ScalarOnTriangles)
  {
    checkScalarLinearForms(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_LinearForm, ScalarOnMixedMesh)
  {
    checkScalarLinearForms(makeMixedMesh());
  }

  TEST(Kokkos_Variational_LinearForm, ScalarOnTetrahedra)
  {
    checkScalarLinearForms(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_LinearForm, VectorOnTriangles)
  {
    checkVectorLinearForms(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_LinearForm, VectorOnMixedMesh)
  {
    checkVectorLinearForms(makeMixedMesh());
  }

  TEST(Kokkos_Variational_LinearForm, VectorOnTetrahedra)
  {
    checkVectorLinearForms(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_LinearForm, VectorFunctionOnTriangles)
  {
    checkVectorFunctionSource(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_LinearForm, CombinedOnMixedMesh)
  {
    checkCombinedLinearForm(makeMixedMesh());
  }
}
