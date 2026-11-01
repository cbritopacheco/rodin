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

  static_assert(std::is_same_v<decltype(P1(std::declval<const DeviceMesh&>(), size_t(2))),
    P1<Math::SpatialVector<Real>, DeviceMesh>>);

  /// @brief Vector mass, Jacobian and divergence forms match the host.
  void checkVectorForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
    ASSERT_EQ(hfes.getSize(), dfes.getSize());
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    expectSameForm(hu, hv, Integral(hu, hv), du, dv, Integral(du, dv));
    expectSameForm(hu, hv, Integral(Jacobian(hu), Jacobian(hv)), du, dv,
      Integral(Jacobian(du), Jacobian(dv)));
    expectSameForm(
      hu, hv, Integral(Div(hu), Div(hv)), du, dv, Integral(Div(du), Div(dv)));
    expectSameForm(
      hu, hv, Integral(-Div(hu), Div(hv)), du, dv, Integral(-Div(du), Div(dv)));
  }

  /// @brief The symmetric gradient form matches the host.
  void checkSymmetricGradient(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    const auto hsu = Real(0.5) * (Jacobian(hu) + Jacobian(hu).T());
    const auto hsv = Real(0.5) * (Jacobian(hv) + Jacobian(hv).T());
    const auto dsu = Real(0.5) * (Jacobian(du) + Jacobian(du).T());
    const auto dsv = Real(0.5) * (Jacobian(dv) + Jacobian(dv).T());
    expectSameForm(hu, hv, Integral(hsu, hsv), du, dv, Integral(dsu, dsv));
  }

  /// @brief Device coefficients of vector forms match the same host coefficients.
  void checkVectorCoefficient(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });
    expectSameForm(hu, hv, Integral(hf * Jacobian(hu), Jacobian(hv)), du, dv,
      Integral(df * Jacobian(du), Jacobian(dv)));
    expectSameForm(hu, hv, Integral(Jacobian(hu) * hf, Jacobian(hv)), du, dv,
      Integral(Jacobian(du) * df, Jacobian(dv)));
    expectSameForm(hu, hv, Integral(hf * hu, hv), du, dv, Integral(df * du, dv));
  }

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
      for (size_t j = 0; j < d; ++j)
        a(i, j) = 1 + Real(i) - 2 * Real(j) + Real(i == j);
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

  TEST(Kokkos_Variational_VectorBilinearForm, FormsOnTriangles)
  {
    checkVectorForms(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, FormsOnMixedMesh)
  {
    checkVectorForms(makeMixedMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, FormsOnTetrahedra)
  {
    checkVectorForms(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, FormsOnQuadrilaterals)
  {
    checkVectorForms(makeQuadrilateralMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, FormsOnHexahedra)
  {
    checkVectorForms(makeHexahedronMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, SymmetricGradientOnTriangles)
  {
    checkSymmetricGradient(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, SymmetricGradientOnTetrahedra)
  {
    checkSymmetricGradient(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, SymmetricGradientOnHexahedra)
  {
    checkSymmetricGradient(makeHexahedronMesh());
  }

  TEST(Kokkos_Variational_VectorBilinearForm, CoefficientOnMixedMesh)
  {
    checkVectorCoefficient(makeMixedMesh());
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
