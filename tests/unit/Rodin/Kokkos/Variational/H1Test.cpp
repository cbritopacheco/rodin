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

  /// @brief Degree of the device H1 spaces.
  constexpr std::integral_constant<size_t, 2> Quadratic{};

  /**
   * @brief Compares a device assembled problem with its host reference.
   * @param host Problem assembled on the local mesh.
   * @param device Problem assembled on the resident mesh.
   */
  template <class HostProblem, class DeviceProblem>
  void expectSameSystem(HostProblem& host, DeviceProblem& device)
  {
    host.assemble();
    device.assemble();
    const auto& h = host.getLinearSystem();
    const auto& d = device.getLinearSystem();
    expectSameOperator(h.getOperator(), d.getOperator());
    expectSameVector(h.getVector(), d.getVector());
  }

  /// @brief Scalar quadratic forms match the host.
  template <size_t K = 2>
  void checkScalarForms(const HostMesh& source)
  {
    constexpr std::integral_constant<size_t, K> degree{};
    const DeviceMesh mesh(TestContext{}, source);
    H1 hfes(degree, source);
    H1 dfes(degree, mesh);
    ASSERT_EQ(hfes.getSize(), dfes.getSize());
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() * p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() * p.y();
    });
    GridFunction hg(hfes);
    fillCoefficients(hg.getData());
    GridFunction dg(dfes);
    dg.getData() = hg.getData();

    expectSameForm(hu, hv, Integral(hu, hv), du, dv, Integral(du, dv));
    expectSameForm(
      hu, hv, Integral(Grad(hu), Grad(hv)), du, dv, Integral(Grad(du), Grad(dv)));
    expectSameForm(hu, hv, Integral(hf * Grad(hu), Grad(hv)), du, dv,
      Integral(df * Grad(du), Grad(dv)));
    expectSameForm(hu, hv, Integral(hg * hu, hv), du, dv, Integral(dg * du, dv));
    expectSameForm(hu, hv, Integral(Dot(Grad(hg), Grad(hu)), hv), du, dv,
      Integral(Dot(Grad(dg), Grad(du)), dv));

    expectSameLinearForm(hv, Integral(hf, hv), dv, Integral(df, dv));
    expectSameLinearForm(hv, Integral(hg, hv), dv, Integral(dg, dv));
    expectSameLinearForm(
      hv, Integral(Dot(Grad(hg), Grad(hv))), dv, Integral(Dot(Grad(dg), Grad(dv))));
  }

  /// @brief Vector quadratic forms match the host.
  template <size_t K = 2>
  void checkVectorForms(const HostMesh& source)
  {
    constexpr std::integral_constant<size_t, K> degree{};
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    H1 hfes(degree, source, d);
    H1 dfes(degree, mesh, d);
    ASSERT_EQ(hfes.getSize(), dfes.getSize());
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });
    GridFunction hg(hfes);
    fillCoefficients(hg.getData());
    GridFunction dg(dfes);
    dg.getData() = hg.getData();

    expectSameForm(hu, hv, Integral(hu, hv), du, dv, Integral(du, dv));
    expectSameForm(hu, hv, Integral(Jacobian(hu), Jacobian(hv)), du, dv,
      Integral(Jacobian(du), Jacobian(dv)));
    expectSameForm(hu, hv, Integral(hf * Jacobian(hu), Jacobian(hv)), du, dv,
      Integral(df * Jacobian(du), Jacobian(dv)));
    expectSameForm(
      hu, hv, Integral(Div(hu), Div(hv)), du, dv, Integral(Div(du), Div(dv)));
    expectSameForm(
      hu, hv, Integral(Jacobian(hu) * hg, hv), du, dv, Integral(Jacobian(du) * dg, dv));

    expectSameLinearForm(hv, Integral(Dot(hg, hv)), dv, Integral(Dot(dg, dv)));
    expectSameLinearForm(
      hv, Integral(Div(hg) * Div(hv)), dv, Integral(Div(dg) * Div(dv)));
    expectSameLinearForm(hv, Integral(Dot(Jacobian(hg), Jacobian(hv))), dv,
      Integral(Dot(Jacobian(dg), Jacobian(dv))));
  }

  /// @brief Linear elasticity on vector quadratic spaces matches the host.
  void checkLinearElasticity(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    H1 hfes(Quadratic, source, d);
    H1 dfes(Quadratic, mesh, d);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hlambda([](const Geometry::Point& p) { return 2 + p.x() * p.y(); });
    RealFunction dlambda(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() * p.y();
    });

    expectSameForm(hu, hv, LinearElasticityIntegral(hu, hv)(Real(1.5), Real(0.7)), du, dv,
      LinearElasticityIntegral(du, dv)(Real(1.5), Real(0.7)));
    expectSameForm(hu, hv, LinearElasticityIntegral(hu, hv)(hlambda, Real(0.7)), du, dv,
      LinearElasticityIntegral(du, dv)(dlambda, Real(0.7)));
  }

  /// @brief Taylor-Hood coupling forms, with a P1 pressure, match the host.
  void checkTaylorHoodForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    H1 hvh(Quadratic, source, d);
    P1 hqh(source);
    H1 dvh(Quadratic, mesh, d);
    P1 dqh(mesh);
    TrialFunction hu(hvh);
    TrialFunction hp(hqh);
    TestFunction hv(hvh);
    TestFunction hq(hqh);
    TrialFunction du(dvh);
    TrialFunction dp(dqh);
    TestFunction dv(dvh);
    TestFunction dq(dqh);

    expectSameForm(hu, hq, Integral(Div(hu), hq), du, dq, Integral(Div(du), dq));
    expectSameForm(hp, hv, Integral(hp, Div(hv)), dp, dv, Integral(dp, Div(dv)));
    expectSameForm(hp, hq, Integral(hp, hq), dp, dq, Integral(dp, dq));
  }

  /// @brief Quadratic face forms match the host.
  template <size_t K = 2>
  void checkFaceForms(const HostMesh& source)
  {
    constexpr std::integral_constant<size_t, K> degree{};
    const DeviceMesh mesh(TestContext{}, source);
    H1 hfes(degree, source);
    H1 dfes(degree, mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });
    GridFunction hg(hfes);
    fillCoefficients(hg.getData());
    GridFunction dg(dfes);
    dg.getData() = hg.getData();

    expectSameForm(hu, hv, BoundaryIntegral(hu, hv), du, dv, BoundaryIntegral(du, dv));
    expectSameForm(
      hu, hv, FaceIntegral(hu, hv).over(2), du, dv, FaceIntegral(du, dv).over(2));
    expectSameLinearForm(hv, BoundaryIntegral(hf, hv), dv, BoundaryIntegral(df, dv));
    expectSameLinearForm(hv, InterfaceIntegral(hg, hv), dv, InterfaceIntegral(dg, dv));
  }

  /// @brief Quadratic problems with boundary data match the host.
  template <size_t K = 2>
  void checkProblems(const HostMesh& source)
  {
    constexpr std::integral_constant<size_t, K> degree{};
    const DeviceMesh mesh(TestContext{}, source);
    H1 hfes(degree, source);
    H1 dfes(degree, mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hg([](const Geometry::Point& p) { return 1 + p.x() * p.x() - p.y(); });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.x() - p.y();
    });
    GridFunction hw(hfes);
    fillCoefficients(hw.getData());
    GridFunction dw(dfes);
    dw.getData() = hw.getData();

    Problem hp(hu, hv);
    const auto hgu = Grad(hu);
    const auto hgv = Grad(hv);
    hp = Integral(hgu, hgv) - Integral(RealFunction(1), hv) + DirichletBC(hu, hg).on(1) +
      DirichletBC(hu, hw).on(2);
    Problem dp(du, dv);
    const auto dgu = Grad(du);
    const auto dgv = Grad(dv);
    dp = Integral(dgu, dgv) - Integral(RealFunction(1), dv) + DirichletBC(du, dg).on(1) +
      DirichletBC(du, dw).on(2);
    expectSameSystem(hp, dp);

    const size_t d = source.getSpaceDimension();
    H1 hvh(degree, source, d);
    H1 dvh(degree, mesh, d);
    TrialFunction hx(hvh);
    TestFunction hy(hvh);
    TrialFunction dx(dvh);
    TestFunction dy(dvh);
    Problem hq(hx, hy);
    const auto hjx = Jacobian(hx);
    const auto hjy = Jacobian(hy);
    hq = Integral(hjx, hjy) - BoundaryIntegral(Dot(BoundaryNormal(source), hy)) +
      DirichletBC(hx, VectorFunction{RealFunction(hg), Real(1), Real(-1)}).on(1);
    Problem dq(dx, dy);
    const auto djx = Jacobian(dx);
    const auto djy = Jacobian(dy);
    dq = Integral(djx, djy) - BoundaryIntegral(Dot(BoundaryNormal(mesh), dy)) +
      DirichletBC(dx, VectorFunction{RealFunction(dg), Real(1), Real(-1)}).on(1);
    if (d == 3)
      expectSameSystem(hq, dq);
  }

  /// @brief Taylor-Hood Stokes with a mean pressure multiplier matches the host.
  void checkStokes(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    H1 hvh(Quadratic, source, 2);
    P1 hqh(source);
    P0g hmh(source);
    H1 dvh(Quadratic, mesh, 2);
    P1 dqh(mesh);
    P0g dmh(mesh);
    TrialFunction hu(hvh);
    TrialFunction hp(hqh);
    TrialFunction hl(hmh);
    TestFunction hv(hvh);
    TestFunction hq(hqh);
    TestFunction hm(hmh);
    TrialFunction du(dvh);
    TrialFunction dp(dqh);
    TrialFunction dl(dmh);
    TestFunction dv(dvh);
    TestFunction dq(dqh);
    TestFunction dm(dmh);

    VectorFunction hf{[](const Geometry::Point& p) { return p.x() * p.y(); }, Real(1)};
    VectorFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() * p.y(); },
      Real(1));

    Problem hs(hu, hp, hl, hv, hq, hm);
    hs = Integral(Jacobian(hu), Jacobian(hv)) - Integral(hp, Div(hv)) +
      Integral(Div(hu), hq) + Integral(hl, hq) + Integral(hp, hm) - Integral(hf, hv) +
      DirichletBC(hu, VectorFunction{0, 0});
    Problem ds(du, dp, dl, dv, dq, dm);
    ds = Integral(Jacobian(du), Jacobian(dv)) - Integral(dp, Div(dv)) +
      Integral(Div(du), dq) + Integral(dl, dq) + Integral(dp, dm) - Integral(df, dv) +
      DirichletBC(du, VectorFunction{0, 0});
    expectSameSystem(hs, ds);
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOnTriangles)
  {
    checkScalarForms(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOnTetrahedra)
  {
    checkScalarForms(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, VectorFormsOnTriangles)
  {
    checkVectorForms(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, VectorFormsOnTetrahedra)
  {
    checkVectorForms(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, LinearElasticityOnTriangles)
  {
    checkLinearElasticity(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, LinearElasticityOnTetrahedra)
  {
    checkLinearElasticity(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, TaylorHoodFormsOnTriangles)
  {
    checkTaylorHoodForms(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, TaylorHoodFormsOnTetrahedra)
  {
    checkTaylorHoodForms(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, FaceFormsOnTriangles)
  {
    checkFaceForms(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, FaceFormsOnTetrahedra)
  {
    checkFaceForms(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, ProblemsOnTriangles)
  {
    checkProblems(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, ProblemsOnTetrahedra)
  {
    checkProblems(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, StokesOnTriangles)
  {
    checkStokes(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, MissingIncidencesAreComputed)
  {
    const HostMesh source = makeTriangleMesh();
    const DeviceMesh mesh(TestContext{}, source);
    H1 dfes(Quadratic, mesh);
    const HostMesh closed = withClosure(makeTriangleMesh());
    H1 hfes(Quadratic, closed);
    EXPECT_EQ(dfes.getSize(), hfes.getSize());
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOfDegreeThreeOnTriangles)
  {
    checkScalarForms<3>(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOfDegreeThreeOnTetrahedra)
  {
    checkScalarForms<3>(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, VectorFormsOfDegreeThreeOnTriangles)
  {
    checkVectorForms<3>(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_H1, ProblemsOfDegreeThreeOnTetrahedra)
  {
    checkProblems<3>(withClosure(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOnQuadrilaterals)
  {
    checkScalarForms(withClosure(makeQuadrilateralMesh()));
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOnHexahedra)
  {
    checkScalarForms(withClosure(makeHexahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOnMixedMesh)
  {
    checkScalarForms(withClosure(makeMixedMesh()));
  }

  TEST(Kokkos_Variational_H1, ScalarFormsOfDegreeThreeOnQuadrilaterals)
  {
    checkScalarForms<3>(withClosure(makeQuadrilateralMesh()));
  }

  TEST(Kokkos_Variational_H1, VectorFormsOnQuadrilaterals)
  {
    checkVectorForms(withClosure(makeQuadrilateralMesh()));
  }

  TEST(Kokkos_Variational_H1, VectorFormsOnHexahedra)
  {
    checkVectorForms(withClosure(makeHexahedronMesh()));
  }

  /// One team per cell holds the local matrices beyond the thread capacity.
  TEST(Kokkos_Variational_H1, ScalarFormsOfDegreeThreeOnHexahedra)
  {
    checkScalarForms<3>(withClosure(makeHexahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, VectorFormsOfDegreeThreeOnHexahedra)
  {
    checkVectorForms<3>(withClosure(makeHexahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, FaceFormsOnMixedMesh)
  {
    checkFaceForms(withClosure(makeMixedMesh()));
  }

  TEST(Kokkos_Variational_H1, FaceFormsOnHexahedra)
  {
    checkFaceForms(withClosure(makeHexahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, ProblemsOnMixedMesh)
  {
    checkProblems(withClosure(makeMixedMesh()));
  }

  TEST(Kokkos_Variational_H1, ProblemsOnHexahedra)
  {
    checkProblems(withClosure(makeHexahedronMesh()));
  }

  TEST(Kokkos_Variational_H1, PyramidMeshRaises)
  {
    const HostMesh source = HostMesh::Builder()
                              .initialize(3)
                              .nodes(5)
                              .vertex({0, 0, 0})
                              .vertex({1, 0, 0})
                              .vertex({1, 1, 0})
                              .vertex({0, 1, 0})
                              .vertex({0, 0, 1})
                              .polytope(Polytope::Type::Pyramid, {0, 1, 2, 3, 4})
                              .finalize();
    const DeviceMesh mesh(TestContext{}, source);
    EXPECT_THROW(H1(Quadratic, mesh), Alert::Exception);
  }
}
