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

  /**
   * @brief Compares a device assembled linear system with its host reference.
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

  /// @brief Periodic scalar problems match the host.
  /// @param source Mesh tagged by withTaggedFace().
  void checkPeriodic(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);
    const auto pairs = makePairs(hfes.getSize(), fixedDOFs(hu, 5));

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction hg([](const Geometry::Point& p) { return 2 - p.x(); });
    RealFunction dg(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 2 - p.x(); });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) + Integral(hu, hv) - Integral(hf, hv) +
      DirichletBC(hu, hg).on(5) + PeriodicBC(hu, pairs);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      DirichletBC(du, dg).on(5) + PeriodicBC(du, pairs);
    expectSameSystem(hp, dp);
  }

  /// @brief Self identifications, which fix or couple the components, match the host.
  void checkSelfIdentification(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return p.x() - p.y(); });
    RealFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() - p.y(); });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) - Integral(hf, hv) + DirichletBC(hu, -hu).on(1) +
      DirichletBC(hu, RealFunction(3)).on(2);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) - Integral(df, dv) + DirichletBC(du, -du).on(1) +
      DirichletBC(du, RealFunction(3)).on(2);
    expectSameSystem(hp, dp);

    P1 hvh(source, 2);
    P1 dvh(mesh, 2);
    TrialFunction hw(hvh);
    TestFunction hz(hvh);
    TrialFunction dw(dvh);
    TestFunction dz(dvh);
    Math::SpatialMatrix<Real> t(2, 2);
    t.setZero();
    t(0, 0) = Real(0.5);
    t(0, 1) = Real(0.25);
    t(1, 1) = Real(0.5);
    Problem hq(hw, hz);
    hq = Integral(Jacobian(hw), Jacobian(hz)) - Integral(VectorFunction{hf, 1}, hz) +
      DirichletBC(hw, MatrixFunction(t) * hw, VectorFunction{1, 0}).on(1);
    Problem dq(dw, dz);
    dq = Integral(Jacobian(dw), Jacobian(dz)) - Integral(VectorFunction{df, 1}, dz) +
      DirichletBC(dw, MatrixFunction(t) * dw, VectorFunction{1, 0}).on(1);
    expectSameSystem(hq, dq);
  }

  /// @brief Identifications and periodic conditions between fields match the host.
  /// @param source Mesh tagged by withTaggedFace().
  void checkMixed(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hvh(source, 2);
    P1 hqh(source);
    P1 dvh(mesh, 2);
    P1 dqh(mesh);
    TrialFunction hu(hvh);
    TrialFunction hp(hqh);
    TestFunction hv(hvh);
    TestFunction hq(hqh);
    TrialFunction du(dvh);
    TrialFunction dp(dqh);
    TestFunction dv(dvh);
    TestFunction dq(dqh);
    const auto pairs = makePairs(hqh.getSize(), fixedDOFs(hp, 5));

    VectorFunction hf{[](const Geometry::Point& p) { return p.x() * p.y(); }, Real(1)};
    VectorFunction df(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() * p.y(); },
      Real(1));
    RealFunction hs([](const Geometry::Point& p) { return p.x() - p.y(); });
    RealFunction ds(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() - p.y(); });
    const Real eps = 1e-3;

    Problem hstokes(hu, hp, hv, hq);
    hstokes = Integral(Jacobian(hu), Jacobian(hv)) - Integral(hp, Div(hv)) +
      Integral(Div(hu), hq) + Integral(eps * hp, hq) - Integral(hf, hv) -
      Integral(hs, hq) + DirichletBC(hu, VectorFunction{0, 0}).on(1) +
      DirichletBC(hp, RealFunction(1)).on(5) + PeriodicBC(hp, pairs);
    Problem dstokes(du, dp, dv, dq);
    dstokes = Integral(Jacobian(du), Jacobian(dv)) - Integral(dp, Div(dv)) +
      Integral(Div(du), dq) + Integral(eps * dp, dq) - Integral(df, dv) -
      Integral(ds, dq) + DirichletBC(du, VectorFunction{0, 0}).on(1) +
      DirichletBC(dp, RealFunction(1)).on(5) + PeriodicBC(dp, pairs);
    expectSameSystem(hstokes, dstokes);

    P1 hah(source);
    P1 dah(mesh);
    TrialFunction ha(hah);
    TestFunction hb(hah);
    TrialFunction da(dah);
    TestFunction db(dah);
    Problem hcoupled(ha, hp, hb, hq);
    hcoupled = Integral(Grad(ha), Grad(hb)) + Integral(hp, hb) +
      Integral(Grad(hp), Grad(hq)) - Integral(hs, hb) - Integral(hs, hq) +
      DirichletBC(ha, Real(0.5) * hp, RealFunction(2)).on(1) +
      DirichletBC(hp, RealFunction(1)).on(2);
    Problem dcoupled(da, dp, db, dq);
    dcoupled = Integral(Grad(da), Grad(db)) + Integral(dp, db) +
      Integral(Grad(dp), Grad(dq)) - Integral(ds, db) - Integral(ds, dq) +
      DirichletBC(da, Real(0.5) * dp, RealFunction(2)).on(1) +
      DirichletBC(dp, RealFunction(1)).on(2);
    expectSameSystem(hcoupled, dcoupled);
  }

  TEST(Kokkos_Variational_Constraint, PeriodicOnTriangles)
  {
    checkPeriodic(withTaggedFace(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_Constraint, PeriodicOnTetrahedra)
  {
    checkPeriodic(withTaggedFace(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_Constraint, SelfIdentificationOnTriangles)
  {
    checkSelfIdentification(withClosure(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_Constraint, SelfIdentificationOnMixedMesh)
  {
    checkSelfIdentification(withClosure(makeMixedMesh()));
  }

  TEST(Kokkos_Variational_Constraint, MixedOnTriangles)
  {
    checkMixed(withTaggedFace(makeTriangleMesh()));
  }
}
