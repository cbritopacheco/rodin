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

  /**
   * @brief Tags the faces of a mesh with attributes 1, 2 and 3.
   * @param mesh Host mesh, with its faces computed.
   */
  void tagFaces(HostMesh& mesh)
  {
    const size_t d = mesh.getDimension();
    for (Index i = 0; i < mesh.getFaceCount(); ++i)
      mesh.setAttribute({d - 1, i}, Geometry::Attribute(1 + i % 3));
  }

  /**
   * @brief Prepares the faces of a mesh for boundary conditions.
   * @param mesh Host mesh.
   * @return Mesh with its face to cell incidence and tagged faces.
   */
  HostMesh withFaces(HostMesh mesh)
  {
    const size_t d = mesh.getDimension();
    mesh.getConnectivity().compute(d - 1, d);
    tagFaces(mesh);
    return mesh;
  }

  /// @brief Poisson problems with boundary data match the host.
  void checkPoisson(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction hg([](const Geometry::Point& p) { return 2 * p.x() - p.y() + 1; });
    RealFunction dg(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) - Integral(hf, hv) + DirichletBC(hu, hg);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) - Integral(df, dv) + DirichletBC(du, dg);
    expectSameSystem(hp, dp);

    Problem ha(hu, hv);
    ha = Integral(Grad(hu), Grad(hv)) + Integral(hu, hv) - Integral(hf, hv) +
      DirichletBC(hu, hg).on(1) + DirichletBC(hu, RealFunction(3)).on(2);
    Problem da(du, dv);
    da = Integral(Grad(du), Grad(dv)) + Integral(du, dv) - Integral(df, dv) +
      DirichletBC(du, dg).on(1) + DirichletBC(du, RealFunction(3)).on(2);
    expectSameSystem(ha, da);
  }

  /// @brief Boundary data given by a grid function match the host.
  void checkGridFunctionData(const HostMesh& source)
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

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) - Integral(hg, hv) + DirichletBC(hu, hg);
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) - Integral(dg, dv) + DirichletBC(du, dg);
    expectSameSystem(hp, dp);
  }

  /// @brief Vector problems with boundary data match the host.
  void checkVectorProblem(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source, 2);
    P1 dfes(mesh, 2);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    VectorFunction hf{[](const Geometry::Point& p) { return 1 + p.x(); }, Real(-1)};
    VectorFunction df{
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real{return 1 + p.x();
  }
  , Real(-1)
};
VectorFunction hg{[](const Geometry::Point& p) { return p.y(); }, Real(2)};
VectorFunction dg{KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real{return p.y();
}
, Real(2)
}
;

Problem hp(hu, hv);
hp = Integral(Jacobian(hu), Jacobian(hv)) + Integral(Div(hu), Div(hv)) -
  Integral(hf, hv) + DirichletBC(hu, hg).on(FlatSet<Geometry::Attribute>{1, 3});
Problem dp(du, dv);
dp = Integral(Jacobian(du), Jacobian(dv)) + Integral(Div(du), Div(dv)) -
  Integral(df, dv) + DirichletBC(du, dg).on(FlatSet<Geometry::Attribute>{1, 3});
expectSameSystem(hp, dp);
}

/// @brief Mixed velocity and pressure problems match the host.
void checkMixedProblem(const HostMesh& source)
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

  VectorFunction hf{[](const Geometry::Point& p) { return p.x() * p.y(); }, Real(1)};
  VectorFunction df{
    KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real{return p.x() * p.y();
}
, Real(1)
}
;
RealFunction hs([](const Geometry::Point& p) { return p.x() - p.y(); });
RealFunction ds(
  KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() - p.y(); });
const Real eps = 1e-3;

Problem hstokes(hu, hp, hv, hq);
hstokes = Integral(Jacobian(hu), Jacobian(hv)) - Integral(hp, Div(hv)) +
  Integral(Div(hu), hq) + Integral(eps * hp, hq) - Integral(hf, hv) - Integral(hs, hq) +
  DirichletBC(hu, VectorFunction{0, 0}).on(1) + DirichletBC(hp, RealFunction(1)).on(2);
Problem dstokes(du, dp, dv, dq);
dstokes = Integral(Jacobian(du), Jacobian(dv)) - Integral(dp, Div(dv)) +
  Integral(Div(du), dq) + Integral(eps * dp, dq) - Integral(df, dv) - Integral(ds, dq) +
  DirichletBC(du, VectorFunction{0, 0}).on(1) + DirichletBC(dp, RealFunction(1)).on(2);
expectSameSystem(hstokes, dstokes);
}

/// @brief Boundary conditions without face to cell incidence raise.
void checkMissingIncidenceRaises()
{
  auto source = makeMixedMesh();
  const DeviceMesh mesh(TestContext{}, source);
  if (mesh.getConnectivity().hasIncidence(1, 2))
    GTEST_SKIP() << "The mixed mesh already has its face to cell incidence.";
  P1 dfes(mesh);
  TrialFunction du(dfes);
  TestFunction dv(dfes);
  Problem dp(du, dv);
  dp = Integral(Grad(du), Grad(dv)) + DirichletBC(du, RealFunction(0));
  EXPECT_ANY_THROW(dp.assemble());
}

TEST(Kokkos_Variational_Problem, PoissonOnTriangles)
{
  checkPoisson(withFaces(makeTriangleMesh()));
}

TEST(Kokkos_Variational_Problem, PoissonOnMixedMesh)
{
  checkPoisson(withFaces(makeMixedMesh()));
}

TEST(Kokkos_Variational_Problem, PoissonOnTetrahedra)
{
  checkPoisson(withFaces(makeTetrahedronMesh()));
}

TEST(Kokkos_Variational_Problem, GridFunctionDataOnMixedMesh)
{
  checkGridFunctionData(withFaces(makeMixedMesh()));
}

TEST(Kokkos_Variational_Problem, VectorOnTriangles)
{
  checkVectorProblem(withFaces(makeTriangleMesh()));
}

TEST(Kokkos_Variational_Problem, VectorOnMixedMesh)
{
  checkVectorProblem(withFaces(makeMixedMesh()));
}

TEST(Kokkos_Variational_Problem, MixedOnTriangles)
{
  checkMixedProblem(withFaces(makeTriangleMesh()));
}

TEST(Kokkos_Variational_Problem, MixedOnMixedMesh)
{
  checkMixedProblem(withFaces(makeMixedMesh()));
}

TEST(Kokkos_Variational_Problem, MissingIncidenceRaises)
{
  checkMissingIncidenceRaises();
}
}
