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

  static_assert(std::is_same_v<P1<Real, DeviceMesh>::MeshType, DeviceMesh>);

  /// @brief Mass and stiffness matrices match the host on a resident mesh.
  void checkMassAndStiffness(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    BilinearForm hm(hu, hv);
    hm = Integral(hu, hv);
    hm.assemble();
    BilinearForm dm(du, dv);
    dm = Integral(du, dv);
    dm.assemble();
    expectSameOperator(hm.getOperator(), dm.getOperator());

    BilinearForm hk(hu, hv);
    hk = Integral(Grad(hu), Grad(hv));
    hk.assemble();
    BilinearForm dk(du, dv);
    dk = Integral(Grad(du), Grad(dv));
    dk.assemble();
    expectSameOperator(hk.getOperator(), dk.getOperator());

    BilinearForm hs(hu, hv);
    hs = Integral(Grad(hu), Grad(hv)) + Integral(Real(3) * hu, hv);
    hs.assemble();
    BilinearForm ds(du, dv);
    ds = Integral(Grad(du), Grad(dv)) + Integral(Real(3) * du, dv);
    ds.assemble();
    expectSameOperator(hs.getOperator(), ds.getOperator());
  }

  /// @brief A device coefficient matches the same host coefficient.
  void checkCoefficient(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() + 2 * p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() + 2 * p.y();
    });

    BilinearForm ha(hu, hv);
    ha = Integral(hf * hu, hv) + Integral(hf * Grad(hu), Grad(hv));
    ha.assemble();
    BilinearForm da(du, dv);
    da = Integral(df * du, dv) + Integral(df * Grad(du), Grad(dv));
    da.assemble();
    expectSameOperator(ha.getOperator(), da.getOperator());
  }

  /// @brief Attribute restriction matches the host.
  void checkAttributes(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    BilinearForm ha(hu, hv);
    auto hak = Integral(Grad(hu), Grad(hv));
    hak.over(1);
    auto ham = Integral(hu, hv);
    ham.over(2);
    ha = hak + ham;
    ha.assemble();
    BilinearForm da(du, dv);
    auto dak = Integral(Grad(du), Grad(dv));
    dak.over(1);
    auto dam = Integral(du, dv);
    dam.over(2);
    da = dak + dam;
    da.assemble();
    expectSameOperator(ha.getOperator(), da.getOperator());
  }

  /// @brief Host evaluation of a device integral raises.
  void checkHostEvaluationRaises(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 dfes(mesh);
    TrialFunction du(dfes);
    TestFunction dv(dfes);
    Integral integral(du, dv);
    EXPECT_ANY_THROW(integral.integrate(0, 0));
  }

  TEST(Kokkos_Variational_BilinearForm, MassAndStiffnessOnTriangles)
  {
    checkMassAndStiffness(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_BilinearForm, MassAndStiffnessOnMixedMesh)
  {
    checkMassAndStiffness(makeMixedMesh());
  }

  TEST(Kokkos_Variational_BilinearForm, MassAndStiffnessOnTetrahedra)
  {
    checkMassAndStiffness(HostMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 4, 3}));
  }

  TEST(Kokkos_Variational_BilinearForm, CoefficientOnTriangles)
  {
    checkCoefficient(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_BilinearForm, CoefficientOnMixedMesh)
  {
    checkCoefficient(makeMixedMesh());
  }

  TEST(Kokkos_Variational_BilinearForm, Attributes)
  {
    checkAttributes(makeMixedMesh());
  }

  TEST(Kokkos_Variational_BilinearForm, HostEvaluationRaises)
  {
    checkHostEvaluationRaises(makeTriangleMesh());
  }
}
