/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Kokkos.h>

#include "../TestContext.h"

namespace Rodin::Tests::Unit
{
  namespace DeviceGeometry = Rodin::Kokkos::Geometry;
  using HostMesh = Geometry::Mesh<Context::Local>;
  using DeviceMesh = Geometry::Mesh<TestContext>;
  using Polytope = Geometry::Polytope;
  using namespace Rodin::Variational;

  static_assert(std::is_same_v<P1<Real, DeviceMesh>::MeshType, DeviceMesh>);

  /**
   * @brief Compares a device assembled operator with its host reference.
   * @param host Operator assembled on the local mesh.
   * @param device Operator assembled on the resident mesh.
   */
  void expectSameOperator(
    const Math::SparseMatrix<Real>& host, const Math::SparseMatrix<Real>& device)
  {
    ASSERT_EQ(host.rows(), device.rows());
    ASSERT_EQ(host.cols(), device.cols());
    ASSERT_GT(host.norm(), Real(0));
    const Math::SparseMatrix<Real> difference = host - device;
    EXPECT_LE(difference.norm(), Real(1e-12) * host.norm());
  }

  /**
   * @brief Two quadrilaterals and two triangles with tagged cells.
   * @return Mixed mesh with cell attributes 1 and 2.
   */
  HostMesh makeMixedMesh()
  {
    auto mesh = HostMesh::Builder()
                  .initialize(2)
                  .nodes(8)
                  .vertex({0, 0})
                  .vertex({1, 0})
                  .vertex({0, 1})
                  .vertex({1.2, 1.1})
                  .vertex({2, 0})
                  .vertex({2.1, 1})
                  .vertex({0.4, 2})
                  .vertex({1.5, 2.2})
                  .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2})
                  .polytope(Polytope::Type::Triangle, {1, 4, 3})
                  .polytope(Polytope::Type::Triangle, {4, 5, 3})
                  .polytope(Polytope::Type::Quadrilateral, {2, 3, 7, 6})
                  .finalize();
    mesh.getConnectivity().compute(2, 0);
    for (Index i = 0; i < mesh.getCellCount(); ++i)
      mesh.setAttribute({2, i}, Geometry::Attribute(1 + i % 2));
    return mesh;
  }

  /**
   * @brief Distorted triangle grid.
   * @return Uniform triangle grid with displaced interior vertices.
   */
  HostMesh makeTriangleMesh()
  {
    auto mesh = HostMesh::UniformGrid(Polytope::Type::Triangle, {6, 5});
    mesh.scale(Real(1) / 5);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
    {
      const auto x = mesh.getVertexCoordinates(i);
      Math::SpatialPoint y(2);
      y(0) = x(0) + Real(0.03) * std::sin(Real(3) * x(1));
      y(1) = x(1) + Real(0.02) * std::cos(Real(2) * x(0));
      mesh.setVertexCoordinates(i, y);
    }
    return mesh;
  }

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
