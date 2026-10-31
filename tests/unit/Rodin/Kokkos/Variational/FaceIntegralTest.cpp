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
   * @brief Computes the faces of a mesh and tags them.
   * @param source Mesh whose faces are computed.
   * @return Mesh with face attributes 1, 2 and 3 in turn.
   */
  HostMesh withFaces(HostMesh source)
  {
    const size_t d = source.getDimension();
    source.getConnectivity().compute(d - 1, d);
    source.getConnectivity().compute(d - 1, 0);
    for (Index i = 0; i < source.getFaceCount(); ++i)
      source.setAttribute({d - 1, i}, Geometry::Attribute(1 + i % 3));
    return source;
  }

  /// @brief Scalar face bilinear forms match the host.
  void checkScalarFaceForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });

    expectSameForm(hu, hv, BoundaryIntegral(hu, hv), du, dv, BoundaryIntegral(du, dv));
    expectSameForm(hu, hv, restricted(BoundaryIntegral(hu, hv), 2), du, dv,
      restricted(BoundaryIntegral(du, dv), 2));
    expectSameForm(
      hu, hv, BoundaryIntegral(hf * hu, hv), du, dv, BoundaryIntegral(df * du, dv));
    expectSameForm(hu, hv, FaceIntegral(hu, hv), du, dv, FaceIntegral(du, dv));
    expectSameForm(hu, hv, restricted(FaceIntegral(hu, hv), 1, 3), du, dv,
      restricted(FaceIntegral(du, dv), 1, 3));
    expectSameForm(hu, hv, InterfaceIntegral(hu, hv), du, dv, InterfaceIntegral(du, dv));
  }

  /// @brief Vector face bilinear forms match the host.
  void checkVectorFaceForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    expectSameForm(hu, hv, BoundaryIntegral(hu, hv), du, dv, BoundaryIntegral(du, dv));
    expectSameForm(hu, hv, restricted(FaceIntegral(hu, hv), 2), du, dv,
      restricted(FaceIntegral(du, dv), 2));
  }

  /// @brief Face linear forms, including the boundary normal, match the host.
  void checkFaceLinearForms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source);
    P1 dfes(mesh);
    P1 hvh(source, d);
    P1 dvh(mesh, d);
    TestFunction hv(hfes);
    TestFunction dv(dfes);
    TestFunction hw(hvh);
    TestFunction dw(dvh);

    GridFunction hg(hfes);
    fill(hg, source, 1);
    GridFunction dg(dfes);
    dg.getData() = hg.getData();

    RealFunction hf([](const Geometry::Point& p) { return 2 + p.x() - p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 + p.x() - p.y();
    });

    expectSameLinearForm(hv, BoundaryIntegral(hf, hv), dv, BoundaryIntegral(df, dv));
    expectSameLinearForm(hv, restricted(BoundaryIntegral(hg, hv), 1), dv,
      restricted(BoundaryIntegral(dg, dv), 1));
    expectSameLinearForm(hv, FaceIntegral(hf, hv), dv, FaceIntegral(df, dv));
    expectSameLinearForm(hv, InterfaceIntegral(hg, hv), dv, InterfaceIntegral(dg, dv));
    expectSameLinearForm(hw, BoundaryIntegral(Dot(BoundaryNormal(source), hw)), dw,
      BoundaryIntegral(Dot(BoundaryNormal(mesh), dw)));
  }

  /// @brief Linear elasticity matches the host P1 integrator.
  void checkLinearElasticity(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
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

  TEST(Kokkos_Variational_FaceIntegral, ScalarFormsOnTriangles)
  {
    checkScalarFaceForms(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, ScalarFormsOnMixedMesh)
  {
    checkScalarFaceForms(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, ScalarFormsOnTetrahedra)
  {
    checkScalarFaceForms(withFaces(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, VectorFormsOnTriangles)
  {
    checkVectorFaceForms(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, VectorFormsOnTetrahedra)
  {
    checkVectorFaceForms(withFaces(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, LinearFormsOnTriangles)
  {
    checkFaceLinearForms(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, LinearFormsOnMixedMesh)
  {
    checkFaceLinearForms(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, LinearFormsOnTetrahedra)
  {
    checkFaceLinearForms(withFaces(makeTetrahedronMesh()));
  }

  TEST(Kokkos_Variational_FaceIntegral, MissingFacesRaise)
  {
    const HostMesh source = makeTriangleMesh();
    const DeviceMesh mesh(TestContext{}, source);
    P1 dfes(mesh);
    TrialFunction du(dfes);
    TestFunction dv(dfes);
    BilinearForm form(du, dv);
    form = BoundaryIntegral(du, dv);
    EXPECT_THROW(form.assemble(), Alert::Exception);
  }

  TEST(Kokkos_Variational_LinearElasticity, MatchesHostOnTriangles)
  {
    checkLinearElasticity(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_LinearElasticity, MatchesHostOnMixedMesh)
  {
    checkLinearElasticity(makeMixedMesh());
  }

  TEST(Kokkos_Variational_LinearElasticity, MatchesHostOnTetrahedra)
  {
    checkLinearElasticity(makeTetrahedronMesh());
  }
}
