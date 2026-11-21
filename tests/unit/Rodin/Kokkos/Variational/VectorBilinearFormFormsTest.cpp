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
}
