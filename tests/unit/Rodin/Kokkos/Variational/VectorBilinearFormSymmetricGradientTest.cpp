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
}
