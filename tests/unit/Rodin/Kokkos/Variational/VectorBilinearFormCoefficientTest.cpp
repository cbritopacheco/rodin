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

  TEST(Kokkos_Variational_VectorBilinearForm, CoefficientOnMixedMesh)
  {
    checkVectorCoefficient(makeMixedMesh());
  }
}
