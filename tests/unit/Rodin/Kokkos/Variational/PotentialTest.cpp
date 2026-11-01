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
   * @brief Assembles a host and a device form and compares their operators.
   * @param hu Host trial function.
   * @param hv Host test function.
   * @param host Assigns the host integrals to a form.
   * @param du Device trial function.
   * @param dv Device test function.
   * @param device Assigns the same integrals to a device form.
   */
  template <class HostTrial, class HostTest, class Host, class DeviceTrial,
    class DeviceTest, class Device>
  void expectSameGlobalForm(const HostTrial& hu, const HostTest& hv, const Host& host,
    const DeviceTrial& du, const DeviceTest& dv, const Device& device)
  {
    BilinearForm hf(hu, hv);
    host(hf);
    hf.assemble();
    BilinearForm df(du, dv);
    device(df);
    df.assemble();
    expectSameOperator(hf.getOperator(), df.getOperator());
  }

  /**
   * @brief Scalar potentials match the host.
   * @param source Host mesh.
   * @param singular Also checks a kernel singular on the diagonal, which only
   * the coincident triangles integrate.
   */
  void checkScalarPotential(const HostMesh& source, bool singular)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    const auto hk = [](const Geometry::Point& x, const Geometry::Point& y) {
      return 1 + x.x() * y.y() + Real(0.5) * (x.y() - y.x());
    };
    const auto dk = KOKKOS_LAMBDA(
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
                      ->Real
    {
      return 1 + x.x() * y.y() + Real(0.5) * (x.y() - y.x());
    };
    const auto hs = [](const Geometry::Point& x, const Geometry::Point& y) {
      const Real dx = x.x() - y.x();
      const Real dy = x.y() - y.y();
      return 1 / std::sqrt(dx * dx + dy * dy);
    };
    const auto ds = KOKKOS_LAMBDA(
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
                      ->Real
    {
      const Real dx = x.x() - y.x();
      const Real dy = x.y() - y.y();
      return 1 / ::Kokkos::sqrt(dx * dx + dy * dy);
    };

    expectSameGlobalForm(
      hu, hv, [&](auto& f) { f = Integral(Potential(hk, hu), hv); }, du, dv,
      [&](auto& f) { f = Integral(Potential(dk, du), dv); });
    expectSameGlobalForm(
      hu, hv,
      [&](auto& f) {
        f = Integral(hu, hv);
        f -= Integral(Potential(hk, hu), hv);
      },
      du, dv,
      [&](auto& f) {
        f = Integral(du, dv);
        f -= Integral(Potential(dk, du), dv);
      });
    if (singular)
    {
      expectSameGlobalForm(
        hu, hv, [&](auto& f) { f += Integral(Potential(hs, hu), hv); }, du, dv,
        [&](auto& f) { f += Integral(Potential(ds, du), dv); });
    }
  }

  /**
   * @brief Potentials restricted to tagged cells match the host.
   * @param source Host mesh with cell attributes 1 and 2.
   */
  void checkRestrictedPotential(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    const auto hk = [](const Geometry::Point& x, const Geometry::Point& y) {
      return 2 + x.x() * y.x() - y.y();
    };
    const auto dk = KOKKOS_LAMBDA(
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
                      ->Real
    {
      return 2 + x.x() * y.x() - y.y();
    };
    auto hi = Integral(Potential(hk, hu), hv);
    auto di = Integral(Potential(dk, du), dv);
    hi.setTrialAttributes({1});
    di.setTrialAttributes({1});
    expectSameGlobalForm(
      hu, hv, [&](auto& f) { f += hi; }, du, dv, [&](auto& f) { f += di; });
    hi.setTestAttributes({2});
    di.setTestAttributes({2});
    expectSameGlobalForm(
      hu, hv,
      [&](auto& f) {
        f = Integral(hu, hv);
        f += hi;
      },
      du, dv,
      [&](auto& f) {
        f = Integral(du, dv);
        f += di;
      });
  }

  /**
   * @brief Vector potentials with a matrix kernel match the host.
   * @param source Host mesh.
   */
  void checkVectorPotential(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    const size_t d = source.getSpaceDimension();
    P1 hfes(source, d);
    P1 dfes(mesh, d);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    const auto hk = [d](Math::SpatialMatrix<Real>& m, const Geometry::Point& x,
                      const Geometry::Point& y) {
      m.resize(d, d);
      for (Index i = 0; i < d; ++i)
      {
        for (Index j = 0; j < d; ++j)
          m(i, j) = (i == j ? 2 : 0) + (i + 1) * x.x() * y.y() - Real(0.25) * j * y.x();
      }
    };
    const auto dk = KOKKOS_LAMBDA(DeviceGeometry::Matrix & m,
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
    {
      for (Index i = 0; i < m.rows(); ++i)
      {
        for (Index j = 0; j < m.cols(); ++j)
          m(i, j) = (i == j ? 2 : 0) + (i + 1) * x.x() * y.y() - Real(0.25) * j * y.x();
      }
    };
    expectSameGlobalForm(
      hu, hv, [&](auto& f) { f = Integral(Potential(hk, hu), hv); }, du, dv,
      [&](auto& f) { f = Integral(Potential(dk, du), dv); });
    expectSameGlobalForm(
      hu, hv,
      [&](auto& f) {
        f = Integral(Jacobian(hu), Jacobian(hv));
        f -= Integral(Potential(hk, hu), hv);
      },
      du, dv,
      [&](auto& f) {
        f = Integral(Jacobian(du), Jacobian(dv));
        f -= Integral(Potential(dk, du), dv);
      });
  }

  /**
   * @brief Problems with a potential term and boundary data match the host.
   * @param source Host mesh with its face to cell incidence.
   */
  void checkPotentialProblem(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 dfes(mesh);
    TrialFunction hu(hfes);
    TestFunction hv(hfes);
    TrialFunction du(dfes);
    TestFunction dv(dfes);

    const auto hk = [](const Geometry::Point& x, const Geometry::Point& y) {
      return 1 + x.x() * y.y();
    };
    const auto dk = KOKKOS_LAMBDA(
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
                      ->Real
    {
      return 1 + x.x() * y.y();
    };
    RealFunction hf([](const Geometry::Point& p) { return 1 + p.x() * p.y(); });
    RealFunction df(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });

    Problem hp(hu, hv);
    hp = Integral(Grad(hu), Grad(hv)) + Integral(Potential(hk, hu), hv) -
      Integral(hf, hv) + DirichletBC(hu, RealFunction(1));
    Problem dp(du, dv);
    dp = Integral(Grad(du), Grad(dv)) + Integral(Potential(dk, du), dv) -
      Integral(df, dv) + DirichletBC(du, RealFunction(1));
    hp.assemble();
    dp.assemble();
    const auto& h = hp.getLinearSystem();
    const auto& s = dp.getLinearSystem();
    expectSameOperator(h.getOperator(), s.getOperator());
    expectSameVector(h.getVector(), s.getVector());
  }

  TEST(Kokkos_Variational_Potential, ScalarOnTriangles)
  {
    checkScalarPotential(makeTriangleMesh(), true);
  }

  TEST(Kokkos_Variational_Potential, ScalarOnTetrahedra)
  {
    checkScalarPotential(makeTetrahedronMesh(), false);
  }

  TEST(Kokkos_Variational_Potential, ScalarOnMixedMesh)
  {
    checkScalarPotential(makeMixedMesh(), false);
    checkRestrictedPotential(makeMixedMesh());
  }

  TEST(Kokkos_Variational_Potential, VectorOnTriangles)
  {
    checkVectorPotential(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_Potential, VectorOnTetrahedra)
  {
    checkVectorPotential(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_Potential, VectorOnMixedMesh)
  {
    checkVectorPotential(makeMixedMesh());
  }

  TEST(Kokkos_Variational_Potential, Problem)
  {
    checkPotentialProblem(withClosure(makeTriangleMesh()));
  }
}
