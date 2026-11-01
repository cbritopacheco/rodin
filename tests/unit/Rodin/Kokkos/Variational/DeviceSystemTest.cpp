/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Configure.h>
#include <Rodin/Kokkos/LinearSystem.h>

#ifdef RODIN_USE_CUDSS
#include <Rodin/Solver/CuDSS.h>
#endif

#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /// @brief System whose matrix stays in the memory of the test backend.
  using DeviceLinearSystem =
    Math::LinearSystem<Rodin::Kokkos::SparseMatrix<Real, TestBackend>,
      Math::Vector<Real>>;

  static_assert(
    std::is_same_v<Rodin::Kokkos::SparseMatrix<Real, TestBackend>::MemorySpace,
      TestMemorySpace>);

  /**
   * @brief Checks that the device matrix equals the host copy of the same
   * device assembly, coefficient by coefficient.
   * @param host Problem assembled into a host matrix.
   * @param device Same problem assembled into a device matrix.
   */
  template <class HostProblem, class DeviceProblem>
  void expectSameMatrix(HostProblem& host, DeviceProblem& device)
  {
    host.assemble();
    device.assemble();
    const auto& h = host.getLinearSystem();
    const auto& d = device.getLinearSystem();
    const Math::SparseMatrix<Real> a = d.getOperator().toHost();
    ASSERT_EQ(a.rows(), h.getOperator().rows());
    ASSERT_EQ(a.cols(), h.getOperator().cols());
    ASSERT_EQ(a.nonZeros(), h.getOperator().nonZeros());
    ASSERT_GT(a.nonZeros(), 0);
    // The coefficients are summed in the same order, so they are equal.
    const Math::SparseMatrix<Real> difference = a - h.getOperator();
    EXPECT_EQ(difference.norm(), Real(0));
    ASSERT_EQ(d.getVector().size(), h.getVector().size());
    EXPECT_EQ((d.getVector() - h.getVector()).norm(), Real(0));

    // The rows are sorted, and a second assembly keeps the arrays.
    const auto& m = d.getOperator();
    const auto starts =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), m.getRowStarts());
    const auto columns =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace(), m.getColumns());
    ASSERT_EQ(starts.extent(0), m.getRowCount() + 1);
    EXPECT_EQ(static_cast<size_t>(starts(m.getRowCount())), m.getNonZeroCount());
    for (size_t i = 0; i < m.getRowCount(); ++i)
    {
      for (auto k = starts(i) + 1; k < starts(i + 1); ++k)
        EXPECT_LT(columns(k - 1), columns(k));
    }
    const Real* values = m.getValues().data();
    device.assemble();
    EXPECT_EQ(d.getOperator().getValues().data(), values);
  }

  /// @brief Dirichlet conditions on two attributes.
  void checkPoisson()
  {
    auto source = makeMixedMesh();
    source.getConnectivity().compute(1, 2);
    for (Index i = 0; i < source.getFaceCount(); ++i)
      source.setAttribute({1, i}, Geometry::Attribute(1 + i % 3));
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction g(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 2 * p.x() + 1; });
    Problem host(u, v);
    host = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, g).on(1) + DirichletBC(u, RealFunction(3)).on(2);
    Problem<DeviceLinearSystem, decltype(u), decltype(v)> device(u, v);
    device = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, g).on(1) + DirichletBC(u, RealFunction(3)).on(2);
    expectSameMatrix(host, device);
  }

  /// @brief Periodic and Dirichlet conditions.
  void checkPeriodic()
  {
    const auto source = withTaggedFace(makeTriangleMesh());
    const DeviceMesh mesh(TestContext{}, source);
    P1 hfes(source);
    P1 fes(mesh);
    TrialFunction hu(hfes);
    TrialFunction u(fes);
    TestFunction v(fes);
    const auto pairs = makePairs(fes.getSize(), fixedDOFs(hu, 5));
    RealFunction g(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 2 - p.x(); });
    Problem host(u, v);
    host = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1), v) +
      DirichletBC(u, g).on(5) + PeriodicBC(u, pairs);
    Problem<DeviceLinearSystem, decltype(u), decltype(v)> device(u, v);
    device = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1), v) +
      DirichletBC(u, g).on(5) + PeriodicBC(u, pairs);
    expectSameMatrix(host, device);
  }

  /// @brief Mixed velocity and pressure problem.
  void checkStokes()
  {
    auto source = makeTriangleMesh();
    source.getConnectivity().compute(1, 2);
    for (Index i = 0; i < source.getFaceCount(); ++i)
      source.setAttribute({1, i}, Geometry::Attribute(1 + i % 3));
    const DeviceMesh mesh(TestContext{}, source);
    P1 vh(mesh, 2);
    P1 qh(mesh);
    TrialFunction u(vh);
    TrialFunction p(qh);
    TestFunction v(vh);
    TestFunction q(qh);
    const Real eps = 1e-3;
    Problem host(u, p, v, q);
    host = Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
      Integral(Div(u), q) + Integral(eps * p, q) - Integral(VectorFunction{1, 0}, v) +
      DirichletBC(u, VectorFunction{0, 0}).on(1) + DirichletBC(p, RealFunction(1)).on(2);
    Problem<DeviceLinearSystem, decltype(u), decltype(p), decltype(v), decltype(q)>
      device(u, p, v, q);
    device = Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
      Integral(Div(u), q) + Integral(eps * p, q) - Integral(VectorFunction{1, 0}, v) +
      DirichletBC(u, VectorFunction{0, 0}).on(1) + DirichletBC(p, RealFunction(1)).on(2);
    expectSameMatrix(host, device);
  }

#ifdef RODIN_USE_CUDSS
  /// @brief Solution by CuDSS from the device matrix equals the host one.
  void checkCuDSS()
  {
    auto source = makeTriangleMesh();
    source.getConnectivity().compute(1, 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    Problem host(u, v);
    host = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, RealFunction(1));
    Problem<DeviceLinearSystem, decltype(u), decltype(v)> device(u, v);
    device = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, RealFunction(1));
    Solver::CuDSS hs(host);
    hs.solve();
    ASSERT_TRUE(hs.success());
    const Math::Vector<Real> reference = u.getSolution().getData();
    Solver::CuDSS ds(device);
    static_assert(std::is_same_v<decltype(ds), Solver::CuDSS<DeviceLinearSystem>>);
    ds.solve();
    ASSERT_TRUE(ds.success());
    const Math::Vector<Real> solution = u.getSolution().getData();
    ASSERT_GT(reference.norm(), Real(0));
    EXPECT_LE((solution - reference).norm(), Real(1e-12) * reference.norm());

    // A refactorization after a reassembly keeps the analysis.
    device.assemble();
    ds.factorize(device.getLinearSystem());
    ASSERT_TRUE(ds.success());
    ds.solve(device.getLinearSystem());
    ASSERT_TRUE(ds.success());
    const auto& x = device.getLinearSystem().getSolution();
    EXPECT_LE((x - reference).norm(), Real(1e-12) * reference.norm());
  }
#endif
  TEST(Kokkos_Variational_DeviceSystem, Poisson)
  {
    checkPoisson();
  }

  TEST(Kokkos_Variational_DeviceSystem, Periodic)
  {
    checkPeriodic();
  }

  TEST(Kokkos_Variational_DeviceSystem, Stokes)
  {
    checkStokes();
  }

#ifdef RODIN_USE_CUDSS
  TEST(Kokkos_Variational_DeviceSystem, CuDSS)
  {
    checkCuDSS();
  }
#endif
}
