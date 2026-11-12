/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Configure.h>
#include <Rodin/Kokkos/Math/LinearSystem.h>

#ifdef RODIN_USE_CUDSS
#include <Rodin/Solver/CuDSS.h>
#include <Rodin/Solver/SparseLU.h>
#endif

#ifdef RODIN_USE_ROCSOLVER
#include <Rodin/Solver/RocSolverRF.h>
#include <Rodin/Solver/SparseLU.h>
#endif

#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /// @brief System whose matrix stays in the memory of the test backend.
  using DeviceLinearSystem =
    Math::LinearSystem<Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>,
      Rodin::Kokkos::Math::Vector<Real, TestBackend>>;

  static_assert(
    std::is_same_v<Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>::MemorySpace,
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
    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(u), decltype(v)> host(u, v);
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
    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(u), decltype(v)> host(u, v);
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
    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(u), decltype(p), decltype(v), decltype(q)> host(
      u, p, v, q);
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
    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(u), decltype(v)> host(u, v);
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

    // Exercise cuDSS invalidation against an independent host factorization.
    // These CSR graphs are all nonsingular, including the in-place edits.
    using Matrix = typename DeviceLinearSystem::MatrixType;
    const auto execution = mesh.getContext().getExecutionSpace();
    const auto offsets = Rodin::Kokkos::Internal::upload(
      execution, std::vector<Math::SparseIndex>{0, 2, 3, 5}, "cuDSS test offsets");
    const auto columns = Rodin::Kokkos::Internal::upload(
      execution, std::vector<Math::SparseIndex>{0, 1, 2, 0, 2}, "cuDSS test columns");
    const auto values = Rodin::Kokkos::Internal::upload(
      execution, std::vector<Real>{4, 1, 2, 1, 5}, "cuDSS test values");
    auto& system = device.getLinearSystem();
    system.getOperator() = Matrix(3, 3, offsets, columns, values);
    Math::Vector<Real> rhs(3);
    rhs << 1, 2, 3;
    system.getVector() = rhs;
    Solver::SparseLU referenceSolver(host);
    const auto check = [&]() {
      auto& referenceSystem = host.getLinearSystem();
      referenceSystem.getOperator() = system.getOperator().toHost();
      referenceSystem.getVector() = rhs;
      referenceSolver.solve(referenceSystem);
      ASSERT_TRUE(referenceSolver.success());
      ds.factorize(system);
      ASSERT_TRUE(ds.success());
      ds.solve(system);
      ASSERT_TRUE(ds.success());
      EXPECT_LE((system.getSolution().toHost() - referenceSystem.getSolution()).norm(),
        Real(1e-12) * referenceSystem.getSolution().norm());
      EXPECT_TRUE(ds.getResources().pattern.matches(system.getOperator()));
    };
    check();
    const auto* snapshot = ds.getResources().data;
    system.getOperator() *= Real(2);
    check();
    EXPECT_EQ(ds.getResources().data, snapshot);
    Matrix relocated(system.getOperator());
    system.getOperator() = std::move(relocated);
    check();
    EXPECT_EQ(ds.getResources().data, snapshot);
    const auto currentColumns = system.getOperator().getColumns();
    ::Kokkos::parallel_for(
      "cuDSS edit column", ::Kokkos::RangePolicy<typename Matrix::ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { currentColumns(3) = 1; });
    EXPECT_FALSE(ds.getResources().pattern.matches(system.getOperator()));
    check();
    const auto currentOffsets = system.getOperator().getRowStarts();
    ::Kokkos::parallel_for(
      "cuDSS edit offset", ::Kokkos::RangePolicy<typename Matrix::ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { currentOffsets(1) = 1; });
    EXPECT_FALSE(ds.getResources().pattern.matches(system.getOperator()));
    check();
    ds.clear(Solver::Factorization::Symbolic);
    EXPECT_FALSE(ds.getResources().pattern.matches(system.getOperator()));
    check();
  }
#endif

#ifdef RODIN_USE_ROCSOLVER
  /// @brief Solution by RocSolverRF from the device matrix equals the host one.
  void checkRocSolverRF()
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
    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(u), decltype(v)> host(u, v);
    host = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, RealFunction(1));
    Problem<DeviceLinearSystem, decltype(u), decltype(v)> device(u, v);
    device = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, RealFunction(1));
    Solver::RocSolverRF hs(host);
    hs.solve();
    ASSERT_TRUE(hs.success());
    const Math::Vector<Real> reference = u.getSolution().getData();
    Solver::RocSolverRF ds(device);
    static_assert(std::is_same_v<decltype(ds), Solver::RocSolverRF<DeviceLinearSystem>>);
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

    using Matrix = typename DeviceLinearSystem::MatrixType;
    const auto execution = mesh.getContext().getExecutionSpace();
    const auto offsets = Rodin::Kokkos::Internal::upload(
      execution, std::vector<Math::SparseIndex>{0, 2, 3, 5}, "rocSOLVER test offsets");
    const auto columns = Rodin::Kokkos::Internal::upload(
      execution, std::vector<Math::SparseIndex>{0, 1, 2, 0, 2}, "rocSOLVER test columns");
    const auto values = Rodin::Kokkos::Internal::upload(
      execution, std::vector<Real>{4, 1, 2, 1, 5}, "rocSOLVER test values");
    auto& system = device.getLinearSystem();
    system.getOperator() = Matrix(3, 3, offsets, columns, values);
    Math::Vector<Real> rhs(3);
    rhs << 1, 2, 3;
    system.getVector() = rhs;
    Solver::SparseLU referenceSolver(host);
    const auto check = [&]() {
      auto& referenceSystem = host.getLinearSystem();
      referenceSystem.getOperator() = system.getOperator().toHost();
      referenceSystem.getVector() = rhs;
      referenceSolver.solve(referenceSystem);
      ASSERT_TRUE(referenceSolver.success());
      ds.factorize(system);
      ASSERT_TRUE(ds.success());
      ds.solve(system);
      ASSERT_TRUE(ds.success());
      EXPECT_LE((system.getSolution().toHost() - referenceSystem.getSolution()).norm(),
        Real(1e-12) * referenceSystem.getSolution().norm());
      EXPECT_TRUE(ds.getResources().pattern.matches(system.getOperator()));
    };
    check();
    const auto metadata = ds.getResources().rfinfo;
    system.getOperator() *= Real(2);
    check();
    EXPECT_EQ(ds.getResources().rfinfo, metadata);
    Matrix relocated(system.getOperator());
    system.getOperator() = std::move(relocated);
    check();
    EXPECT_EQ(ds.getResources().rfinfo, metadata);
    const auto currentColumns = system.getOperator().getColumns();
    ::Kokkos::parallel_for(
      "rocSOLVER edit column",
      ::Kokkos::RangePolicy<typename Matrix::ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { currentColumns(3) = 1; });
    EXPECT_FALSE(ds.getResources().pattern.matches(system.getOperator()));
    check();
    const auto currentOffsets = system.getOperator().getRowStarts();
    ::Kokkos::parallel_for(
      "rocSOLVER edit offset",
      ::Kokkos::RangePolicy<typename Matrix::ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { currentOffsets(1) = 1; });
    EXPECT_FALSE(ds.getResources().pattern.matches(system.getOperator()));
    check();

    // Stored zeros keep the graph but must invalidate numeric factors.
    const auto previous = system.getSolution().toHost();
    ::Kokkos::deep_copy(system.getOperator().getValues(), Real(0));
    ds.factorize(system);
    EXPECT_FALSE(ds.success());
    EXPECT_EQ(ds.getInfo().status, decltype(ds)::Breakdown);
    EXPECT_TRUE(ds.getResources().symbolic);
    EXPECT_FALSE(ds.getResources().numeric);
    EXPECT_TRUE(ds.getResources().pattern.matches(system.getOperator()));
    ds.solve(system);
    EXPECT_FALSE(ds.success());
    EXPECT_EQ((system.getSolution().toHost() - previous).norm(), Real(0));
    // Recover using the retained analysis and the same graph.
    ::Kokkos::deep_copy(system.getOperator().getValues(), values);
    check();
    ds.clear(Solver::Factorization::Symbolic);
    EXPECT_FALSE(ds.getResources().pattern.matches(system.getOperator()));
    check();
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

#ifdef RODIN_USE_ROCSOLVER
  TEST(Kokkos_Variational_DeviceSystem, RocSolverRF)
  {
    checkRocSolverRF();
  }
#endif
}
