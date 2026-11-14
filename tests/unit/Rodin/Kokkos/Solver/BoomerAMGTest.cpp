/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Kokkos/Solver/BoomerAMG.h>

#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /// @brief System whose matrix stays in the memory of the test backend.
  using DeviceLinearSystem =
    Math::LinearSystem<Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>,
      Rodin::Kokkos::Math::Vector<Real, TestBackend>>;

  /// @brief Unit square split into triangles, on the test backend.
  inline DeviceMesh makeSquare(size_t n)
  {
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
    source.getConnectivity().compute(1, 2);
    return DeviceMesh(TestContext{}, source);
  }

  /// @brief Relative distance between two vectors.
  inline Real distance(const Math::Vector<Real>& a, const Math::Vector<Real>& b)
  {
    return (a - b).norm() / b.norm();
  }

  /// @brief Scales the values of a device matrix.
  void scale(DeviceLinearSystem& system, Real factor)
  {
    const auto values = system.getOperator().getValues();
    ::Kokkos::parallel_for(
      ::Kokkos::RangePolicy<
        Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>::ExecutionSpace>(
        0, values.extent(0)),
      KOKKOS_LAMBDA(const size_t i) { values(i) *= factor; });
    ::Kokkos::fence();
  }

  /// @brief BoomerAMG preconditioned CG matches a direct solve, restarts from
  /// its solution, reuses its hierarchy, and reports a solve that does not
  /// converge.
  void checkBoomerAMG()
  {
    const auto mesh = makeSquare(32);
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
    Solver::SparseLU(host).solve();
    const Math::Vector<Real> reference = u.getSolution().getData();
    // The assembly seeds the initial guess with the current solution.
    u.getSolution() = Real(0);

    Problem<DeviceLinearSystem, decltype(u), decltype(v)> device(u, v);
    device = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(f, v) +
      DirichletBC(u, RealFunction(1));
    Solver::BoomerAMG solver(device);
    static_assert(
      std::is_same_v<decltype(solver), Solver::BoomerAMG<DeviceLinearSystem>>);
    solver.setTolerance(1e-12).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_GT(solver.getIterationNumber(), 0u);
    // Multigrid converges in far fewer iterations than Jacobi CG.
    EXPECT_LT(solver.getIterationNumber(), 30u);
    EXPECT_LE(solver.getError(), 1e-12);
    EXPECT_LE(distance(u.getSolution().getData(), reference), 1e-10);

    // The solution of the system is the initial guess of the next solve, which
    // hypre checks after one iteration.
    auto& axb = device.getLinearSystem();
    solver.solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 1u);
    EXPECT_LE(distance(axb.getSolution(), reference), 1e-10);

    // The retained hierarchy preconditions the conjugate gradient on a
    // matrix that changed since its setup. Doubling the values is exact, so
    // the solution is half the reference.
    ASSERT_TRUE(solver.isSetup());
    scale(axb, 2);
    axb.getSolution().resize(0);
    solver.solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_LT(solver.getIterationNumber(), 30u);
    EXPECT_LE(distance(2 * axb.getSolution(), reference), 1e-10);

    // An explicit setup builds the hierarchy again from the current matrix.
    scale(axb, 0.5);
    solver.setup(axb);
    ASSERT_TRUE(solver.isSetup());
    axb.getSolution().resize(0);
    solver.solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_LT(solver.getIterationNumber(), 30u);
    EXPECT_LE(distance(axb.getSolution(), reference), 1e-10);

    // A copy leaves the hierarchy to the original, which releases it on
    // request.
    {
      Solver::BoomerAMG<DeviceLinearSystem> copy(solver);
      EXPECT_FALSE(copy.isSetup());
    }
    solver.clear();
    EXPECT_FALSE(solver.isSetup());

    // An iteration limit that is too small fails without raising.
    axb.getSolution().resize(0);
    solver.setTolerance(1e-14).setMaxIterations(1).solve(axb);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 1u);

    // The coarsening settings apply, since changing them releases the
    // hierarchy.
    ASSERT_TRUE(solver.isSetup());
    solver.setTolerance(1e-12).setMaxIterations(100).setStrongThreshold(0.5);
    EXPECT_FALSE(solver.isSetup());
    solver.setAggressiveLevels(1).setMaxInterpolationElements(2).solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_LE(distance(axb.getSolution(), reference), 1e-10);

    // A zero right-hand side gives a zero solution.
    axb.getVector().setZero();
    solver.solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_EQ(axb.getSolution().norm(), Real(0));
  }

  TEST(Kokkos_Solver_BoomerAMG, CG)
  {
    checkBoomerAMG();
  }
}
