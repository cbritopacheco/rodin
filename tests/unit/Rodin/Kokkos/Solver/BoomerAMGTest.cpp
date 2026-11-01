/*
 *          Copyright Carlos BRITO PACHECO 2026.
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
    Math::LinearSystem<Rodin::Kokkos::SparseMatrix<Real, TestBackend>,
      Math::Vector<Real>>;

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

  /// @brief BoomerAMG preconditioned CG matches a direct solve, restarts from
  /// its solution, and reports a solve that does not converge.
  void checkBoomerAMG()
  {
    const auto mesh = makeSquare(32);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    Problem host(u, v);
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

    // An iteration limit that is too small fails without raising.
    axb.getSolution().resize(0);
    solver.setTolerance(1e-14).setMaxIterations(1).solve(axb);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 1u);

    // The coarsening settings apply.
    solver.setTolerance(1e-12)
      .setMaxIterations(100)
      .setStrongThreshold(0.5)
      .setAggressiveLevels(1)
      .setMaxInterpolationElements(2)
      .solve(axb);
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
