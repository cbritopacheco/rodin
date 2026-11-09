/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Kokkos/Solver.h>

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

  /// @brief CG on the device matches a direct solve, and restarts from its
  /// solution without iterating.
  void checkCG()
  {
    const auto mesh = makeSquare(24);
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
    Solver::CG solver(device);
    static_assert(std::is_same_v<decltype(solver), Solver::CG<DeviceLinearSystem>>);
    solver.setTolerance(1e-12).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_GT(solver.getIterationNumber(), 0u);
    EXPECT_LE(solver.getError(), 1e-12);
    EXPECT_LE(distance(u.getSolution().getData(), reference), 1e-10);

    // The solution of the system is the initial guess of the next solve.
    auto& axb = device.getLinearSystem();
    solver.solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_LE(distance(axb.getSolution(), reference), 1e-10);

    // An iteration limit that is too small fails without raising.
    axb.getSolution().resize(0);
    solver.setMaxIterations(2).solve(axb);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 2u);

    // A zero right-hand side gives a zero solution.
    axb.getVector().setZero();
    solver.solve(axb);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(axb.getSolution().norm(), Real(0));
  }

  /// @brief GMRES on the device matches a direct solve of an advection
  /// diffusion problem, whose matrix is not symmetric.
  void checkGMRES()
  {
    const auto mesh = makeSquare(24);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    const VectorFunction beta{20, 10};
    Problem host(u, v);
    host = Integral(Grad(u), Grad(v)) + Integral(Dot(beta, Grad(u)), v) -
      Integral(RealFunction(1), v) + DirichletBC(u, RealFunction(0));
    Solver::SparseLU(host).solve();
    const Math::Vector<Real> reference = u.getSolution().getData();
    // The assembly seeds the initial guess with the current solution.
    u.getSolution() = Real(0);

    Problem<DeviceLinearSystem, decltype(u), decltype(v)> device(u, v);
    device = Integral(Grad(u), Grad(v)) + Integral(Dot(beta, Grad(u)), v) -
      Integral(RealFunction(1), v) + DirichletBC(u, RealFunction(0));
    Solver::GMRES solver(device);
    static_assert(std::is_same_v<decltype(solver), Solver::GMRES<DeviceLinearSystem>>);
    solver.setTolerance(1e-12).setRestart(40).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_GT(solver.getIterationNumber(), 0u);
    EXPECT_LE(solver.getError(), 1e-12);
    EXPECT_LE(distance(u.getSolution().getData(), reference), 1e-9);
  }

  TEST(Kokkos_Solver_Krylov, CG)
  {
    checkCG();
  }

  TEST(Kokkos_Solver_Krylov, GMRES)
  {
    checkGMRES();
  }
}
