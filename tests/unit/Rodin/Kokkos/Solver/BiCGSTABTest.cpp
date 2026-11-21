/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Kokkos/Solver.h>
#include <Rodin/Solver/SparseLU.h>

#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;
  using BiCGMatrix = Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>;
  using BiCGVector = Rodin::Kokkos::Math::Vector<Real, TestBackend>;

  template <class Vector = BiCGVector>
  struct BiCGFixture
  {
    using Space = decltype(P1(std::declval<const DeviceMesh&>()));
    using Trial = decltype(TrialFunction(std::declval<Space&>()));
    using Test = decltype(TestFunction(std::declval<Space&>()));
    using System = Math::LinearSystem<BiCGMatrix, Vector>;
    DeviceMesh mesh{TestContext{}, makeTriangleMesh()};
    Space fes{mesh};
    Trial u{fes};
    Test v{fes};
    Problem<System, Trial, Test> problem{u, v};
  };

  /// @brief Small independent CSR fixtures, uploaded without FEM assembly.
  BiCGMatrix biCGMatrix(size_t rows, size_t cols,
    const std::vector<Math::SparseIndex>& starts,
    const std::vector<Math::SparseIndex>& columns, const std::vector<Real>& values)
  {
    const auto execution = TestContext{}.getExecutionSpace();
    return BiCGMatrix(rows, cols,
      Rodin::Kokkos::Internal::upload(execution, starts, "BiCGSTAB test rows"),
      Rodin::Kokkos::Internal::upload(execution, columns, "BiCGSTAB test columns"),
      Rodin::Kokkos::Internal::upload(execution, values, "BiCGSTAB test values"));
  }

  /// @brief A nonsymmetric, strictly diagonally dominant operator.
  BiCGMatrix biCGNonsymmetric()
  {
    return biCGMatrix(4, 4, {0, 2, 5, 8, 10}, {0, 1, 0, 1, 2, 1, 2, 3, 0, 3},
      {5, 1, -2, 6, 1, -1, 7, 2, 1, 8});
  }

  template <class Vector>
  Math::Vector<Real> biCGHost(const Vector& v)
  {
    if constexpr (requires { v.toHost(); })
      return v.toHost();
    else
      return v;
  }

  /// @brief Both vector adapters solve accurately after matrix and size changes.
  template <class Vector>
  void checkBiCGAdapters()
  {
    BiCGFixture<Vector> f;
    typename BiCGFixture<Vector>::System system;
    system.getOperator() = biCGNonsymmetric();
    Math::Vector<Real> exact(4);
    exact << 1, -2, 3, 0.5;
    const auto a = system.getOperator().toHost();
    Math::Vector<Real> b = a * exact;
    system.getVector() = b;
    Solver::BiCGSTAB solver(f.problem);
    static_assert(std::is_same_v<decltype(solver), Solver::BiCGSTAB<decltype(system)>>);
    EXPECT_FALSE(solver.success());
    solver.setTolerance(1e-12).setMaxIterations(100).solve(system);
    ASSERT_TRUE(solver.success());
    EXPECT_GT(solver.getIterationNumber(), 0u);
    const auto x = biCGHost(system.getSolution());
    const Real residual = (b - a * x).norm() / b.norm();
    EXPECT_LE(residual, 1e-12);
    EXPECT_NEAR(solver.getError(), residual, 1e-14);
    EXPECT_LE((x - exact).norm(), 1e-10);
    EXPECT_EQ((biCGHost(system.getVector()) - b).norm(), Real(0));
    const auto* storage = system.getSolution().data();
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_EQ(system.getSolution().data(), storage);

    system.getOperator() *= Real(3);
    b = Real(3) * a * exact;
    system.getVector() = b;
    system.getSolution().setZero();
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    EXPECT_LE((biCGHost(system.getSolution()) - exact).norm(), 1e-10);
    EXPECT_EQ(system.getSolution().data(), storage);
    system.getVector().setZero();
    solver.solve(system);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_EQ(solver.getError(), Real(0));
    EXPECT_EQ(system.getSolution().norm(), Real(0));
    system.getOperator() = biCGMatrix(2, 2, {0, 1, 2}, {0, 1}, {3, 4});
    Math::Vector<Real> smaller(2);
    smaller << 3, -8;
    system.getVector() = smaller;
    solver.solve(system);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(system.getSolution().size(), 2u);
    Math::Vector<Real> smallerExact(2);
    smallerExact << 1, -2;
    EXPECT_LE((biCGHost(system.getSolution()) - smallerExact).norm(), 1e-12);
  }

  TEST(Kokkos_Solver_BiCGSTAB, BackendVectors)
  {
    checkBiCGAdapters<BiCGVector>();
  }

  TEST(Kokkos_Solver_BiCGSTAB, HostVectors)
  {
    checkBiCGAdapters<Math::Vector<Real>>();
  }

  TEST(Kokkos_Solver_BiCGSTAB, SharedHostAdapterSupportsCGAndGMRES)
  {
    BiCGFixture<Math::Vector<Real>> f;
    auto& s = f.problem.getLinearSystem();
    s.getOperator() = biCGMatrix(3, 3, {0, 1, 2, 3}, {0, 1, 2}, {2, 3, 4});
    s.getVector().resize(3);
    s.getVector() << 2, 3, 4;
    Solver::CG cg(f.problem);
    cg.setTolerance(1e-12).solve(s);
    ASSERT_TRUE(cg.success());
    EXPECT_LE((s.getSolution().array() - Real(1)).matrix().norm(), 1e-12);
    s.getSolution().setZero();
    Solver::GMRES gmres(f.problem);
    gmres.setTolerance(1e-12).solve(s);
    ASSERT_TRUE(gmres.success());
    EXPECT_LE((s.getSolution().array() - Real(1)).matrix().norm(), 1e-12);
  }

  TEST(Kokkos_Solver_BiCGSTAB, ExactGuessAndZeroIterationLimit)
  {
    BiCGFixture<> f;
    auto& s = f.problem.getLinearSystem();
    s.getOperator() = biCGNonsymmetric();
    Math::Vector<Real> exact(4);
    exact << 1, 2, 3, 4;
    s.getVector() = Math::Vector<Real>(s.getOperator().toHost() * exact);
    s.getSolution() = exact;
    Solver::BiCGSTAB solver(f.problem);
    solver.setTolerance(1e-12).setMaxIterations(0).solve(s);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    s.getSolution().setZero();
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_NEAR(solver.getError(), Real(1), 1e-15);
    solver.setMaxIterations(1).solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 1u);
    const auto b = s.getVector().toHost();
    EXPECT_NEAR(solver.getError(),
      (b - s.getOperator().toHost() * s.getSolution().toHost()).norm() / b.norm(), 1e-14);
  }

  TEST(Kokkos_Solver_BiCGSTAB, DiagonalHappyBreakdownAndEmptySystem)
  {
    BiCGFixture<> f;
    auto& s = f.problem.getLinearSystem();
    // Exact Jacobi makes the intermediate residual zero after one alpha step.
    s.getOperator() = biCGMatrix(3, 3, {0, 1, 2, 3}, {0, 1, 2}, {2, 3, 4});
    Math::Vector<Real> b(3);
    b << 2, 6, 12;
    s.getVector() = b;
    Solver::BiCGSTAB solver(f.problem);
    solver.setTolerance(1e-12).solve(s);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 1u);
    Math::Vector<Real> exact(3);
    exact << 1, 2, 3;
    EXPECT_LE((s.getSolution().toHost() - exact).norm(), 1e-12);
    typename BiCGFixture<>::System empty;
    solver.solve(empty);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_EQ(solver.getError(), 0);
    EXPECT_EQ(empty.getSolution().size(), 0u);
  }

  TEST(Kokkos_Solver_BiCGSTAB, VerifiesTheTrueResidual)
  {
    BiCGFixture<> f;
    auto& s = f.problem.getLinearSystem();
    // Around 3.3/49, no representable double makes 49*x round to 3.3.
    // The recursive intermediate residual can nevertheless round to zero.
    s.getOperator() = biCGMatrix(1, 1, {0, 1}, {0}, {49});
    Math::Vector<Real> b(1);
    b << 3.3;
    s.getVector() = b;
    Solver::BiCGSTAB solver(f.problem);
    solver.setTolerance(1e-30).setMaxIterations(6).solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_GT(solver.getError(), 1e-30);
    EXPECT_LE(solver.getIterationNumber(), 6u);
    EXPECT_TRUE(s.getSolution().allFinite());
  }

  TEST(Kokkos_Solver_BiCGSTAB, BreakdownAndNonfiniteInputsFail)
  {
    BiCGFixture<> f;
    auto& s = f.problem.getLinearSystem();
    // A nonsingular rotation has r.dot(A*r) == 0: BiCGSTAB breaks down.
    s.getOperator() = biCGMatrix(2, 2, {0, 1, 2}, {1, 0}, {1, -1});
    Math::Vector<Real> b(2);
    b << 1, 0;
    s.getVector() = b;
    Solver::BiCGSTAB solver(f.problem);
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_NEAR(solver.getError(), 1, 1e-15);
    EXPECT_TRUE(s.getSolution().allFinite());
    // A rank-one operator makes the stabilization product vanish.
    s.getOperator() = biCGMatrix(2, 2, {0, 2, 2}, {0, 1}, {1, 1});
    b << 1, 1;
    s.getVector() = b;
    s.getSolution().setZero();
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_TRUE(s.getSolution().allFinite());
    EXPECT_NEAR(solver.getError(), 1, 1e-15);
    // This nonsingular operator gives a zero omega despite a nonzero product.
    s.getOperator() = biCGMatrix(2, 2, {0, 2, 3}, {0, 1, 0}, {1, 1, 1});
    b << 1, 0;
    s.getVector() = b;
    s.getSolution().setZero();
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 1u);
    EXPECT_TRUE(s.getSolution().allFinite());
    EXPECT_NEAR(solver.getError(), 1, 1e-15);
    s.getOperator() = biCGMatrix(2, 2, {0, 0, 0}, {}, {});
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_TRUE(s.getSolution().allFinite());
    s.getOperator() = biCGMatrix(
      2, 2, {0, 1, 2}, {0, 1}, {std::numeric_limits<Real>::quiet_NaN(), Real(1)});
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_TRUE(std::isinf(solver.getError()));
    s.getOperator() = biCGMatrix(2, 2, {0, 1, 2}, {0, 1}, {1, 1});
    b(0) = std::numeric_limits<Real>::infinity();
    s.getVector() = b;
    solver.solve(s);
    EXPECT_FALSE(solver.success());
    EXPECT_TRUE(std::isinf(solver.getError()));
  }

  TEST(Kokkos_Solver_BiCGSTAB, InvalidRequestsAndIndependentCopies)
  {
    BiCGFixture<> f;
    auto& s = f.problem.getLinearSystem();
    Solver::BiCGSTAB solver(f.problem);
    EXPECT_ANY_THROW(solver.setTolerance(-1));
    EXPECT_ANY_THROW(solver.setTolerance(std::numeric_limits<Real>::infinity()));
    EXPECT_ANY_THROW(solver.setTolerance(std::numeric_limits<Real>::quiet_NaN()));
    s.getOperator() = biCGMatrix(1, 2, {0, 1}, {0}, {1});
    s.getVector().setZero(1);
    EXPECT_ANY_THROW(solver.solve(s));
    s.getOperator() = biCGNonsymmetric();
    EXPECT_ANY_THROW(solver.solve(s));
    Math::Vector<Real> b(4);
    b << 1, 2, 3, 4;
    s.getVector() = b;
    solver.setTolerance(1e-12).setMaxIterations(0);
    Solver::BiCGSTAB copy(solver);
    solver.setMaxIterations(100).solve(s);
    ASSERT_TRUE(solver.success());
    const auto previous = s.getSolution().toHost();
    typename BiCGFixture<>::System other;
    other.getOperator() = s.getOperator();
    other.getVector() = b;
    copy.solve(other);
    EXPECT_FALSE(copy.success());
    EXPECT_EQ(copy.getIterationNumber(), 0u);
    Solver::BiCGSTAB moved(std::move(copy));
    moved.setMaxIterations(100).solve(other);
    ASSERT_TRUE(moved.success());
    EXPECT_LE((other.getSolution().toHost() - previous).norm(), 1e-10);
    EXPECT_EQ((s.getSolution().toHost() - previous).norm(), 0);
    std::unique_ptr<decltype(solver)> cloned(solver.copy());
    other.getSolution().setZero();
    cloned->solve(other);
    EXPECT_TRUE(cloned->success());
  }

  /// @brief The unchanged form-language API solves nonsymmetric FEM operators.
  void checkBiCGForm()
  {
    constexpr size_t Nodes = 14;
    auto source = HostMesh::UniformGrid(Polytope::Type::Triangle, {Nodes, Nodes});
    source.scale(Real(1) / (Nodes - 1));
    source.getConnectivity().compute(1, 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    const VectorFunction beta{20, 10};
    Problem pb(u, v);
    pb = Integral(Grad(u), Grad(v)) + Integral(Dot(beta, Grad(u)), v) + Integral(u, v) -
      Integral(RealFunction(1), v) + DirichletBC(u, RealFunction(0));
    Solver::BiCGSTAB solver(pb);
    solver.setTolerance(1e-11).setMaxIterations(1000).solve();
    ASSERT_TRUE(solver.success());
    auto& s = pb.getLinearSystem();
    const auto a = s.getOperator().toHost();
    const auto b = s.getVector().toHost();
    const auto x = s.getSolution().toHost();
    EXPECT_LE((b - a * x).norm() / b.norm(), 1e-11);
    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(u), decltype(v)> referenceProblem(u, v);
    auto& referenceSystem = referenceProblem.getLinearSystem();
    referenceSystem.getOperator() = a;
    referenceSystem.getVector() = b;
    Solver::SparseLU reference(referenceProblem);
    reference.solve(referenceSystem);
    ASSERT_TRUE(reference.success());
    const auto exact = referenceSystem.getSolution();
    EXPECT_LE((x - exact).norm() / exact.norm(), 1e-9);
    EXPECT_LE((u.getSolution().getData().toHost() - exact).norm() / exact.norm(), 1e-9);
    const auto* storage = s.getSolution().data();
    solver.solve(s);
    EXPECT_TRUE(solver.success());
    EXPECT_EQ(solver.getIterationNumber(), 0u);
    EXPECT_EQ(s.getSolution().data(), storage);
  }

  TEST(Kokkos_Solver_BiCGSTAB, VariationalFormAndCTAD)
  {
    checkBiCGForm();
  }
}
