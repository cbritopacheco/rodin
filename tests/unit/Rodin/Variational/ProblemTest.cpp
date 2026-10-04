/*
 * @file ProblemTest.cpp
 * @brief Unit tests for Problem, SparseProblem, and DenseProblem classes.
 */
#include <gtest/gtest.h>
#include <type_traits>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

TEST(Rodin_Variational_Problem, BodylessAdaptersRejectCompoundAssignment)
{
  using System = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
  class Adapter final : public ProblemBase<System>
  {
    public:
      Adapter& operator=(const ProblemBodyType&) override
      {
        return *this;
      }
      Adapter& assemble() override
      {
        return *this;
      }
      void solve(Solver::LinearSolverBase<System>&) override {}
      System& getLinearSystem() override
      {
        return m_system;
      }
      const System& getLinearSystem() const override
      {
        return m_system;
      }
      Adapter* copy() const noexcept override
      {
        return new Adapter(*this);
      }

    private:
      System m_system;
  };
  Adapter adapter;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  P1 fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  EXPECT_THROW(adapter.getBody(), Alert::Exception);
  EXPECT_THROW(adapter += Integral(u, v), Alert::Exception);
  EXPECT_THROW(adapter -= Integral(u, v), Alert::Exception);
}

TEST(Rodin_Variational_Problem, CompoundAssignmentAcceptsPreassembledMetric)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  P1 fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  BilinearForm metric(u, v);
  metric = Integral(u, v);
  metric.assemble();
  Problem expected(u, v), actual(u, v);
  expected = metric - Integral(RealFunction(1), v);
  actual += metric;
  actual -= Integral(RealFunction(1), v);
  expected.assemble();
  actual.assemble();
  Math::SparseMatrix<Real> difference =
    actual.getLinearSystem().getOperator() - expected.getLinearSystem().getOperator();
  EXPECT_LT(difference.norm(), Real(1e-12));
  EXPECT_LT(
    (actual.getLinearSystem().getVector() - expected.getLinearSystem().getVector())
      .norm(),
    Real(1e-12));
}

TEST(Rodin_Variational_Problem, CompoundAssignmentMatchesBodyP1P2P3)
{
  const auto check = []<size_t Order>() {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    for (size_t from = 0; from <= 2; ++from)
      for (size_t to = 0; to <= 2; ++to)
        if (from != to)
          mesh.getConnectivity().compute(from, to);
    H1 fes(std::integral_constant<size_t, Order>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem expected(u, v), actual(u, v);
    RealFunction load(2), boundary(1);
    expected = Integral(u, v) + Integral(Grad(u), Grad(v)) - Integral(load, v) +
      DirichletBC(u, boundary);
    using System = std::remove_reference_t<decltype(actual.getLinearSystem())>;
    ProblemBase<System>& base = actual;
    actual += Integral(u, v);
    base += Integral(Grad(u), Grad(v));
    base -= Integral(load, v);
    actual += Integral(load, v);
    actual -= Integral(load, v);
    actual += Integral(u, v);
    actual -= Integral(u, v);
    base += DirichletBC(u, boundary);
    expected.assemble();
    actual.assemble();
    Math::SparseMatrix<Real> difference =
      actual.getLinearSystem().getOperator() - expected.getLinearSystem().getOperator();
    EXPECT_LT(difference.norm(), Real(1e-12));
    EXPECT_LT(
      (actual.getLinearSystem().getVector() - expected.getLinearSystem().getVector())
        .norm(),
      Real(1e-12));
    const auto condition = DirichletBC(u, boundary);
    static_assert(!requires(ProblemBase<System>& p) { p -= condition; });
  };
  check.template operator()<1>();
  check.template operator()<2>();
  check.template operator()<3>();
}

TEST(Rodin_Variational_Problem, CompoundAssignmentInvalidatesSolveAndRetainsCopies)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  P1 fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  RealFunction one(1);
  Problem problem(u, v);
  problem += Integral(u, v);
  problem -= Integral(one, v);
  Solver::CG solver(problem);
  problem.solve(solver);
  EXPECT_LT((u.getSolution().getData().array() - Real(1)).matrix().norm(), Real(1e-10));
  auto copy = problem;
  problem += Integral(u, v);
  problem.solve(solver);
  EXPECT_LT((u.getSolution().getData().array() - Real(0.5)).matrix().norm(), Real(1e-10));
  problem -= Integral(u, v);
  problem.solve(solver);
  EXPECT_LT((u.getSolution().getData().array() - Real(1)).matrix().norm(), Real(1e-10));
  copy.assemble();
  Math::SparseMatrix<Real> difference =
    copy.getLinearSystem().getOperator() - problem.getLinearSystem().getOperator();
  EXPECT_LT(difference.norm(), Real(1e-12));
}

/// @brief Verifies construction from trial test for variational problem.
TEST(Rodin_Variational_Problem, ConstructionFromTrialTest)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  Problem problem(u, v);
  // Construction succeeds
  SUCCEED();
}

/// @brief Verifies assemble poisson for variational problem by checking form assembly.
TEST(Rodin_Variational_Problem, AssemblePoisson)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  RealFunction f(1.0);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, RealFunction(0.0));
  problem.assemble();

  // Assemble succeeds
  SUCCEED();
}

/// @brief Verifies solve poisson for variational problem by checking true predicates, solver behavior.
TEST(Rodin_Variational_Problem, SolvePoisson)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  RealFunction f(1.0);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, RealFunction(0.0));

  Solver::CG(problem).solve();

  // Solution should be non-trivial
  const auto& sol = u.getSolution();
  bool hasNonZero = false;
  for (Index i = 0; i < sol.getSize(); ++i)
  {
    if (std::abs(sol.getData()(i)) > 1e-15)
    {
      hasNonZero = true;
      break;
    }
  }
  EXPECT_TRUE(hasNonZero);
}

/// @brief Verifies solution size for variational problem by checking exact expected values, solver behavior.
TEST(Rodin_Variational_Problem, SolutionSize)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
  mesh.getConnectivity().compute(1, 2);
  P1 Vh(mesh);
  TrialFunction u(Vh);
  TestFunction v(Vh);

  RealFunction f(1.0);

  Problem problem(u, v);
  problem = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, RealFunction(0.0));

  Solver::CG(problem).solve();

  const auto& sol = u.getSolution();
  EXPECT_EQ(sol.getSize(), Vh.getSize());
}
