/*
 * @file ProblemTest.cpp
 * @brief Unit tests for Problem, SparseProblem, and DenseProblem classes.
 */
#include <gtest/gtest.h>
#include <type_traits>
#include <array>
#include <algorithm>
#include <utility>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/**
 * @brief Checks mixed problem construction without legacy tuple assemblers.
 *
 * A cellwise P0 field and an H1 field have different block sizes. Constant
 * manufactured fields solve a positive-definite reaction system with diagonal
 * coefficients one and two. Complex cases use conjugate off-diagonal
 * coefficients: @f$R=\begin{pmatrix}1&\beta\\\overline\beta&2\end{pmatrix}@f$,
 * with @f$\beta=1/4@f$ or @f$\beta=1/4+i/8@f$. Since
 * @f$\det R=2-|\beta|^2>0@f$, CG is applicable to both coefficient fields.
 * Three H1 degrees and all seven cell families exercise both
 * sparse and dense system representations, including repeated assembly and
 * solution scatter. Copy/move construction must retain the
 * owned formulation, block offsets and assembled system; explicit reassembly
 * must reproduce the original operator and load exactly.
 * Separate scalar mass forms determine every matrix block.
 * The dimensionless true-residual and coefficient budgets are @f$10^{-11}@f$;
 * CG's requested tolerance and iteration cap are @f$10^{-13}@f$ and 20,000,
 * as in native convergence solves. Order twelve integrates degree-three basis products
 * and retains the established rational-pyramid quadrature policy.
 *
 * @tparam Scalar Real or complex coefficient field.
 * @tparam Matrix Sparse or dense native operator representation.
 */
template <class Scalar, class Matrix>
void checkMixedBlockProblem()
{
  using System = Math::LinearSystem<Matrix, Math::Vector<Scalar>>;
  constexpr Real algebraicTolerance = 1e-11, solverTolerance = 1e-13;
  constexpr size_t maximumIterations = 20000, quadratureOrder = 12;
  const std::array<std::pair<Polytope::Type, size_t>, 7> geometries{
    {{Polytope::Type::Segment, 1}, {Polytope::Type::Triangle, 2},
      {Polytope::Type::Quadrilateral, 2}, {Polytope::Type::Tetrahedron, 3},
      {Polytope::Type::Pyramid, 3}, {Polytope::Type::Hexahedron, 3},
      {Polytope::Type::Wedge, 3}}};
  const auto checkDegree = [&]<size_t K>() {
    for (const auto& [geometry, dimension] : geometries)
    {
      SCOPED_TRACE(
        ::testing::Message() << "geometry=" << int(geometry) << " degree=" << K);
      Array<size_t> shape(dimension);
      std::fill(shape.begin(), shape.end(), 2);
      auto mesh = LocalMesh::UniformGrid(geometry, shape);
      for (size_t from = 0; from <= dimension; ++from)
      {
        for (size_t to = 0; to <= dimension; ++to)
        {
          if (from != to)
            mesh.getConnectivity().compute(from, to);
        }
      }
      P0<Scalar> low(mesh);
      H1<K, Scalar> high(std::integral_constant<size_t, K>{}, mesh);
      TrialFunction u(low);
      TrialFunction w(high);
      TestFunction v(low);
      TestFunction z(high);
      const auto constant = [](Scalar value) {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(value);
        else
          return RealFunction(value);
      };
      const Scalar first = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(1, 1);
        else
          return Scalar(1);
      }();
      const Scalar second = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(2, -1);
        else
          return Scalar(2);
      }();
      const Scalar coupling = [] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(0.25, 0.125);
        else
          return Scalar(0.25);
      }();
      const Scalar reverse = [&] {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return std::conj(coupling);
        else
          return coupling;
      }();
      u.getSolution() = constant(first);
      w.getSolution() = constant(second);
      const Math::Vector<Scalar> exactU = u.getSolution().getData();
      const Math::Vector<Scalar> exactW = w.getSolution().getData();
      Problem<System, decltype(u), decltype(v), decltype(w), decltype(z)> problem(
        u, v, w, z);
      auto aUU = Integral(u, v);
      auto aWW = Integral(Real(2) * w, z);
      auto aUW = Integral(coupling * w, v);
      auto aWU = Integral(reverse * u, z);
      auto bU = Integral(constant(first + coupling * second), v);
      auto bW = Integral(constant(reverse * first + Scalar(2) * second), z);
      aUU.setOrder(quadratureOrder);
      aWW.setOrder(quadratureOrder);
      aUW.setOrder(quadratureOrder);
      aWU.setOrder(quadratureOrder);
      bU.setOrder(quadratureOrder);
      bW.setOrder(quadratureOrder);
      problem = aUU + aUW + aWU + aWW - bU - bW;
      problem.assemble();
      const auto& trial = problem.getTrialOffsets();
      const auto& test = problem.getTestOffsets();
      EXPECT_EQ(trial[0], 0);
      EXPECT_EQ(trial[1], low.getSize());
      EXPECT_EQ(test, trial);
      auto& system = problem.getLinearSystem();
      ASSERT_EQ(system.getOperator().rows(), low.getSize() + high.getSize());
      const Math::Vector<Scalar> exact = system.getSolution();
      EXPECT_EQ((exact.segment(trial[0], low.getSize()) - exactU).norm(), 0);
      EXPECT_EQ((exact.segment(trial[1], high.getSize()) - exactW).norm(), 0);
      BilinearForm lowMass(u, v);
      BilinearForm highMass(w, z);
      BilinearForm crossMass(w, v);
      auto mUU = Integral(u, v);
      auto mWW = Integral(w, z);
      auto mUW = Integral(w, v);
      mUU.setOrder(quadratureOrder);
      mWW.setOrder(quadratureOrder);
      mUW.setOrder(quadratureOrder);
      lowMass = mUU;
      highMass = mWW;
      crossMass = mUW;
      lowMass.assemble();
      highMass.assemble();
      crossMass.assemble();
      Math::Matrix<Scalar> expected =
        Math::Matrix<Scalar>::Zero(exact.size(), exact.size());
      expected.block(test[0], trial[0], low.getSize(), low.getSize()) =
        lowMass.getOperator();
      expected.block(test[1], trial[1], high.getSize(), high.getSize()) =
        Scalar(2) * highMass.getOperator();
      expected.block(test[0], trial[1], low.getSize(), high.getSize()) =
        coupling * crossMass.getOperator();
      expected.block(test[1], trial[0], high.getSize(), low.getSize()) =
        reverse * crossMass.getOperator().adjoint();
      const Math::Matrix<Scalar> assembled = system.getOperator();
      const Real matrixScale = std::max(Real(1), expected.norm());
      EXPECT_LT((assembled - expected).norm() / matrixScale, algebraicTolerance);
      EXPECT_LT(
        (assembled - assembled.adjoint()).norm() / matrixScale, algebraicTolerance);
      const Math::Vector<Scalar> load = system.getVector();
      const Real rhsScale = std::max(Real(1), load.norm());
      EXPECT_LT(
        (system.getOperator() * exact - load).norm() / rhsScale, algebraicTolerance);
      problem.assemble();
      EXPECT_EQ((Math::Matrix<Scalar>(system.getOperator()) - assembled).norm(), 0);
      EXPECT_EQ((system.getVector() - load).norm(), 0);
      const auto checkTransferred = [&](auto& transferred) {
        ASSERT_EQ(transferred.getBody().getLocalBFIs().size(),
          problem.getBody().getLocalBFIs().size());
        ASSERT_EQ(
          transferred.getBody().getLFIs().size(), problem.getBody().getLFIs().size());
        EXPECT_EQ(transferred.getTrialOffsets(), trial);
        EXPECT_EQ(transferred.getTestOffsets(), test);
        auto& transferredSystem = transferred.getLinearSystem();
        EXPECT_EQ(
          (Math::Matrix<Scalar>(transferredSystem.getOperator()) - assembled).norm(), 0);
        EXPECT_EQ((transferredSystem.getVector() - load).norm(), 0);
        transferred.assemble();
        EXPECT_EQ(
          (Math::Matrix<Scalar>(transferredSystem.getOperator()) - assembled).norm(), 0);
        EXPECT_EQ((transferredSystem.getVector() - load).norm(), 0);
      };
      auto copied(problem);
      {
        SCOPED_TRACE("copy construction");
        checkTransferred(copied);
      }
      auto moved(std::move(copied));
      {
        SCOPED_TRACE("move construction");
        checkTransferred(moved);
      }
      system.getSolution().setZero();
      Solver::CG solver(problem);
      solver.setTolerance(solverTolerance).setMaxIterations(maximumIterations);
      problem.solve(solver);
      ASSERT_TRUE(solver.success());
      EXPECT_LT((system.getOperator() * system.getSolution() - load).norm() / rhsScale,
        algebraicTolerance);
      EXPECT_LT(
        (u.getSolution().getData() - exactU).norm() / std::max(Real(1), exactU.norm()),
        algebraicTolerance);
      EXPECT_LT(
        (w.getSolution().getData() - exactW).norm() / std::max(Real(1), exactW.norm()),
        algebraicTolerance);
    }
  };
  checkDegree.template operator()<1>();
  checkDegree.template operator()<2>();
  checkDegree.template operator()<3>();
}

TEST(Rodin_Variational_Problem, MixedRealSparseP0H1)
{
  checkMixedBlockProblem<Real, Math::SparseMatrix<Real>>();
}

TEST(Rodin_Variational_Problem, MixedComplexSparseP0H1)
{
  checkMixedBlockProblem<Complex, Math::SparseMatrix<Complex>>();
}

TEST(Rodin_Variational_Problem, MixedRealDenseP0H1)
{
  checkMixedBlockProblem<Real, Math::Matrix<Real>>();
}

TEST(Rodin_Variational_Problem, MixedComplexDenseP0H1)
{
  checkMixedBlockProblem<Complex, Math::Matrix<Complex>>();
}

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

TEST(Rodin_Variational_Problem, PreassembledFormSignsP1P2P3SparseAndDense)
{
  const auto check = []<size_t Order, class Matrix>() {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    for (size_t from = 0; from <= 2; ++from)
    {
      for (size_t to = 0; to <= 2; ++to)
      {
        if (from != to)
          mesh.getConnectivity().compute(from, to);
      }
    }
    H1 fes(std::integral_constant<size_t, Order>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    using System = Math::LinearSystem<Matrix, Math::Vector<Real>>;
    BilinearForm<std::remove_reference_t<decltype(u.getSolution())>, decltype(fes),
      decltype(fes), Matrix>
      metric(u, v);
    LinearForm load(v);
    metric = Integral(u, v);
    load = Integral(RealFunction(1), v);
    metric.assemble();
    load.assemble();
    const Matrix matrix = metric.getOperator();
    const Math::Vector<Real> vector = load.getVector();
    Problem<System, decltype(u), decltype(v)> problem(u, v);
    ProblemBase<System>& base = problem;
    base += metric;
    base -= load;
    problem.assemble();
    EXPECT_LT((problem.getLinearSystem().getOperator() - matrix).norm(), Real(1e-12));
    EXPECT_LT((problem.getLinearSystem().getVector() - vector).norm(), Real(1e-12));
    base -= metric;
    base += load;
    problem.assemble();
    EXPECT_LT(problem.getLinearSystem().getOperator().norm(), Real(1e-12));
    EXPECT_LT(problem.getLinearSystem().getVector().norm(), Real(1e-12));
    base -= metric;
    base += load;
    metric.getOperator() *= 7;
    load.getVector() *= 9;
    problem.assemble();
    EXPECT_LT((problem.getLinearSystem().getOperator() + matrix).norm(), Real(1e-12));
    EXPECT_LT((problem.getLinearSystem().getVector() + vector).norm(), Real(1e-12));
  };
  check.template operator()<1, Math::SparseMatrix<Real>>();
  check.template operator()<2, Math::SparseMatrix<Real>>();
  check.template operator()<3, Math::SparseMatrix<Real>>();
  check.template operator()<1, Math::Matrix<Real>>();
  check.template operator()<2, Math::Matrix<Real>>();
  check.template operator()<3, Math::Matrix<Real>>();
}

TEST(Rodin_Variational_Problem, PreassembledFormsInvalidateEverySolve)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  P1 fes(mesh);
  TrialFunction u(fes);
  TestFunction v(fes);
  BilinearForm metric(u, v);
  LinearForm load(v);
  metric = Integral(u, v);
  load = Integral(RealFunction(1), v);
  metric.assemble();
  load.assemble();
  Problem problem(u, v);
  problem += metric;
  problem -= load;
  Solver::CG solver(problem);
  const auto solve = [&](Real expected) {
    problem.solve(solver);
    EXPECT_LT(
      (u.getSolution().getData().array() - expected).matrix().norm(), Real(1e-10));
  };
  solve(1);
  problem += metric;
  solve(0.5);
  problem -= metric;
  solve(1);
  problem += load;
  solve(0);
  problem -= load;
  solve(1);
}

TEST(Rodin_Variational_Problem, CompoundAssignmentMatchesBodyP1P2P3)
{
  const auto check = []<size_t Order>() {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    for (size_t from = 0; from <= 2; ++from)
    {
      for (size_t to = 0; to <= 2; ++to)
      {
        if (from != to)
          mesh.getConnectivity().compute(from, to);
      }
    }
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
