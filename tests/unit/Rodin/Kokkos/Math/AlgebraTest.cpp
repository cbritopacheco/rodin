/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Kokkos/Solver.h>
#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;
  using Coefficients = Rodin::Kokkos::Math::Vector<Real, TestBackend>;
  using Matrix = Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>;

  TEST(Kokkos_Math, VectorValueSemantics)
  {
    Math::Vector<Real> source(4);
    source << 1, 2, 3, 4;
    Coefficients a(source);
    Coefficients b(a);
    EXPECT_NE(a.getView().data(), b.getView().data());
    ::Kokkos::deep_copy(b.getView(), Real(7));
    EXPECT_EQ(a.toHost()(0), 1);
    EXPECT_EQ(b.toHost()(0), 7);
    b.segment(1, 2) = static_cast<const Coefficients&>(a).segment(2, 2);
    EXPECT_EQ(b.toHost()(1), 3);
    EXPECT_EQ(b.toHost()(2), 4);
    a *= Real(2);
    EXPECT_NEAR(a.norm(), 2 * source.norm(), 1e-12);
    Index i = 0;
    EXPECT_EQ(a.minCoeff(&i), 2);
    EXPECT_EQ(i, 0);
    EXPECT_EQ(a.maxCoeff(&i), 8);
    EXPECT_EQ(i, 3);
    EXPECT_TRUE(a.allFinite());
    Coefficients c(std::move(b));
    EXPECT_EQ(c.size(), 4);
  }

  TEST(Kokkos_Math, CoefficientsAndSolveStayInBackendStorage)
  {
    const auto source = makeTriangleMesh();
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    GridFunction coefficient(fes);
    static_assert(std::is_same_v<typename decltype(coefficient)::DataType, Coefficients>);
    coefficient = Real(1);
    const auto native = coefficient.getData().getView();
    const auto borrowed =
      Rodin::Kokkos::Variational::GridFunctionCoefficients<Coefficients>::getView(
        fes, coefficient.getData());
    EXPECT_EQ(native.data(), borrowed.data());
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem pb(u, v);
    static_assert(std::is_same_v<typename decltype(pb)::OperatorType, Matrix>);
    static_assert(std::is_same_v<typename decltype(pb)::VectorType, Coefficients>);
    pb = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(coefficient, v);
    pb.assemble();
    auto& system = pb.getLinearSystem();
    static_assert(::KokkosSparse::is_crs_matrix<typename Matrix::View>::value);
    EXPECT_EQ(system.getOperator().getView().values.data(),
      system.getOperator().getValues().data());
    const auto* fieldStorage = u.getSolution().getData().data();
    Solver::CG solver(pb);
    solver.setTolerance(1e-12).setMaxIterations(1000).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_LE((u.getSolution().getData().toHost().array() - 1).matrix().norm(), 1e-10);
    const auto* rhs = system.getVector().data();
    const auto* solution = system.getSolution().data();
    ::Kokkos::deep_copy(native, Real(3));
    pb.assemble();
    solver.solve();
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(u.getSolution().getData().data(), fieldStorage);
    EXPECT_EQ(system.getVector().data(), rhs);
    EXPECT_EQ(system.getSolution().data(), solution);
    EXPECT_LE((u.getSolution().getData().toHost().array() - 3).matrix().norm(), 1e-10);
    const auto previous = system.getSolution().toHost();
    typename decltype(pb)::LinearSystemType other;
    other.getOperator() = system.getOperator();
    other.getVector() = system.getVector();
    solver.solve(other);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ((system.getSolution().toHost() - previous).norm(), 0);
    EXPECT_LE((other.getSolution().toHost() - previous).norm(), 1e-10);
    Matrix independent(system.getOperator());
    EXPECT_NE(independent.getValues().data(), system.getOperator().getValues().data());
    independent *= Real(2);
    EXPECT_LE((independent.toHost() - 2 * system.getOperator().toHost()).norm(), 1e-12);
  }

  TEST(Kokkos_Math, PreassembledFormsAndTargetedAssembly)
  {
    const auto source = makeTriangleMesh();
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    BilinearForm a(u, v);
    a = Integral(Grad(u), Grad(v)) + Integral(u, v);
    a.assemble();
    LinearForm b(v);
    b = Integral(RealFunction(2), v);
    b.assemble();
    Problem pb(u, v);
    using Body = ProblemBody<Matrix, Coefficients, Real>;
    pb = Body(a) - b;
    pb.assemble();
    auto& system = pb.getLinearSystem();
    expectSameOperator(a.getOperator(), system.getOperator());
    expectSameVector(b.getVector(), system.getVector());
    const auto rhs = system.getVector().toHost();
    pb.assemble(AssemblyTarget::LHS);
    EXPECT_EQ((system.getVector().toHost() - rhs).norm(), 0);
    Solver::GMRES solver(pb);
    solver.setTolerance(1e-12).setMaxIterations(1000).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_LE((u.getSolution().getData().toHost().array() - 2).matrix().norm(), 1e-9);
  }
}
