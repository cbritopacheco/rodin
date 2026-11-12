/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Kokkos/Solver.h>
#include <Rodin/Solver/CuDSSOperator.h>
#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;
  using Coefficients = Rodin::Kokkos::Math::Vector<Real, TestBackend>;
  using Matrix = Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>;

  /// @brief A solver snapshot detects edits even when sizes and pointers stay fixed.
  void checkCuDSSPatternExactComparison()
  {
    using Pattern = Solver::Internal::CuDSSOperator<Matrix>::Pattern;
    const auto execution = TestContext{}.getExecutionSpace();
    const auto make = [&](size_t rows, size_t columns,
                        const std::vector<Math::SparseIndex>& offsets,
                        const std::vector<Math::SparseIndex>& indices) {
      auto starts = Rodin::Kokkos::Internal::upload(execution, offsets, "test offsets");
      auto entries = Rodin::Kokkos::Internal::upload(execution, indices, "test columns");
      Matrix::Array<Real> values("test values", indices.size());
      ::Kokkos::deep_copy(values, Real(2));
      return Matrix(rows, columns, starts, entries, values);
    };
    auto matrix = make(3, 3, {0, 1, 2, 3}, {0, 1, 2});
    Pattern pattern;
    EXPECT_FALSE(pattern.matches(matrix));
    pattern.assign(matrix);
    EXPECT_TRUE(pattern.matches(matrix));
    EXPECT_TRUE(pattern.matches(matrix));

    // Numeric zeros do not change a stored CSR graph.
    ::Kokkos::deep_copy(matrix.getValues(), Real(0));
    EXPECT_TRUE(pattern.matches(matrix));
    Matrix relocated(matrix);
    EXPECT_NE(matrix.getColumns().data(), relocated.getColumns().data());
    EXPECT_TRUE(pattern.matches(relocated));

    const auto columns = matrix.getColumns();
    ::Kokkos::parallel_for(
      "edit column", ::Kokkos::RangePolicy<Matrix::ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { columns(0) = 1; });
    EXPECT_FALSE(pattern.matches(matrix));
    EXPECT_TRUE(pattern.matches(relocated));
    pattern.assign(matrix);
    EXPECT_TRUE(pattern.matches(matrix));
    EXPECT_FALSE(pattern.matches(relocated));

    // Identical column arrays and counts can still have different row boundaries.
    pattern.assign(relocated);
    const auto starts = relocated.getRowStarts();
    ::Kokkos::parallel_for(
      "edit offset", ::Kokkos::RangePolicy<Matrix::ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { starts(1) = 2; });
    EXPECT_FALSE(pattern.matches(relocated));
    pattern.assign(relocated);
    EXPECT_TRUE(pattern.matches(relocated));

    EXPECT_FALSE(pattern.matches(make(3, 3, {0, 1, 1, 2}, {0, 2})));
    EXPECT_FALSE(pattern.matches(make(2, 3, {0, 1, 3}, {0, 1, 2})));
    EXPECT_FALSE(pattern.matches(make(3, 4, {0, 2, 2, 3}, {0, 1, 2})));
    pattern.clear();
    EXPECT_FALSE(pattern.matches(relocated));

    auto emptyRows = make(3, 3, {0, 0, 0, 0}, {});
    pattern.assign(emptyRows);
    EXPECT_TRUE(pattern.matches(emptyRows));
    EXPECT_FALSE(pattern.matches(matrix));
    auto empty = make(0, 0, {0}, {});
    pattern.assign(empty);
    EXPECT_TRUE(pattern.matches(empty));
    pattern.clear();
    EXPECT_FALSE(pattern.matches(empty));

    // A mismatch must survive reduction across several GPU blocks.
    constexpr size_t Count = 4097;
    Matrix::Array<Math::SparseIndex> manyStarts("many offsets", Count + 1);
    Matrix::Array<Math::SparseIndex> manyColumns("many columns", Count);
    Matrix::Array<Real> manyValues("many values", Count);
    ::Kokkos::parallel_for(
      "initialize diagonal graph",
      ::Kokkos::RangePolicy<Matrix::ExecutionSpace>(0, Count + 1), KOKKOS_LAMBDA(int i) {
        manyStarts(i) = i;
        if (i < static_cast<int>(Count))
          manyColumns(i) = i;
      });
    Matrix large(Count, Count, manyStarts, manyColumns, manyValues);
    pattern.assign(large);
    EXPECT_TRUE(pattern.matches(large));
    for (const size_t position : {size_t(0), Count / 2, Count - 1})
    {
      ::Kokkos::parallel_for(
        "change diagonal graph", ::Kokkos::RangePolicy<Matrix::ExecutionSpace>(0, 1),
        KOKKOS_LAMBDA(int) { manyColumns(position) = (position + 1) % Count; });
      EXPECT_FALSE(pattern.matches(large));
      ::Kokkos::parallel_for(
        "restore diagonal graph", ::Kokkos::RangePolicy<Matrix::ExecutionSpace>(0, 1),
        KOKKOS_LAMBDA(int) { manyColumns(position) = position; });
      EXPECT_TRUE(pattern.matches(large));
    }
  }

  TEST(Kokkos_Math, CuDSSPatternExactComparison)
  {
    checkCuDSSPatternExactComparison();
  }

  /// @brief Host matrices retain the same exact graph semantics without CUDA.
  TEST(Kokkos_Math, CuDSSHostPatternComparison)
  {
    using HostMatrix = Math::SparseMatrix<Real>;
    using Arrays = Solver::Internal::CuDSSOperator<HostMatrix>;
    HostMatrix matrix(3, 3);
    matrix.insert(0, 0) = 2;
    matrix.insert(1, 1) = 3;
    matrix.insert(2, 2) = 4;
    Arrays::prepare(matrix);
    Arrays::Pattern pattern;
    EXPECT_FALSE(pattern.matches(matrix));
    pattern.assign(matrix);
    EXPECT_TRUE(pattern.matches(matrix));
    matrix.valuePtr()[0] = 0;
    EXPECT_TRUE(pattern.matches(matrix));
    HostMatrix relocated(matrix);
    EXPECT_TRUE(pattern.matches(relocated));
    matrix.innerIndexPtr()[0] = 1;
    EXPECT_FALSE(pattern.matches(matrix));
    EXPECT_TRUE(pattern.matches(relocated));
    pattern.assign(matrix);
    EXPECT_TRUE(pattern.matches(matrix));
    EXPECT_FALSE(pattern.matches(relocated));
    pattern.clear();
    EXPECT_FALSE(pattern.matches(matrix));
  }

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
