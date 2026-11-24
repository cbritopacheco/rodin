/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Kokkos/Solver.h>
#include <Rodin/Kokkos/Solver/CuDSSOperator.h>
#include <Rodin/Kokkos/Solver/CuDSSVector.h>
#include <Rodin/Kokkos/Solver/CompressedPattern.h>
#include <Rodin/Kokkos/Solver/RocSolverRF.h>
#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;
  using Coefficients = Rodin::Kokkos::Math::Vector<Real, TestBackend>;
  using Matrix = Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend>;

  template <class Destination, class Source>
  concept CanUpload = requires(
    Destination& destination, const Source& source) { destination.upload(source); };

  static_assert(CanUpload<Coefficients, Math::Vector<Real>>);
  static_assert(!CanUpload<const Coefficients, Math::Vector<Real>>);
  static_assert(!CanUpload<Coefficients, Coefficients>);
  static_assert(CanUpload<Matrix, Math::SparseMatrix<Real>>);
  static_assert(!CanUpload<const Matrix, Math::SparseMatrix<Real>>);
  static_assert(!CanUpload<Matrix, Matrix>);
  static_assert(!std::is_convertible_v<Math::Vector<Real>, Coefficients>);
  static_assert(!std::is_assignable_v<Coefficients&, Math::Vector<Real>>);
  static_assert(!std::is_convertible_v<Coefficients, Math::Vector<Real>>);

  /// @brief Symbolic constructors preserve transposition maps and ordering-induced fill.
  TEST(Kokkos_Math, CompressedPatternConstructors)
  {
    using Pattern = Rodin::Kokkos::Solver::CompressedPattern;
    using Indices = std::vector<Math::SparseIndex>;
    const Indices starts{0, 2, 5, 7};
    const Indices indices{0, 1, 0, 1, 2, 1, 2};
    Indices positions;
    const auto transposed =
      Pattern::Transpose(3, starts.data(), indices.data(), positions);
    EXPECT_EQ(transposed.starts, starts);
    EXPECT_EQ(transposed.indices, indices);
    EXPECT_EQ(positions, (Indices{0, 2, 1, 3, 5, 4, 6}));
    const Indices permutation{1, 0, 2};
    const auto cholesky =
      Pattern::Cholesky(3, starts.data(), indices.data(), permutation);
    EXPECT_EQ(cholesky.starts, (Indices{0, 1, 3, 6}));
    EXPECT_EQ(cholesky.indices, (Indices{0, 0, 1, 0, 1, 2}));
    const auto lu =
      Pattern::LU(3, starts.data(), indices.data(), permutation, permutation);
    EXPECT_EQ(lu.starts, (Indices{0, 3, 6, 9}));
    EXPECT_EQ(lu.indices, (Indices{0, 1, 2, 0, 1, 2, 0, 1, 2}));
    const Indices emptyStarts{0};
    const Indices empty;
    const auto emptyTranspose =
      Pattern::Transpose(0, emptyStarts.data(), empty.data(), positions);
    EXPECT_EQ(emptyTranspose.starts, emptyStarts);
    EXPECT_TRUE(emptyTranspose.indices.empty());
    EXPECT_TRUE(positions.empty());
    EXPECT_EQ(
      Pattern::Cholesky(0, emptyStarts.data(), empty.data(), empty).starts, emptyStarts);
    EXPECT_EQ(
      Pattern::LU(0, emptyStarts.data(), empty.data(), empty, empty).starts, emptyStarts);
  }

  /// @brief Transfers replace values and preserve independent host/backend ownership.
  TEST(Kokkos_Math, VectorUploadDownloadSnapshots)
  {
    Math::Vector<Real> source(4);
    source << 1, -2, 3, 4;
    const auto expected = source;
    Coefficients backend;
    EXPECT_EQ(&backend.upload(source), &backend);
    source.setConstant(9);
    const auto& read = backend;
    auto snapshot = read.download();
    EXPECT_EQ((snapshot - expected).norm(), 0);
    Math::Vector<Real> destination(11);
    read.download(destination);
    EXPECT_EQ(destination.size(), expected.size());
    EXPECT_EQ((destination - expected).norm(), 0);
    snapshot.setZero();
    destination.setZero();
    EXPECT_EQ((read.download() - expected).norm(), 0);
    ::Kokkos::deep_copy(backend.getView(), Real(7));
    EXPECT_EQ(source(0), 9);
    EXPECT_EQ(snapshot.norm(), 0);
    EXPECT_EQ((read.download().array() - Real(7)).matrix().norm(), 0);
    backend.upload(expected.head(2).array() * Real(2));
    EXPECT_EQ(backend.size(), 2u);
    EXPECT_EQ((backend.download() - Real(2) * expected.head(2)).norm(), 0);
    Math::Vector<Real> empty;
    backend.upload(empty);
    read.download(destination);
    EXPECT_EQ(destination.size(), 0);
  }

  /// @brief Sparse transfers preserve rectangular operators through CSC/CSR conversion.
  TEST(Kokkos_Math, SparseUploadDownloadSnapshots)
  {
    Math::SparseMatrix<Real> source(3, 5);
    source.insert(2, 0) = 4;
    source.insert(0, 4) = -2;
    source.insert(1, 1) = 3;
    source.insert(2, 3) = 0;
    const auto expected = source;
    Matrix backend;
    EXPECT_EQ(&backend.upload(source), &backend);
    source.coeffRef(2, 0) = 11;
    const auto& read = backend;
    auto snapshot = read.download();
    EXPECT_EQ(snapshot.rows(), 3);
    EXPECT_EQ(snapshot.cols(), 5);
    EXPECT_EQ(snapshot.nonZeros(), expected.nonZeros());
    EXPECT_EQ((snapshot - expected).norm(), 0);
    Math::SparseMatrix<Real> destination(1, 1);
    destination.insert(0, 0) = 99;
    read.download(destination);
    EXPECT_EQ((destination - expected).norm(), 0);
    destination.setZero();
    snapshot.setZero();
    EXPECT_EQ((read.download() - expected).norm(), 0);
    ::Kokkos::deep_copy(backend.getValues(), Real(6));
    EXPECT_EQ(source.coeff(2, 0), 11);
    EXPECT_EQ(snapshot.norm(), 0);
    backend.upload(Math::SparseMatrix<Real>(0, 5));
    read.download(destination);
    EXPECT_EQ(destination.rows(), 0);
    EXPECT_EQ(destination.cols(), 5);
    EXPECT_EQ(destination.nonZeros(), 0);
    Matrix empty;
    empty.download(destination);
    EXPECT_EQ(destination.rows(), 0);
    EXPECT_EQ(destination.cols(), 0);
  }

  /// @brief Vendor input borrows preserve constness and use the original arrays.
  TEST(Kokkos_Math, DirectSolverArrayBorrows)
  {
    using Arrays = Solver::CuDSSOperator<Matrix>;
    using HostMatrix = Math::SparseMatrix<Real>;
    using HostArrays = Solver::CuDSSOperator<HostMatrix>;
    static_assert(std::is_same_v<decltype(Arrays::getStarts(std::declval<Matrix&>())),
      Math::SparseIndex*>);
    static_assert(std::is_same_v<decltype(Arrays::getIndices(std::declval<Matrix&>())),
      Math::SparseIndex*>);
    static_assert(
      std::is_same_v<decltype(Arrays::getValues(std::declval<Matrix&>())), Real*>);
    static_assert(
      std::is_same_v<decltype(Arrays::getStarts(std::declval<const Matrix&>())),
        const Math::SparseIndex*>);
    static_assert(
      std::is_same_v<decltype(Arrays::getIndices(std::declval<const Matrix&>())),
        const Math::SparseIndex*>);
    static_assert(
      std::is_same_v<decltype(Arrays::getValues(std::declval<const Matrix&>())),
        const Real*>);
    static_assert(
      std::is_same_v<decltype(HostArrays::getStarts(std::declval<HostMatrix&>())),
        Math::SparseIndex*>);
    static_assert(
      std::is_same_v<decltype(HostArrays::getIndices(std::declval<HostMatrix&>())),
        Math::SparseIndex*>);
    static_assert(
      std::is_same_v<decltype(HostArrays::getValues(std::declval<HostMatrix&>())),
        Real*>);
    static_assert(
      std::is_same_v<decltype(HostArrays::getValues(std::declval<const HostMatrix&>())),
        const Real*>);
    static_assert(
      std::is_same_v<decltype(HostArrays::getStarts(std::declval<const HostMatrix&>())),
        const Math::SparseIndex*>);
    static_assert(
      std::is_same_v<decltype(HostArrays::getIndices(std::declval<const HostMatrix&>())),
        const Math::SparseIndex*>);
    HostMatrix source(3, 3);
    source.insert(0, 0) = 2;
    source.insert(1, 1) = 3;
    source.insert(2, 2) = 4;
    HostArrays::prepare(source);
    EXPECT_EQ(HostArrays::getStarts(source), source.outerIndexPtr());
    EXPECT_EQ(HostArrays::getIndices(source), source.innerIndexPtr());
    EXPECT_EQ(HostArrays::getValues(source), source.valuePtr());
    const auto make = [&] {
      const auto execution = TestContext{}.getExecutionSpace();
      return Matrix(3, 3,
        Rodin::Kokkos::Internal::upload(
          execution, std::vector<Math::SparseIndex>{0, 1, 2, 3}, "borrow offsets"),
        Rodin::Kokkos::Internal::upload(
          execution, std::vector<Math::SparseIndex>{0, 1, 2}, "borrow columns"),
        Rodin::Kokkos::Internal::upload(
          execution, std::vector<Real>{2, 3, 4}, "borrow values"));
    };
    Matrix matrix = make();
    const Matrix& read = matrix;
    const auto check = [&] {
      EXPECT_EQ(Arrays::getStarts(matrix), matrix.getRowStarts().data());
      EXPECT_EQ(Arrays::getIndices(matrix), matrix.getColumns().data());
      EXPECT_EQ(Arrays::getValues(matrix), matrix.getValues().data());
      EXPECT_EQ(Arrays::getStarts(read), Arrays::getStarts(matrix));
      EXPECT_EQ(Arrays::getIndices(read), Arrays::getIndices(matrix));
      EXPECT_EQ(Arrays::getValues(read), Arrays::getValues(matrix));
    };
    check();
    const auto* starts = Arrays::getStarts(matrix);
    const auto* indices = Arrays::getIndices(matrix);
    const auto* values = Arrays::getValues(matrix);
    const auto borrowed =
      ::Kokkos::View<Real*, TestMemorySpace, ::Kokkos::MemoryTraits<::Kokkos::Unmanaged>>(
        Arrays::getValues(matrix), 3);
    ::Kokkos::deep_copy(borrowed, Real(7));
    EXPECT_EQ(matrix.download().coeff(1, 1), Real(7));
    EXPECT_EQ(Arrays::getStarts(matrix), starts);
    EXPECT_EQ(Arrays::getIndices(matrix), indices);
    EXPECT_EQ(Arrays::getValues(matrix), values);
    Matrix replacement = make();
    matrix = std::move(replacement);
    check();
    EXPECT_NE(Arrays::getStarts(matrix), starts);
    EXPECT_NE(Arrays::getIndices(matrix), indices);
    EXPECT_NE(Arrays::getValues(matrix), values);
  }

#if defined(RODIN_USE_CUDSS) || defined(RODIN_USE_ROCSOLVER)
  /**
   * @brief Repeated device factorizations borrow input arrays without modifying them.
   * @tparam Direct Vendor direct solver template.
   */
  template <template <class> class Direct>
  void checkDirectSolverInputBorrows()
  {
    const auto source = makeTriangleMesh();
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem pb(u, v);
    const auto diffusion = Integral(Grad(u), Grad(v));
    pb = diffusion + Integral(u, v) - Integral(RealFunction(1), v);
    pb.assemble();
    auto& system = pb.getLinearSystem();
    auto& matrix = system.getOperator();
    using System = std::decay_t<decltype(system)>;
    Direct<System> solver(pb);
    static_assert(Solver::CuDSSOperator<Matrix>::Resident);
    const auto check = [&](Real expected) {
      const auto before = matrix.download();
      typename Solver::CuDSSOperator<Matrix>::Pattern pattern;
      pattern.assign(matrix);
      const auto* starts = matrix.getRowStarts().data();
      const auto* columns = matrix.getColumns().data();
      const auto* values = matrix.getValues().data();
      solver.factorize(system);
      ASSERT_TRUE(solver.success());
      solver.solve(system);
      ASSERT_TRUE(solver.success());
      EXPECT_LT(
        (system.getSolution().download().array() - expected).matrix().norm(), 1e-10);
      EXPECT_EQ((matrix.download() - before).norm(), Real(0));
      EXPECT_TRUE(pattern.matches(matrix));
      EXPECT_EQ(matrix.getRowStarts().data(), starts);
      EXPECT_EQ(matrix.getColumns().data(), columns);
      EXPECT_EQ(matrix.getValues().data(), values);
      // Native input arrays are borrowed, not duplicated in the host-input buffers.
      EXPECT_EQ(solver.getResources().outer.size, 0u);
      EXPECT_EQ(solver.getResources().inner.size, 0u);
      EXPECT_EQ(solver.getResources().values.size, 0u);
    };
    check(Real(1));
    const auto values = matrix.getValues();
    ::Kokkos::parallel_for(
      "scale vendor input",
      ::Kokkos::RangePolicy<TestExecutionSpace>(0, values.extent(0)),
      KOKKOS_LAMBDA(size_t i) { values(i) *= Real(2); });
    check(Real(0.5));
    const auto* previous = matrix.getValues().data();
    matrix = Matrix(matrix);
    EXPECT_NE(matrix.getValues().data(), previous);
    check(Real(0.5));
  }
#endif

#ifdef RODIN_USE_CUDSS
  /// @brief cuDSS keeps native inputs in place through reuse and reallocation.
  TEST(Kokkos_Math, CuDSSInputArraysAreBorrowed)
  {
    checkDirectSolverInputBorrows<Solver::CuDSS>();
  }
#endif

#ifdef RODIN_USE_ROCSOLVER
  /// @brief rocSOLVER keeps native inputs in place through reuse and reallocation.
  TEST(Kokkos_Math, RocSolverInputArraysAreBorrowed)
  {
    checkDirectSolverInputBorrows<Solver::RocSolverRF>();
  }
#endif

  /// @brief Pivot certification must not reject small or large finite values.
  TEST(Kokkos_Math, RocSolverPivotPredicate)
  {
    using Pivots = Rodin::Kokkos::Solver::RocSolverRFPivots;
    for (const Real value : {Real(1), Real(-1), std::numeric_limits<Real>::min(),
           -std::numeric_limits<Real>::min(), std::numeric_limits<Real>::denorm_min(),
           -std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::max(),
           -std::numeric_limits<Real>::max()})
    {
      EXPECT_TRUE(Pivots::valid(value));
    }
    for (const Real value :
      {Real(0), -Real(0), std::numeric_limits<Real>::infinity(),
        -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
    {
      EXPECT_FALSE(Pivots::valid(value));
    }
  }

  /// @brief Exercises the actual pivot reduction on the selected backend.
  TEST(Kokkos_Math, RocSolverPivotReduction)
  {
    using Pivots = Rodin::Kokkos::Solver::RocSolverRFPivots;
    constexpr size_t Count = 4097;
    Matrix::Array<Real> values("test pivots", Count);
    const auto execution = TestContext{}.getExecutionSpace();
    const auto host = ::Kokkos::create_mirror_view(values);
    for (size_t i = 0; i < Count; ++i)
      host(i) = Real(1);
    host(1) = std::numeric_limits<Real>::max();
    host(2) = -std::numeric_limits<Real>::max();
    host(3) = std::numeric_limits<Real>::denorm_min();
    host(4) = -std::numeric_limits<Real>::denorm_min();
    const auto run = [&](bool expected) {
      ::Kokkos::deep_copy(execution, values, host);
      EXPECT_EQ(Pivots::check(execution, values), expected);
    };
    run(true);
    for (const size_t position : {size_t(0), Count / 2, Count - 1})
    {
      for (const Real invalid :
        {Real(0), -Real(0), std::numeric_limits<Real>::infinity(),
          -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()})
      {
        host(position) = invalid;
        run(false);
        host(position) = Real(1);
        run(true);
      }
    }
    EXPECT_TRUE(Pivots::check(execution, Matrix::Array<Real>()));
  }

  /// @brief A solver snapshot detects edits even when sizes and pointers stay fixed.
  void checkCuDSSPatternExactComparison()
  {
    using Pattern = Solver::CuDSSOperator<Matrix>::Pattern;
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
    using Arrays = Solver::CuDSSOperator<HostMatrix>;
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
    Coefficients a;
    a.upload(source);
    Coefficients b(a);
    EXPECT_NE(a.getView().data(), b.getView().data());
    ::Kokkos::deep_copy(b.getView(), Real(7));
    EXPECT_EQ(a.download()(0), 1);
    EXPECT_EQ(b.download()(0), 7);
    b.segment(1, 2) = static_cast<const Coefficients&>(a).segment(2, 2);
    EXPECT_EQ(b.download()(1), 3);
    EXPECT_EQ(b.download()(2), 4);
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

  TEST(Kokkos_Math, GridFunctionCopiesCrossBackendSlices)
  {
    const auto source = makeTriangleMesh();
    const DeviceMesh mesh(TestContext{}, source);
    P1 hostSpace(source);
    P1 backendSpace(mesh);
    GridFunction host(hostSpace);
    GridFunction backend(backendSpace);
    const auto size = backendSpace.getSize();
    const Math::Vector<Real> values =
      Math::Vector<Real>::LinSpaced(size + 2, 1, size + 2);
    backend.setData(values, 1);
    EXPECT_EQ((backend.getData().download() - values.segment(1, size)).norm(), Real(0));

    Coefficients uploaded;
    uploaded.upload(values);
    host.setData(uploaded, 2);
    EXPECT_EQ((host.getData() - values.segment(2, size)).norm(), Real(0));
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
      Rodin::Kokkos::Variational::GridFunctionData<Coefficients>::prepare(
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
    EXPECT_LE((u.getSolution().getData().download().array() - 1).matrix().norm(), 1e-10);
    const auto* rhs = system.getVector().data();
    const auto* solution = system.getSolution().data();
    ::Kokkos::deep_copy(native, Real(3));
    pb.assemble();
    solver.solve();
    ASSERT_TRUE(solver.success());
    EXPECT_EQ(u.getSolution().getData().data(), fieldStorage);
    EXPECT_EQ(system.getVector().data(), rhs);
    EXPECT_EQ(system.getSolution().data(), solution);
    EXPECT_LE((u.getSolution().getData().download().array() - 3).matrix().norm(), 1e-10);
    const auto previous = system.getSolution().download();
    typename decltype(pb)::LinearSystemType other;
    other.getOperator() = system.getOperator();
    other.getVector() = system.getVector();
    solver.solve(other);
    ASSERT_TRUE(solver.success());
    EXPECT_EQ((system.getSolution().download() - previous).norm(), 0);
    EXPECT_LE((other.getSolution().download() - previous).norm(), 1e-10);
    Matrix independent(system.getOperator());
    EXPECT_NE(independent.getValues().data(), system.getOperator().getValues().data());
    independent *= Real(2);
    EXPECT_LE(
      (independent.download() - 2 * system.getOperator().download()).norm(), 1e-12);
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
    const auto rhs = system.getVector().download();
    pb.assemble(AssemblyTarget::LHS);
    EXPECT_EQ((system.getVector().download() - rhs).norm(), 0);
    Solver::GMRES solver(pb);
    solver.setTolerance(1e-12).setMaxIterations(1000).solve();
    ASSERT_TRUE(solver.success());
    EXPECT_LE((u.getSolution().getData().download().array() - 2).matrix().norm(), 1e-9);
  }
}
