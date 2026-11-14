/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Variational;
  using SymbolicAssembly = Assembly::Internal::SystemAssembly<TestBackend>;

  /**
   * @brief Compares a retained assembly with an independent, empty-policy copy.
   * @param problem Problem whose symbolic state may already be populated.
   */
  template <class ProblemType>
  void expectFreshAssembly(ProblemType& problem)
  {
    ProblemType fresh(problem);
    fresh.assemble();
    problem.assemble();
    const auto a = problem.getLinearSystem().getOperator().toHost();
    const auto b = fresh.getLinearSystem().getOperator().toHost();
    EXPECT_EQ(a.nonZeros(), b.nonZeros());
    EXPECT_EQ((a - b).norm(), Real(0));
    EXPECT_EQ((problem.getLinearSystem().getVector().toHost() -
                fresh.getLinearSystem().getVector().toHost())
                .norm(),
      Real(0));
  }

  /**
   * @brief Uploads a fixture array without retaining a host mirror.
   * @param values Array entries. @returns Independent backend array.
   */
  template <class T>
  typename SymbolicAssembly::template Array<T> symbolicArray(
    std::initializer_list<T> values)
  {
    typename SymbolicAssembly::template Array<T> result(
      "symbolic fixture", values.size());
    const ::Kokkos::View<const T*, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> host(
      values.begin(), values.size());
    ::Kokkos::deep_copy(result, host);
    return result;
  }

  TEST(Kokkos_Variational_SymbolicAssembly, ZeroCrossingsRetainPotentialCouplings)
  {
    const DeviceMesh mesh(TestContext{}, makeTriangleMesh());
    P1 fes(mesh);
    GridFunction coefficient(fes);
    TrialFunction u(fes);
    TestFunction v(fes);
    auto integral = Integral(coefficient * u, v);
    auto load = Integral(coefficient, v);
    SymbolicAssembly::Symbolic symbolic;
    Rodin::Kokkos::Math::SparseMatrix<Real, TestBackend> matrix;
    Rodin::Kokkos::Math::Vector<Real, TestBackend> rhs;
    const auto assemble = [&] {
      SymbolicAssembly engine(
        mesh.getContext().getExecutionSpace(), fes.getSize(), fes.getSize(), &symbolic);
      engine.addMatrix(integral, 0, 0);
      engine.addVector(load, 0, Real(1));
      engine.finish(matrix, rhs);
    };
    coefficient = RealFunction(0);
    assemble();
    EXPECT_EQ(matrix.getNonZeroCount(), 0);
    EXPECT_EQ(rhs.norm(), Real(0));
    ASSERT_EQ(symbolic.patterns.size(), 2);
    const auto vectorOrder = symbolic.patterns[1].ordered;
    const auto order = symbolic.patterns[0].ordered;
    EXPECT_GT(order.extent(0), 0);
    coefficient = RealFunction(2);
    assemble();
    EXPECT_EQ(symbolic.patterns[0].ordered.data(), order.data());
    BilinearForm expected(u, v);
    expected = Integral(Real(2) * u, v);
    expected.assemble();
    expectSameOperator(matrix, expected.getOperator());
    LinearForm expectedLoad(v);
    expectedLoad = Integral(RealFunction(2), v);
    expectedLoad.assemble();
    expectSameVector(expectedLoad.getVector().toHost(), rhs.toHost());
    EXPECT_EQ(symbolic.patterns[1].ordered.data(), vectorOrder.data());
    coefficient = RealFunction(0);
    assemble();
    EXPECT_EQ(matrix.getNonZeroCount(), 0);
    EXPECT_EQ(symbolic.patterns[0].ordered.data(), order.data());
    EXPECT_EQ(symbolic.patterns[1].ordered.data(), vectorOrder.data());
    EXPECT_EQ(rhs.norm(), Real(0));
    SymbolicAssembly::Symbolic copied(symbolic);
    EXPECT_TRUE(copied.patterns.empty());
    copied = symbolic;
    EXPECT_TRUE(copied.patterns.empty());
  }

  TEST(Kokkos_Variational_SymbolicAssembly, OrderedCoordinatesAndRanksAreValidated)
  {
    SymbolicAssembly::Symbolic symbolic;
    SymbolicAssembly::Ranks ranks;
    const auto values = symbolicArray<Real>({1, 2, 4, 8, 16, 32});
    const auto evaluate = [&](bool swapped, const SymbolicAssembly::Ranks& ordering,
                            SymbolicAssembly::Symbolic* cache) {
      SymbolicAssembly engine(TestContext{}.getExecutionSpace(), 2, 2, cache);
      SymbolicAssembly::Listed listed{swapped ? symbolicArray<Index>({1, 0, 1, 0, 1, 0})
                                              : symbolicArray<Index>({0, 1, 0, 1, 0, 1}),
        symbolicArray<Index>({0, 1, 2, 3, 4, 5})};
      const auto groups = engine.groupContributions(listed, ordering, values, 6);
      EXPECT_EQ(groups.size(), 2);
      return ::Kokkos::create_mirror_view_and_copy(
        ::Kokkos::HostSpace(), listed.positions);
    };
    const auto first = evaluate(false, ranks, &symbolic);
    const auto previous = symbolic.patterns[0].ordered;
    const auto changed = evaluate(true, ranks, &symbolic);
    const auto fresh = evaluate(true, ranks, nullptr);
    EXPECT_NE(previous.data(), symbolic.patterns[0].ordered.data());
    for (size_t i = 0; i < changed.extent(0); ++i)
      EXPECT_EQ(changed(i), fresh(i));
    EXPECT_NE(first(0), changed(0));

    // Equal coordinate lists can still require a different summation order.
    ranks.segments = symbolicArray<SymbolicAssembly::Segment>({{0, 6, 2, 0, false}});
    const auto oldOrder = symbolic.patterns[0].ordered;
    const auto reordered = evaluate(true, ranks, &symbolic);
    const auto reference = evaluate(true, ranks, nullptr);
    EXPECT_NE(oldOrder.data(), symbolic.patterns[0].ordered.data());
    for (size_t i = 0; i < reordered.extent(0); ++i)
      EXPECT_EQ(reordered(i), reference(i));
  }

  TEST(Kokkos_Variational_SymbolicAssembly, MovedAttributesAndReplacedFormulation)
  {
    DeviceMesh mesh(TestContext{}, makeMixedMesh());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem problem(u, v);
    problem = Integral(u, v).over(1) - Integral(RealFunction(1), v);
    expectFreshAssembly(problem);
    const auto previous = problem.getLinearSystem().getOperator().toHost();
    for (Index i = 0; i < mesh.getCellCount(); ++i)
      mesh.setAttribute({2, i}, Geometry::Attribute(2 - i % 2));
    expectFreshAssembly(problem);
    EXPECT_GT(
      (problem.getLinearSystem().getOperator().toHost() - previous).norm(), Real(0));
    problem = Integral(Grad(u), Grad(v)) + Integral(u, v);
    expectFreshAssembly(problem);
    auto copied(problem);
    expectFreshAssembly(copied);
    expectFreshAssembly(problem);
  }

  /**
   * @brief Checks numerical refresh after vertex edits with unchanged topology.
   * @param source Mesh whose opposite corner vertices are displaced outward.
   */
  void checkMovedVertexAssembly(const HostMesh& source)
  {
    DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem reused(u, v);
    reused = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1), v);
    reused.assemble();

    // Small outward displacements change geometry without changing connectivity.
    constexpr Real displacement = Real(0.025);
    for (const Index vertex : {Index(0), Index(mesh.getVertexCount() - 1)})
    {
      const auto previous = reused.getLinearSystem().getOperator().toHost();
      const auto previousRHS = reused.getLinearSystem().getVector().toHost();
      auto point = mesh.getVertexCoordinates(vertex);
      const Real direction = vertex == 0 ? Real(-1) : Real(1);
      for (size_t s = 0; s < mesh.getSpaceDimension(); ++s)
        point(s) += direction * displacement * Real(s + 1);
      mesh.setVertexCoordinates(vertex, point);
      reused.assemble();

      // Construct the reference spaces and integrators after the geometry edit.
      P1 freshFES(mesh);
      TrialFunction freshU(freshFES);
      TestFunction freshV(freshFES);
      Problem fresh(freshU, freshV);
      fresh = Integral(Grad(freshU), Grad(freshV)) + Integral(freshU, freshV) -
        Integral(RealFunction(1), freshV);
      fresh.assemble();

      const auto actual = reused.getLinearSystem().getOperator().toHost();
      const auto expected = fresh.getLinearSystem().getOperator().toHost();
      const auto actualRHS = reused.getLinearSystem().getVector().toHost();
      ASSERT_EQ(actual.rows(), expected.rows());
      ASSERT_EQ(actual.cols(), expected.cols());
      EXPECT_EQ(actual.nonZeros(), expected.nonZeros());
      EXPECT_EQ((actual - expected).norm(), Real(0));
      EXPECT_EQ(
        (actualRHS - fresh.getLinearSystem().getVector().toHost()).norm(), Real(0));
      EXPECT_GT((actual - previous).norm(), Real(0));
      EXPECT_GT((actualRHS - previousRHS).norm(), Real(0));
      ASSERT_EQ(actual.nonZeros(), previous.nonZeros());
      for (Index i = 0; i <= static_cast<Index>(actual.outerSize()); ++i)
        EXPECT_EQ(actual.outerIndexPtr()[i], previous.outerIndexPtr()[i]);
      for (Index i = 0; i < static_cast<Index>(actual.nonZeros()); ++i)
        EXPECT_EQ(actual.innerIndexPtr()[i], previous.innerIndexPtr()[i]);
    }
  }

  TEST(Kokkos_Variational_SymbolicAssembly, MovedVerticesOnTriangles)
  {
    checkMovedVertexAssembly(makeTriangleMesh());
  }

  TEST(Kokkos_Variational_SymbolicAssembly, MovedVerticesOnTetrahedra)
  {
    checkMovedVertexAssembly(makeTetrahedronMesh());
  }

  TEST(Kokkos_Variational_SymbolicAssembly, FixedDOFsAndBoundaryValues)
  {
    DeviceMesh mesh(TestContext{}, withClosure(makeTriangleMesh()));
    std::vector<Index> faces;
    for (Index i = 0; i < mesh.getFaceCount(); ++i)
    {
      if (mesh.isBoundary(i))
        faces.push_back(i);
    }
    ASSERT_GT(faces.size(), 1);
    mesh.setAttribute({1, faces.front()}, 5);
    P1 fes(mesh);
    GridFunction value(fes);
    value = RealFunction(1);
    TrialFunction u(fes);
    TestFunction v(fes);
    Problem problem(u, v);
    problem = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1), v) +
      DirichletBC(u, value).on(5);
    expectFreshAssembly(problem);
    const auto previous = problem.getLinearSystem().getVector().toHost();
    value = RealFunction(3);
    expectFreshAssembly(problem);
    EXPECT_GT(
      (problem.getLinearSystem().getVector().toHost() - previous).norm(), Real(0));
    mesh.setAttribute({1, faces.front()}, 0);
    mesh.setAttribute({1, faces.back()}, 5);
    expectFreshAssembly(problem);

    // Replace the relation while keeping the dimensions and pair count fixed.
    IndexMap<IndexSet> pairs;
    pairs[0] = {1};
    problem = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1), v) +
      PeriodicBC(u, pairs);
    expectFreshAssembly(problem);
    pairs[0] = {2};
    problem = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1), v) +
      PeriodicBC(u, pairs);
    expectFreshAssembly(problem);
  }
}
