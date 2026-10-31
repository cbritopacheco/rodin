/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstring>

#include <petsc.h>

#include <Rodin/Kokkos/PETSc.h>
#include <Eigen/SparseLU>

#include "../Variational/Common.h"

namespace Rodin::Tests::Unit
{
  using namespace Rodin::Variational;

  /**
   * @brief Copies a PETSc matrix into a sparse matrix.
   * @param A Assembled PETSc matrix.
   * @return Sparse copy of the matrix.
   */
  Math::SparseMatrix<Real> toSparse(::Mat A)
  {
    PetscInt m = 0, n = 0;
    EXPECT_EQ(MatGetSize(A, &m, &n), PETSC_SUCCESS);
    std::vector<Math::SparseTriplet<Real>> triplets;
    for (PetscInt i = 0; i < m; ++i)
    {
      PetscInt count = 0;
      const PetscInt* cols = nullptr;
      const PetscScalar* values = nullptr;
      EXPECT_EQ(MatGetRow(A, i, &count, &cols, &values), PETSC_SUCCESS);
      for (PetscInt k = 0; k < count; ++k)
        triplets.emplace_back(i, cols[k], values[k]);
      EXPECT_EQ(MatRestoreRow(A, i, &count, &cols, &values), PETSC_SUCCESS);
    }
    Math::SparseMatrix<Real> res(m, n);
    res.setFromTriplets(triplets.begin(), triplets.end());
    return res;
  }

  /**
   * @brief Copies a PETSc vector into a dense vector.
   * @param b PETSc vector.
   * @return Dense copy of the vector.
   */
  Math::Vector<Real> toVector(::Vec b)
  {
    PetscInt n = 0;
    EXPECT_EQ(VecGetSize(b, &n), PETSC_SUCCESS);
    const PetscScalar* values = nullptr;
    EXPECT_EQ(VecGetArrayRead(b, &values), PETSC_SUCCESS);
    Math::Vector<Real> res(n);
    for (PetscInt i = 0; i < n; ++i)
      res(i) = values[i];
    EXPECT_EQ(VecRestoreArrayRead(b, &values), PETSC_SUCCESS);
    return res;
  }

  /**
   * @brief Tests whether a PETSc object has a given type.
   * @param object PETSc object.
   * @param type Expected type name.
   */
  template <class Object>
  bool hasType(Object object, const char* type)
  {
    PetscBool same = PETSC_FALSE;
    EXPECT_EQ(PetscObjectTypeCompare(reinterpret_cast<PetscObject>(object), type, &same),
      PETSC_SUCCESS);
    return same == PETSC_TRUE;
  }

  /**
   * @brief Compares a PETSc device system with the sparse device system.
   * @param reference Problem assembled into sparse matrices.
   * @param candidate Problem assembled into PETSc objects.
   */
  template <class Reference, class Candidate>
  void expectSameSystem(Reference& reference, Candidate& candidate)
  {
    reference.assemble();
    candidate.assemble();
    const auto& r = reference.getLinearSystem();
    auto& c = candidate.getLinearSystem();
    EXPECT_TRUE(hasType(c.getOperator(), MATSEQAIJKOKKOS));
    EXPECT_TRUE(hasType(c.getVector(), VECSEQKOKKOS));
    expectSameOperator(r.getOperator(), toSparse(c.getOperator()));
    expectSameVector(r.getVector(), toVector(c.getVector()));
  }

  /**
   * @brief Tags the faces of a mesh with attributes 1, 2 and 3.
   * @param mesh Host mesh.
   * @return Mesh with its face to cell incidence and tagged faces.
   */
  HostMesh withFaces(HostMesh mesh)
  {
    const size_t d = mesh.getDimension();
    mesh.getConnectivity().compute(d - 1, d);
    for (Index i = 0; i < mesh.getFaceCount(); ++i)
      mesh.setAttribute({d - 1, i}, Geometry::Attribute(1 + i % 3));
    return mesh;
  }

  /// @brief Poisson problems match the sparse device assembly.
  void checkPoisson(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);

    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction g(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });

    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv) - Integral(f, sv) +
      DirichletBC(su, g).on(1) + DirichletBC(su, RealFunction(3)).on(2);
    Problem pp(pu, pv);
    pp = Integral(Grad(pu), Grad(pv)) + Integral(pu, pv) - Integral(f, pv) +
      DirichletBC(pu, g).on(1) + DirichletBC(pu, RealFunction(3)).on(2);
    expectSameSystem(sp, pp);
  }

  /// @brief Mixed problems match the sparse device assembly.
  void checkMixedProblem(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 vh(mesh, 2);
    P1 qh(mesh);
    TrialFunction su(vh);
    TrialFunction sp(qh);
    TestFunction sv(vh);
    TestFunction sq(qh);
    PETSc::Variational::TrialFunction pu(vh);
    PETSc::Variational::TrialFunction pp(qh);
    PETSc::Variational::TestFunction pv(vh);
    PETSc::Variational::TestFunction pq(qh);

    VectorFunction f(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() * p.y(); },
      Real(1));
    RealFunction s(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() - p.y(); });
    const Real eps = 1e-3;

    Problem sstokes(su, sp, sv, sq);
    sstokes = Integral(Jacobian(su), Jacobian(sv)) - Integral(sp, Div(sv)) +
      Integral(Div(su), sq) + Integral(eps * sp, sq) - Integral(f, sv) - Integral(s, sq) +
      DirichletBC(su, VectorFunction{0, 0}).on(1) +
      DirichletBC(sp, RealFunction(1)).on(2);
    Problem pstokes(pu, pp, pv, pq);
    pstokes = Integral(Jacobian(pu), Jacobian(pv)) - Integral(pp, Div(pv)) +
      Integral(Div(pu), pq) + Integral(eps * pp, pq) - Integral(f, pv) - Integral(s, pq) +
      DirichletBC(pu, VectorFunction{0, 0}).on(1) +
      DirichletBC(pp, RealFunction(1)).on(2);
    expectSameSystem(sstokes, pstokes);
  }

  /// @brief Boundary terms and linear elasticity match the sparse device assembly.
  void checkBoundaryTerms(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);
    RealFunction g(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });

    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + BoundaryIntegral(su, sv).over(2) -
      BoundaryIntegral(g, sv).over(3) + DirichletBC(su, g).on(1);
    Problem pp(pu, pv);
    pp = Integral(Grad(pu), Grad(pv)) + BoundaryIntegral(pu, pv).over(2) -
      BoundaryIntegral(g, pv).over(3) + DirichletBC(pu, g).on(1);
    expectSameSystem(sp, pp);

    P1 vfes(mesh, source.getSpaceDimension());
    TrialFunction sw(vfes);
    TestFunction sz(vfes);
    PETSc::Variational::TrialFunction pw(vfes);
    PETSc::Variational::TestFunction pz(vfes);
    Problem se(sw, sz);
    se = LinearElasticityIntegral(sw, sz)(Real(1.5), Real(0.7)) -
      BoundaryIntegral(Dot(BoundaryNormal(mesh), sz)).over(3) +
      DirichletBC(sw, VectorFunction{0, 0}).on(1);
    Problem pe(pw, pz);
    pe = LinearElasticityIntegral(pw, pz)(Real(1.5), Real(0.7)) -
      BoundaryIntegral(Dot(BoundaryNormal(mesh), pz)).over(3) +
      DirichletBC(pw, VectorFunction{0, 0}).on(1);
    expectSameSystem(se, pe);
  }

  /// @brief Assembling again keeps the pattern and the values.
  void checkReassembly(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    PETSc::Variational::TrialFunction u(fes);
    PETSc::Variational::TestFunction v(fes);
    RealFunction f(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() + p.y(); });
    Problem pb(u, v);
    pb = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, RealFunction(1));
    pb.assemble();
    auto& axb = pb.getLinearSystem();
    const auto A = toSparse(axb.getOperator());
    const auto b = toVector(axb.getVector());
    PetscObjectState state = 0;
    ASSERT_EQ(MatGetNonzeroState(axb.getOperator(), &state), PETSC_SUCCESS);

    pb.assemble();
    PetscObjectState again = 0;
    ASSERT_EQ(MatGetNonzeroState(axb.getOperator(), &again), PETSC_SUCCESS);
    EXPECT_EQ(state, again);
    expectSameOperator(A, toSparse(axb.getOperator()));
    expectSameVector(b, toVector(axb.getVector()));

    ASSERT_EQ(VecZeroEntries(axb.getVector()), PETSC_SUCCESS);
    pb.assemble(AssemblyTarget::RHS);
    expectSameVector(b, toVector(axb.getVector()));
    ASSERT_EQ(MatZeroEntries(axb.getOperator()), PETSC_SUCCESS);
    pb.assemble(AssemblyTarget::LHS);
    expectSameOperator(A, toSparse(axb.getOperator()));
  }

  /// @brief A Krylov solve on the device matches a sparse direct solve.
  void checkSolve(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });

    Problem sp(su, sv);
    sp =
      Integral(Grad(su), Grad(sv)) - Integral(f, sv) + DirichletBC(su, RealFunction(0));
    sp.assemble();
    const auto& system = sp.getLinearSystem();
    Eigen::SparseLU<Math::SparseMatrix<Real>> lu(system.getOperator());
    const Math::Vector<Real> expected = lu.solve(system.getVector());

    Problem pp(pu, pv);
    pp =
      Integral(Grad(pu), Grad(pv)) - Integral(f, pv) + DirichletBC(pu, RealFunction(0));
    pp.assemble();
    PETSc::Solver::KSP solver(pp);
    solver.setType(KSPCG);
    ::PC pc = nullptr;
    ASSERT_EQ(KSPGetPC(solver.getHandle(), &pc), PETSC_SUCCESS);
    ASSERT_EQ(PCSetType(pc, PCJACOBI), PETSC_SUCCESS);
    ASSERT_EQ(KSPSetTolerances(solver.getHandle(), 1e-13, 1e-50, PETSC_DEFAULT, 1000),
      PETSC_SUCCESS);
    solver.solve();
    ::KSPConvergedReason reason;
    ASSERT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
    EXPECT_GT(reason, 0);
    const auto solution = toVector(pu.getSolution().getData());
    EXPECT_LE((solution - expected).norm(), Real(1e-10) * expected.norm());
  }

  TEST(Kokkos_PETSc_Problem, PoissonOnTriangles)
  {
    checkPoisson(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, PoissonOnMixedMesh)
  {
    checkPoisson(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, PoissonOnTetrahedra)
  {
    checkPoisson(withFaces(makeTetrahedronMesh()));
  }

  TEST(Kokkos_PETSc_Problem, MixedOnTriangles)
  {
    checkMixedProblem(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, MixedOnMixedMesh)
  {
    checkMixedProblem(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, BoundaryTermsOnTriangles)
  {
    checkBoundaryTerms(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, BoundaryTermsOnMixedMesh)
  {
    checkBoundaryTerms(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, Reassembly)
  {
    checkReassembly(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, SolveOnTriangles)
  {
    checkSolve(withFaces(makeTriangleMesh()));
  }
}

/**
 * @brief Runs tests within the lifetimes of the Kokkos and PETSc runtimes.
 * @param argc Number of command-line arguments.
 * @param argv Command-line argument array.
 * @returns GoogleTest exit status.
 */
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  ::Kokkos::ScopeGuard runtime(argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);
  const int res = RUN_ALL_TESTS();
  PetscFinalize();
  return res;
}
