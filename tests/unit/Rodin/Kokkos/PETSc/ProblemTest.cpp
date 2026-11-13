/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstring>

#include <petsc.h>

#include <Rodin/Kokkos/PETSc.h>
#include <Rodin/Solver/SparseLU.h>

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

  /// @brief Compares the systems of two assembled problems.
  template <class Reference, class Candidate>
  void expectSameAssembled(Reference& reference, Candidate& candidate)
  {
    const auto& r = reference.getLinearSystem();
    auto& c = candidate.getLinearSystem();
    EXPECT_TRUE(hasType(c.getOperator(), MATSEQAIJKOKKOS));
    EXPECT_TRUE(hasType(c.getVector(), VECSEQKOKKOS));
    expectSameOperator(r.getOperator(), toSparse(c.getOperator()));
    expectSameVector(r.getVector(), toVector(c.getVector()));
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
    expectSameAssembled(reference, candidate);
  }

  /**
   * @brief Assembles a device PETSc problem before its reference, and
   * compares their systems.
   *
   * No upload of the reference comes between two assemblies of the
   * candidate, so that the candidate reads again what it kept, if anything.
   */
  template <class Reference, class Candidate>
  void expectSameSystemAssembledFirst(Reference& reference, Candidate& candidate)
  {
    candidate.assemble();
    reference.assemble();
    expectSameAssembled(reference, candidate);
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

    // The block diagonal terms share the full entries of the divergence term.
    Problem sgrad(su, sv);
    sgrad = Integral(su, sv) + Integral(Jacobian(su), Jacobian(sv)) +
      Integral(Div(su), Div(sv)) - Integral(f, sv) +
      DirichletBC(su, VectorFunction{0, 0}).on(1);
    Problem pgrad(pu, pv);
    pgrad = Integral(pu, pv) + Integral(Jacobian(pu), Jacobian(pv)) +
      Integral(Div(pu), Div(pv)) - Integral(f, pv) +
      DirichletBC(pu, VectorFunction{0, 0}).on(1);
    expectSameSystem(sgrad, pgrad);
  }

  /// @brief Potentials, restricted or in a block, match the sparse device assembly.
  void checkPotential(const HostMesh& source)
  {
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);

    const auto k = KOKKOS_LAMBDA(
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
                     ->Real
    {
      return 1 + x.x() * y.y();
    };
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });

    auto si = Integral(Potential(k, su), sv);
    si.setTrialAttributes({1});
    auto pi = Integral(Potential(k, pu), pv);
    pi.setTrialAttributes({1});
    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + si - Integral(f, sv) +
      DirichletBC(su, RealFunction(3)).on(2);
    Problem pp(pu, pv);
    pp = Integral(Grad(pu), Grad(pv)) + pi - Integral(f, pv) +
      DirichletBC(pu, RealFunction(3)).on(2);
    expectSameSystem(sp, pp);

    P1 vh(mesh, 2);
    TrialFunction sw(vh);
    TestFunction sz(vh);
    PETSc::Variational::TrialFunction pw(vh);
    PETSc::Variational::TestFunction pz(vh);
    const auto m = KOKKOS_LAMBDA(DeviceGeometry::Matrix & a,
      const DeviceGeometry::PointView& x, const DeviceGeometry::PointView& y)
    {
      for (Index i = 0; i < a.rows(); ++i)
      {
        for (Index j = 0; j < a.cols(); ++j)
          a(i, j) = (i == j ? 2 : 0) + (i + 1) * x.x() * y.y() - Real(0.25) * j * y.x();
      }
    };
    Problem sb(sw, su, sz, sv);
    sb = Integral(Jacobian(sw), Jacobian(sz)) - Integral(Potential(m, sw), sz) +
      Integral(Grad(su), Grad(sv)) + Integral(Potential(k, su), sv) - Integral(f, sv) +
      DirichletBC(sw, VectorFunction{0, 0}).on(1);
    Problem pb(pw, pu, pz, pv);
    pb = Integral(Jacobian(pw), Jacobian(pz)) - Integral(Potential(m, pw), pz) +
      Integral(Grad(pu), Grad(pv)) + Integral(Potential(k, pu), pv) - Integral(f, pv) +
      DirichletBC(pw, VectorFunction{0, 0}).on(1);
    expectSameSystem(sb, pb);
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

  /// @brief Periodic conditions and identifications match the sparse device assembly.
  /// @param source Mesh tagged by withTaggedFace().
  void checkConstraints(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);
    const auto pairs = makePairs(fes.getSize(), fixedDOFs(su, 5));
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction g(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 2 - p.x(); });

    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv) - Integral(f, sv) +
      DirichletBC(su, g).on(5) + PeriodicBC(su, pairs);
    Problem pp(pu, pv);
    pp = Integral(Grad(pu), Grad(pv)) + Integral(pu, pv) - Integral(f, pv) +
      DirichletBC(pu, g).on(5) + PeriodicBC(pu, pairs);
    expectSameSystem(sp, pp);

    Problem si(su, sv);
    si = Integral(Grad(su), Grad(sv)) - Integral(f, sv) + DirichletBC(su, -su).on(1) +
      DirichletBC(su, RealFunction(3)).on(2);
    Problem pi(pu, pv);
    pi = Integral(Grad(pu), Grad(pv)) - Integral(f, pv) + DirichletBC(pu, -pu).on(1) +
      DirichletBC(pu, RealFunction(3)).on(2);
    expectSameSystem(si, pi);

    P1 vh(mesh, 2);
    TrialFunction sw(vh);
    TrialFunction sa(fes);
    TestFunction sz(vh);
    TestFunction sb(fes);
    PETSc::Variational::TrialFunction pw(vh);
    PETSc::Variational::TrialFunction pa(fes);
    PETSc::Variational::TestFunction pz(vh);
    PETSc::Variational::TestFunction pb(fes);
    const Real eps = 1e-3;
    Problem sstokes(sw, sa, sz, sb);
    sstokes = Integral(Jacobian(sw), Jacobian(sz)) - Integral(sa, Div(sz)) +
      Integral(Div(sw), sb) + Integral(eps * sa, sb) - Integral(f, sb) +
      DirichletBC(sw, VectorFunction{0, 0}).on(1) + DirichletBC(sa, g).on(5) +
      PeriodicBC(sa, pairs);
    Problem pstokes(pw, pa, pz, pb);
    pstokes = Integral(Jacobian(pw), Jacobian(pz)) - Integral(pa, Div(pz)) +
      Integral(Div(pw), pb) + Integral(eps * pa, pb) - Integral(f, pb) +
      DirichletBC(pw, VectorFunction{0, 0}).on(1) + DirichletBC(pa, g).on(5) +
      PeriodicBC(pa, pairs);
    expectSameSystem(sstokes, pstokes);

    Problem scoupled(su, sa, sv, sb);
    scoupled = Integral(Grad(su), Grad(sv)) + Integral(sa, sv) +
      Integral(Grad(sa), Grad(sb)) - Integral(f, sv) - Integral(f, sb) +
      DirichletBC(su, Real(0.5) * sa, RealFunction(2)).on(1) +
      DirichletBC(sa, RealFunction(1)).on(2);
    Problem pcoupled(pu, pa, pv, pb);
    pcoupled = Integral(Grad(pu), Grad(pv)) + Integral(pa, pv) +
      Integral(Grad(pa), Grad(pb)) - Integral(f, pv) - Integral(f, pb) +
      DirichletBC(pu, Real(0.5) * pa, RealFunction(2)).on(1) +
      DirichletBC(pa, RealFunction(1)).on(2);
    expectSameSystem(scoupled, pcoupled);

    pcoupled.assemble();
    expectSameSystem(scoupled, pcoupled);
  }

  /// @brief Device PETSc forms and problems with preassembled forms match the
  /// sparse device assembly.
  /// @param source Mesh tagged by withTaggedFace().
  void checkPreassembled(const HostMesh& source)
  {
    ASSERT_EQ(source.getSpaceDimension(), 2);
    const DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);
    RealFunction f(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 1 + p.x() * p.y();
    });
    RealFunction g(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 2 - p.x(); });

    BilinearForm sm(su, sv);
    sm = Integral(su, sv) + BoundaryIntegral(su, sv).over(2);
    sm.assemble();
    BilinearForm pm(pu, pv);
    pm = Integral(pu, pv) + BoundaryIntegral(pu, pv).over(2);
    pm.assemble();
    EXPECT_TRUE(hasType(pm.getOperator(), MATSEQAIJKOKKOS));
    expectSameOperator(sm.getOperator(), toSparse(pm.getOperator()));
    LinearForm sl(sv);
    sl = Integral(f, sv);
    sl.assemble();
    LinearForm pl(pv);
    pl = Integral(f, pv);
    pl.assemble();
    EXPECT_TRUE(hasType(pl.getVector(), VECSEQKOKKOS));
    expectSameVector(sl.getVector(), toVector(pl.getVector()));

    const auto pairs = makePairs(fes.getSize(), fixedDOFs(su, 5));
    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv) +
      BoundaryIntegral(su, sv).over(2) - Integral(f, sv) + DirichletBC(su, g).on(5) +
      PeriodicBC(su, pairs);
    using Body = ProblemBody<::Mat, ::Vec, PetscScalar>;
    const BilinearFormBase<::Mat>& mass = pm;
    const LinearFormBase<::Vec>& load = pl;
    Problem pp(pu, pv);
    pp = Body(Integral(Grad(pu), Grad(pv))) + mass - load + DirichletBC(pu, g).on(5) +
      PeriodicBC(pu, pairs);
    expectSameSystem(sp, pp);
    pp.assemble();
    expectSameSystem(sp, pp);

    P1 vh(mesh, 2);
    TrialFunction sw(vh);
    TestFunction sz(vh);
    PETSc::Variational::TrialFunction pw(vh);
    PETSc::Variational::TestFunction pz(vh);
    const Real eps = 1e-3;
    BilinearForm pe(pu, pv);
    pe = Integral(eps * pu, pv);
    pe.assemble();
    const BilinearFormBase<::Mat>& penalty = pe;
    Problem sstokes(sw, su, sz, sv);
    sstokes = Integral(Jacobian(sw), Jacobian(sz)) - Integral(su, Div(sz)) +
      Integral(Div(sw), sv) + Integral(eps * su, sv) - Integral(f, sv) +
      DirichletBC(sw, VectorFunction{0, 0}).on(1) + DirichletBC(su, g).on(5);
    Problem pstokes(pw, pu, pz, pv);
    pstokes = Body(Integral(Jacobian(pw), Jacobian(pz))) - Integral(pu, Div(pz)) +
      Integral(Div(pw), pv) + penalty - load +
      DirichletBC(pw, VectorFunction{0, 0}).on(1) + DirichletBC(pu, g).on(5);
    expectSameSystem(sstokes, pstokes);
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

    // The right-hand side, eliminated values included, is summed in a fixed
    // order.
    for (int k = 0; k < 4; ++k)
    {
      pb.assemble();
      const auto again = toVector(axb.getVector());
      ASSERT_EQ(again.size(), b.size());
      for (Eigen::Index i = 0; i < b.size(); ++i)
        EXPECT_EQ(again(i), b(i)) << "entry " << i;
    }

    ASSERT_EQ(VecZeroEntries(axb.getVector()), PETSC_SUCCESS);
    pb.assemble(AssemblyTarget::RHS);
    expectSameVector(b, toVector(axb.getVector()));
    ASSERT_EQ(MatZeroEntries(axb.getOperator()), PETSC_SUCCESS);
    pb.assemble(AssemblyTarget::LHS);
    expectSameOperator(A, toSparse(axb.getOperator()));
  }

  /**
   * @brief Integrators restricted to attributes follow the attributes of
   * the cells between assemblies, while the pattern is kept when they do
   * not move.
   */
  void checkMovedAttributes(HostMesh source)
  {
    const size_t d = source.getDimension();
    const Index cells = source.getCellCount();
    for (Index i = 0; i < cells; ++i)
      source.setAttribute({d, i}, Geometry::Attribute(1 + i % 2));
    DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);
    RealFunction f(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() + p.y(); });
    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv).over(1) +
      BoundaryIntegral(su, sv).over(2) - Integral(f, sv) +
      DirichletBC(su, RealFunction(1)).on(1);
    Problem pp(pu, pv);
    pp = Integral(Grad(pu), Grad(pv)) + Integral(pu, pv).over(1) +
      BoundaryIntegral(pu, pv).over(2) - Integral(f, pv) +
      DirichletBC(pu, RealFunction(1)).on(1);
    expectSameSystem(sp, pp);
    auto& axb = pp.getLinearSystem();
    PetscObjectState state = 0;
    ASSERT_EQ(MatGetNonzeroState(axb.getOperator(), &state), PETSC_SUCCESS);
    expectSameSystem(sp, pp);
    PetscObjectState again = 0;
    ASSERT_EQ(MatGetNonzeroState(axb.getOperator(), &again), PETSC_SUCCESS);
    EXPECT_EQ(state, again);

    // As many cells keep the attribute 1, so that on a mesh of one cell
    // type the restricted entries keep their count but not their
    // coordinates.
    for (Index i = 0; i < cells; ++i)
      mesh.setAttribute({d, i}, Geometry::Attribute(2 - i % 2));
    expectSameSystem(sp, pp);
    expectSameSystem(sp, pp);

    // Back to the first attributes, then to fewer cells with the attribute 1.
    for (Index i = 0; i < cells; ++i)
      mesh.setAttribute({d, i}, Geometry::Attribute(1 + i % 2));
    expectSameSystem(sp, pp);
    for (Index i = 0; i < cells; ++i)
      mesh.setAttribute({d, i}, Geometry::Attribute(i % 4 == 0 ? 1 : 2));
    expectSameSystem(sp, pp);
    expectSameSystem(sp, pp);
  }

  /**
   * @brief The eliminated entries follow the fixed degrees of freedom when
   * they move between assemblies while the pattern is kept.
   */
  void checkMovedConstraints(const HostMesh& source)
  {
    const size_t d = source.getDimension();
    DeviceMesh mesh(TestContext{}, source);
    P1 fes(mesh);
    TrialFunction su(fes);
    TestFunction sv(fes);
    PETSc::Variational::TrialFunction pu(fes);
    PETSc::Variational::TestFunction pv(fes);
    RealFunction f(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return p.x() + p.y(); });
    RealFunction g(
      KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real { return 2 - p.x(); });
    Problem sp(su, sv);
    sp = Integral(Grad(su), Grad(sv)) + Integral(su, sv) - Integral(f, sv) +
      DirichletBC(su, g).on(1);
    Problem pp(pu, pv);
    pp = Integral(Grad(pu), Grad(pv)) + Integral(pu, pv) - Integral(f, pv) +
      DirichletBC(pu, g).on(1);
    expectSameSystem(sp, pp);
    expectSameSystem(sp, pp);
    auto& axb = pp.getLinearSystem();
    PetscObjectState state = 0;
    ASSERT_EQ(MatGetNonzeroState(axb.getOperator(), &state), PETSC_SUCCESS);

    // Other faces carry the condition, then the first ones again.
    const Index faces = mesh.getFaceCount();
    for (Index i = 0; i < faces; ++i)
      mesh.setAttribute({d - 1, i}, Geometry::Attribute(1 + (i + 1) % 3));
    expectSameSystem(sp, pp);
    expectSameSystem(sp, pp);
    for (Index i = 0; i < faces; ++i)
      mesh.setAttribute({d - 1, i}, Geometry::Attribute(1 + i % 3));
    expectSameSystem(sp, pp);
    PetscObjectState again = 0;
    ASSERT_EQ(MatGetNonzeroState(axb.getOperator(), &again), PETSC_SUCCESS);
    EXPECT_EQ(state, again);
  }

  /**
   * @brief Unconstrained right-hand sides, summed from the entries of the
   * linear integrators in place, match the sparse device assembly when
   * assembled again and when the linear integrators change.
   */
  void checkUnconstrained(const HostMesh& source)
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
    RealFunction g(KOKKOS_LAMBDA(const DeviceGeometry::PointView& p)->Real {
      return 2 * p.x() - p.y() + 1;
    });
    const Real eps = 1e-3;

    Problem sstokes(su, sp, sv, sq);
    Problem pstokes(pu, pp, pv, pq);
    const auto stokes = [&] {
      sstokes = Integral(Jacobian(su), Jacobian(sv)) - Integral(sp, Div(sv)) +
        Integral(Div(su), sq) + Integral(eps * sp, sq) - Integral(f, sv) -
        Integral(s, sq);
      pstokes = Integral(Jacobian(pu), Jacobian(pv)) - Integral(pp, Div(pv)) +
        Integral(Div(pu), pq) + Integral(eps * pp, pq) - Integral(f, pv) -
        Integral(s, pq);
    };
    stokes();
    expectSameSystem(sstokes, pstokes);
    const auto b = toVector(pstokes.getLinearSystem().getVector());
    for (int k = 0; k < 3; ++k)
    {
      pstokes.assemble();
      const auto again = toVector(pstokes.getLinearSystem().getVector());
      ASSERT_EQ(again.size(), b.size());
      for (Eigen::Index i = 0; i < b.size(); ++i)
        EXPECT_EQ(again(i), b(i)) << "entry " << i;
    }

    // Constrained, then unconstrained again with the same integrators.
    sstokes = Integral(Jacobian(su), Jacobian(sv)) - Integral(sp, Div(sv)) +
      Integral(Div(su), sq) + Integral(eps * sp, sq) - Integral(f, sv) - Integral(s, sq) +
      DirichletBC(su, VectorFunction{0, 0}).on(1);
    pstokes = Integral(Jacobian(pu), Jacobian(pv)) - Integral(pp, Div(pv)) +
      Integral(Div(pu), pq) + Integral(eps * pp, pq) - Integral(f, pv) - Integral(s, pq) +
      DirichletBC(pu, VectorFunction{0, 0}).on(1);
    expectSameSystem(sstokes, pstokes);
    stokes();
    expectSameSystem(sstokes, pstokes);

    // Other linear integrators, one of them restricted to attributes.
    sstokes = Integral(Jacobian(su), Jacobian(sv)) - Integral(sp, Div(sv)) +
      Integral(Div(su), sq) + Integral(eps * sp, sq) - Integral(s, sq) - Integral(f, sv) +
      BoundaryIntegral(g, sq).over(2);
    pstokes = Integral(Jacobian(pu), Jacobian(pv)) - Integral(pp, Div(pv)) +
      Integral(Div(pu), pq) + Integral(eps * pp, pq) - Integral(s, pq) - Integral(f, pv) +
      BoundaryIntegral(g, pq).over(2);
    expectSameSystem(sstokes, pstokes);
    expectSameSystem(sstokes, pstokes);
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

    using HostSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    Problem<HostSystem, decltype(su), decltype(sv)> sp(su, sv);
    sp =
      Integral(Grad(su), Grad(sv)) - Integral(f, sv) + DirichletBC(su, RealFunction(0));
    Solver::SparseLU(sp).solve();
    const Math::Vector<Real> expected = su.getSolution().getData();

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

    // The PETSc solution enters a problem as a coefficient, through its
    // values and its gradient, as the same coefficients in an Eigen vector.
    GridFunction gs(fes);
    gs.getData() = solution;
    const auto& gp = pu.getSolution();
    Problem sq(su, sv);
    sq = Integral(Grad(su), Grad(sv)) + Integral(gs * su, sv) - Integral(gs, sv) -
      Integral(Dot(Grad(gs), Grad(sv)));
    Problem pq(pu, pv);
    pq = Integral(Grad(pu), Grad(pv)) + Integral(gp * pu, pv) - Integral(gp, pv) -
      Integral(Dot(Grad(gp), Grad(pv)));
    expectSameSystem(sq, pq);

    // Coefficients written again through PETSc are read again.
    ASSERT_EQ(VecScale(gp.getData(), 2), PETSC_SUCCESS);
    ASSERT_EQ(VecShift(gp.getData(), 1), PETSC_SUCCESS);
    gs.getData() = (2 * solution).array() + 1;
    expectSameSystemAssembledFirst(sq, pq);
    PetscScalar* values = nullptr;
    ASSERT_EQ(VecGetArray(gp.getData(), &values), PETSC_SUCCESS);
    values[0] = -3;
    ASSERT_EQ(VecRestoreArray(gp.getData(), &values), PETSC_SUCCESS);
    gs.getData()(0) = -3;
    expectSameSystemAssembledFirst(sq, pq);
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

  TEST(Kokkos_PETSc_Problem, PotentialOnMixedMesh)
  {
    checkPotential(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, ConstraintsOnTriangles)
  {
    checkConstraints(withTaggedFace(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, PreassembledOnTriangles)
  {
    checkPreassembled(withTaggedFace(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, Reassembly)
  {
    checkReassembly(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, MovedAttributes)
  {
    checkMovedAttributes(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, MovedAttributesOnMixedMesh)
  {
    checkMovedAttributes(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, MovedConstraints)
  {
    checkMovedConstraints(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, MovedConstraintsOnMixedMesh)
  {
    checkMovedConstraints(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, UnconstrainedOnTriangles)
  {
    checkUnconstrained(withFaces(makeTriangleMesh()));
  }

  TEST(Kokkos_PETSc_Problem, UnconstrainedOnMixedMesh)
  {
    checkUnconstrained(withFaces(makeMixedMesh()));
  }

  TEST(Kokkos_PETSc_Problem, SolveOnTriangles)
  {
    checkSolve(withFaces(makeTriangleMesh()));
  }

  /// @brief COO duplicates follow integrator rank for every row bucket size and permutation width.
  template <class Position>
  void checkCOODuplicateOrder()
  {
    using Assembler = Assembly::Internal::COOAssembly<TestBackend>;
    using Segment = Assembly::Internal::COOSegment;
    using Memory = typename Assembler::MemorySpace;
    using Indices = ::Kokkos::View<Index*, Memory>;
    using Order = ::Kokkos::View<Position*, Memory>;
    const TestExecutionSpace execution;
    constexpr size_t count = 768;
    constexpr size_t columns = 11;
    constexpr size_t smallBucket = 16;
    constexpr Index coordinatePeriod = 23 * columns;
    std::vector<Segment> segments(2);
    segments[0].spacing = 4;
    segments[1].spacing = 2;
    const std::vector<size_t> offsets{0, count / 3, 2 * count / 3};
    const auto rows =
      Kokkos::Internal::upload<Index>(execution, count, "test rows", [](Index* out) {
        for (Index i = 0; i < count; ++i)
          out[i] = (17 * i) % 23;
      });
    const auto cols =
      Kokkos::Internal::upload<Index>(execution, count, "test columns", [](Index* out) {
        for (Index i = 0; i < count; ++i)
          out[i] = (7 * i) % columns;
      });
    const auto hostRows =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, rows);
    const auto hostCols =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cols);
    std::vector<Index> expected;
    for (size_t s = 0; s <= segments.size(); ++s)
    {
      const size_t begin = offsets[s];
      const size_t size = (s < segments.size() ? offsets[s + 1] : count) - begin;
      const size_t spacing = s < segments.size() ? segments[s].spacing : 1;
      for (size_t r = 0; r < size; ++r)
      {
        expected.push_back(
          begin + Kokkos::Variational::getEntryPosition(r, size, spacing));
      }
    }
    std::stable_sort(expected.begin(), expected.end(), [&](Index lhs, Index rhs) {
      return std::pair{hostRows(lhs), hostCols(lhs)} <
        std::pair{hostRows(rhs), hostCols(rhs)};
    });
    for (const size_t bucket : {count, smallBucket})
    {
      SCOPED_TRACE(bucket);
      const Order order("test COO order", count);
      Assembler::sortPattern(
        execution, columns, segments, offsets, rows, cols, order, bucket);
      const auto hostOrder =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, order);
      for (size_t i = 0; i < count; ++i)
        ASSERT_EQ(hostOrder(i), expected[i]) << "entry " << i;
      const Indices starts = Assembler::distinctStarts(execution, rows, cols, order);
      const auto values = Kokkos::Internal::upload<PetscScalar>(
        execution, count, "test values", [](PetscScalar* out) {
          // Equal coordinates repeat every 23 * 11 entries. Their values
          // expose a change of accumulation order through cancellation.
          for (Index i = 0; i < count; ++i)
          {
            out[i] = i / coordinatePeriod == 0 ? 1e16
              : i / coordinatePeriod == 1      ? 1
                                               : -1e16;
          }
        });
      const ::Kokkos::View<PetscScalar*, Memory> merged(
        "test merged", starts.extent(0) - 1);
      Assembler::merge(execution, order, starts, values, merged);
      const auto hostValues =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, values);
      const auto hostMerged =
        ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, merged);
      size_t group = 0;
      for (size_t begin = 0; begin < count; ++group)
      {
        size_t end = begin + 1;
        const Index first = expected[begin];
        while (end < count && hostRows(expected[end]) == hostRows(first) &&
          hostCols(expected[end]) == hostCols(first))
          ++end;
        PetscScalar sum = 0;
        for (size_t i = begin; i < end; ++i)
          sum += hostValues(expected[i]);
        EXPECT_EQ(hostMerged(group), sum);
        begin = end;
      }
    }
  }

  TEST(Kokkos_PETSc_Problem, COODuplicateOrder)
  {
    checkCOODuplicateOrder<std::uint32_t>();
    checkCOODuplicateOrder<Index>();
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
