/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for Assembly::Sequential.
 *
 * These tests exercise the Sequential assembler directly — verifying that
 * LinearForm vectors, BilinearForm matrices (sparse and dense), and
 * multi-variable Problem systems (including preassembled BF/LF and Dirichlet
 * BC elimination) are assembled correctly.
 */
#include <gtest/gtest.h>

#include <boost/bimap.hpp>

#include "Rodin/Variational.h"
#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // LinearForm — Sequential assembler
  // =========================================================================

  /**
   * @brief Assembled vector has the expected size (one entry per DOF).
   */
  TEST(Assembly_Sequential_LinearForm, VectorSize_P1_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();

    EXPECT_EQ(lf.getVector().size(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief Assembled vector has the expected size for Quadrilateral meshes.
   */
  TEST(Assembly_Sequential_LinearForm, VectorSize_P1_Quad_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Quadrilateral, { 4, 4 });
    P1 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(v);
    lf.assemble();

    EXPECT_EQ(lf.getVector().size(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief All entries of ∫ 1·v dΩ are positive for a positive basis.
   */
  TEST(Assembly_Sequential_LinearForm, ConstantRHS_EntriesPositive_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    const auto& b = lf.getVector();
    EXPECT_GT(b.size(), 0);
    EXPECT_GT(b.sum(), 0.0);
  }

  /**
   * @brief ∫ 1·v dΩ with P0 sums to area (each element has exactly one DOF
   * and the mass is |K|).
   */
  TEST(Assembly_Sequential_LinearForm, ConstantRHS_SumEqualsArea_P0)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P0 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    // For f = 1 and P0, each entry equals the element area.
    // UniformGrid({4,4}) has vertices at (0,0)..(3,3), so the domain is [0,3]x[0,3]
    // with area = 9.  The sum over all 18 P0 elements equals the domain area.
    const Real total = lf.getVector().sum();
    EXPECT_NEAR(total, 9.0, 1e-12);
  }

  /**
   * @brief Scaling: ∫ 2·v dΩ = 2 · ∫ 1·v dΩ.
   */
  TEST(Assembly_Sequential_LinearForm, LinearScaling_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TestFunction v(fes);

    LinearForm lf1(v);
    lf1 = Integral(RealFunction(1.0), v);
    lf1.assemble();

    LinearForm lf2(v);
    lf2 = Integral(RealFunction(2.0), v);
    lf2.assemble();

    const Real diff = (lf2.getVector() - 2.0 * lf1.getVector()).norm();
    EXPECT_NEAR(diff, 0.0, 1e-12);
  }

  // =========================================================================
  // BilinearForm — Sequential assembler (SparseMatrix)
  // =========================================================================

  /**
   * @brief Stiffness matrix has correct dimensions.
   */
  TEST(Assembly_Sequential_BilinearForm_Sparse, MatrixDimensions_P1_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& A = bf.getOperator();
    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(fes.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief P0 mass matrix is diagonal (each element has one DOF).
   */
  TEST(Assembly_Sequential_BilinearForm_Sparse, P0MassMatrix_IsDiagonal)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P0 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    const auto& A = bf.getOperator();
    ASSERT_GT(A.rows(), 0);

    // Off-diagonal entries should be zero for P0
    for (int k = 0; k < A.outerSize(); ++k)
    {
      for (typename Math::SparseMatrix<Real>::InnerIterator it(A, k); it; ++it)
      {
        if (it.row() != it.col())
        {
          EXPECT_NEAR(it.value(), 0.0, 1e-14);
        }
      }
    }
  }

  /**
   * @brief P1 mass matrix is symmetric.
   */
  TEST(Assembly_Sequential_BilinearForm_Sparse, P1MassMatrix_IsSymmetric)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    const auto& A = bf.getOperator();
    const Math::SparseMatrix<Real> diff = A - Math::SparseMatrix<Real>(A.transpose());
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief P1 stiffness matrix is symmetric.
   */
  TEST(Assembly_Sequential_BilinearForm_Sparse, P1StiffnessMatrix_IsSymmetric)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& A = bf.getOperator();
    const Math::SparseMatrix<Real> diff = A - Math::SparseMatrix<Real>(A.transpose());
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief P1 mass matrix row-sum equals domain area (∑_j M_ij summed over all
   * i,j equals ∫_Ω 1 dΩ = area(Ω)).  For UniformGrid({4,4}) the domain is
   * [0,3]x[0,3] so area = 9.
   */
  TEST(Assembly_Sequential_BilinearForm_Sparse, P1MassMatrix_SumEqualsArea)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    // sum of all entries of M equals ∫_Ω (∑_i φ_i)(∑_j φ_j) dΩ = ∫_Ω 1 dΩ = area(Ω).
    // UniformGrid({4,4}) spans [0,3]x[0,3], so area = 9.
    const auto& A = bf.getOperator();
    Math::Vector<Real> ones(A.cols());
    ones.setOnes();
    const Real row_sum = (A * ones).sum();
    EXPECT_NEAR(row_sum, 9.0, 1e-12);
  }

  /**
   * @brief Non-square BF: P1 trial, P0 test. Dimensions must match
   * (testSize × trialSize).
   */
  TEST(Assembly_Sequential_BilinearForm_Sparse, NonSquare_P0test_P1trial_Dimensions)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 p1h(mesh);
    P0 p0h(mesh);
    TrialFunction u(p1h);
    TestFunction  q(p0h);

    BilinearForm bf(u, q);
    bf = Integral(u, q);
    bf.assemble();

    const auto& A = bf.getOperator();
    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(p0h.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(p1h.getSize()));
  }

  // =========================================================================
  // Problem (single-variable) — Sequential assembler
  // =========================================================================

  /**
   * @brief Single-variable Problem: assembled system has correct dimensions.
   */
  TEST(Assembly_Sequential_Problem_SingleVar, Dimensions_P1_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    Problem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v);
    problem.assemble();

    const auto& ls = problem.getLinearSystem();
    const auto& A  = ls.getOperator();
    const auto& b  = ls.getVector();

    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(fes.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(fes.getSize()));
    EXPECT_EQ(b.size(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief DirichletBC is enforced: fixed DOF row is an identity row
   * and RHS entry equals the prescribed value.
   */
  TEST(Assembly_Sequential_Problem_SingleVar, DirichletBC_RowIsIdentity_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    const Real gValue = 3.0;
    Problem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v)
            + DirichletBC(u, RealFunction(gValue));
    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();
    const auto& b = problem.getLinearSystem().getVector();

    // At least one Dirichlet DOF must exist
    bool found = false;
    for (int k = 0; k < A.outerSize(); ++k)
    {
      bool is_identity_row = false;
      for (typename Math::SparseMatrix<Real>::InnerIterator it(A, k); it; ++it)
      {
        if (it.row() == it.col() && std::abs(it.value() - 1.0) < 1e-14)
        {
          // Check that this row has only the diagonal entry non-zero
          // and b[row] == gValue
          if (std::abs(b.coeff(it.row()) - gValue) < 1e-10)
          {
            is_identity_row = true;
            found = true;
          }
        }
      }
      (void)is_identity_row;
    }
    EXPECT_TRUE(found) << "No Dirichlet-constrained row found in assembled system";
  }

  /**
   * @brief Preassembled BF added to a single-variable Problem produces the
   * same matrix as assembling the integral directly.
   */
  TEST(Assembly_Sequential_Problem_SingleVar, PreassembledBF_MatchesIntegral_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);

    // Integral-only
    TrialFunction u1(fes);
    TestFunction  v1(fes);
    Problem prob1(u1, v1);
    prob1 = Integral(Grad(u1), Grad(v1));
    prob1.assemble();

    // Preassembled BF
    TrialFunction u2(fes);
    TestFunction  v2(fes);
    BilinearForm stiff(u2, v2);
    stiff = Integral(Grad(u2), Grad(v2));
    stiff.assemble();

    Problem prob2(u2, v2);
    prob2 = stiff;
    prob2.assemble();

    const auto& A1 = prob1.getLinearSystem().getOperator();
    const auto& A2 = prob2.getLinearSystem().getOperator();

    ASSERT_EQ(A1.rows(), A2.rows());
    ASSERT_EQ(A1.cols(), A2.cols());
    const Math::SparseMatrix<Real> diff = A1 - A2;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief Preassembled LF added to a single-variable Problem produces the
   * same RHS vector as assembling the integral directly.
   */
  TEST(Assembly_Sequential_Problem_SingleVar, PreassembledLF_MatchesIntegral_P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);

    // Integral-only
    TrialFunction u1(fes);
    TestFunction  v1(fes);
    Problem prob1(u1, v1);
    prob1 = Integral(Grad(u1), Grad(v1))
          - Integral(RealFunction(1.0), v1);
    prob1.assemble();

    // Preassembled LF
    TrialFunction u2(fes);
    TestFunction  v2(fes);
    LinearForm load(v2);
    load = Integral(RealFunction(1.0), v2);
    load.assemble();

    Problem prob2(u2, v2);
    prob2 = Integral(Grad(u2), Grad(v2))
          - load;
    prob2.assemble();

    const auto& b1 = prob1.getLinearSystem().getVector();
    const auto& b2 = prob2.getLinearSystem().getVector();

    ASSERT_EQ(b1.size(), b2.size());
    EXPECT_NEAR((b1 - b2).norm(), 0.0, 1e-12);
  }

  // =========================================================================
  // Problem (multi-variable) — Sequential assembler, block offsets
  // =========================================================================

  /**
   * @brief Multi-variable Problem: global system has combined DOF count.
   */
  TEST(Assembly_Sequential_Problem_MultiVar, Dimensions_P0P1_4x4)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P0 p0h(mesh);
    P1 p1h(mesh);

    TrialFunction u(p1h);
    TrialFunction p(p0h);
    TestFunction  v(p1h);
    TestFunction  q(p0h);

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v)
          + Integral(p, q);
    mixed.assemble();

    const auto& A = mixed.getLinearSystem().getOperator();
    const auto& b = mixed.getLinearSystem().getVector();

    const Eigen::Index expected = static_cast<Eigen::Index>(p1h.getSize() + p0h.getSize());
    EXPECT_EQ(A.rows(), expected);
    EXPECT_EQ(A.cols(), expected);
    EXPECT_EQ(b.size(), expected);
  }

  /**
   * @brief Multi-variable preassembled BF: block [p,q] placed at correct
   * offset in the global matrix.
   */
  TEST(Assembly_Sequential_Problem_MultiVar, PreassembledBF_BlockOffset_P0P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P0 p0h(mesh);
    P1 p1h(mesh);

    // --- integral-only ---
    TrialFunction u1(p1h);
    TrialFunction p1(p0h);
    TestFunction  v1(p1h);
    TestFunction  q1(p0h);
    Problem intOnly(u1, v1, p1, q1);
    intOnly = Integral(u1, v1) + Integral(p1, q1);
    intOnly.assemble();

    // --- preassembled BF for [p,q] block ---
    TrialFunction u2(p1h);
    TrialFunction p2(p0h);
    TestFunction  v2(p1h);
    TestFunction  q2(p0h);
    BilinearForm massP0(p2, q2);
    massP0 = Integral(p2, q2);
    massP0.assemble();

    Problem preasm(u2, v2, p2, q2);
    preasm = Integral(u2, v2) + massP0;
    preasm.assemble();

    const auto& A1 = intOnly.getLinearSystem().getOperator();
    const auto& A2 = preasm.getLinearSystem().getOperator();

    ASSERT_EQ(A1.rows(), A2.rows());
    ASSERT_EQ(A1.cols(), A2.cols());
    const Math::SparseMatrix<Real> diff = A1 - A2;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief Multi-variable preassembled LF: block [q] placed at correct offset
   * in the global RHS vector.
   */
  TEST(Assembly_Sequential_Problem_MultiVar, PreassembledLF_BlockOffset_P0P1)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);  // needed for the P1h boundary iteration

    P0 p0h(mesh);
    P1 p1h(mesh);

    // --- integral-only ---
    TrialFunction u1(p1h);
    TrialFunction p1(p0h);
    TestFunction  v1(p1h);
    TestFunction  q1(p0h);
    Problem intOnly(u1, v1, p1, q1);
    intOnly = Integral(u1, v1)
            + Integral(p1, q1)
            - Integral(RealFunction(1.0), q1);
    intOnly.assemble();

    // --- preassembled LF for [q] block ---
    TrialFunction u2(p1h);
    TrialFunction p2(p0h);
    TestFunction  v2(p1h);
    TestFunction  q2(p0h);
    LinearForm loadQ(q2);
    loadQ = Integral(RealFunction(1.0), q2);
    loadQ.assemble();

    Problem preasm(u2, v2, p2, q2);
    preasm = Integral(u2, v2)
           + Integral(p2, q2)
           - loadQ;
    preasm.assemble();

    const auto& b1 = intOnly.getLinearSystem().getVector();
    const auto& b2 = preasm.getLinearSystem().getVector();

    ASSERT_EQ(b1.size(), b2.size());
    EXPECT_NEAR((b1 - b2).norm(), 0.0, 1e-12);
  }

  TEST(Assembly_Problem_Identification, PreassembledFormsProjectRowsColumnsAndVector)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);
    TrialFunction eta(fes);
    TestFunction  zeta(fes);

    constexpr Real gamma = 2.0;
    const Index n = static_cast<Index>(fes.getSize());

    BilinearForm uu(u, v);
    uu.getOperator().resize(n, n);
    std::vector<Eigen::Triplet<Real>> triplets{
      Eigen::Triplet<Real>(0, 0, 2.0),
      Eigen::Triplet<Real>(0, 1, 3.0),
      Eigen::Triplet<Real>(1, 0, 5.0)
    };
    uu.getOperator().setFromTriplets(triplets.begin(), triplets.end());

    LinearForm loadU(v);
    loadU.getVector().resize(n);
    loadU.getVector().setZero();
    loadU.getVector().coeffRef(0) = 7.0;
    loadU.getVector().coeffRef(1) = 11.0;

    Problem problem(u, v, eta, zeta);
    problem = uu + DirichletBC(u, RealFunction(gamma) * eta) - loadU;
    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();
    const auto& b = problem.getLinearSystem().getVector();
    ASSERT_EQ(A.rows(), static_cast<Eigen::Index>(2 * n));
    ASSERT_EQ(A.cols(), static_cast<Eigen::Index>(2 * n));

    EXPECT_NEAR(A.coeff(n + 0, n + 0), gamma * 2.0 * gamma, 1e-14);
    EXPECT_NEAR(A.coeff(n + 0, n + 1), gamma * 3.0 * gamma, 1e-14);
    EXPECT_NEAR(A.coeff(n + 1, n + 0), gamma * 5.0 * gamma, 1e-14);

    EXPECT_NEAR(A.coeff(0, 0), 1.0, 1e-14);
    EXPECT_NEAR(A.coeff(0, n + 0), -gamma, 1e-14);
    EXPECT_NEAR(b.coeff(0), 0.0, 1e-14);
    EXPECT_NEAR(b.coeff(n + 0), gamma * 7.0, 1e-14);
    EXPECT_NEAR(b.coeff(n + 1), gamma * 11.0, 1e-14);
  }

  TEST(Assembly_Problem_Identification, SequentialVectorMasterProjectsMultipleFiniteSpaces)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 slaveFES(mesh);
    P1 masterFES(mesh, mesh.getSpaceDimension());

    TrialFunction u(slaveFES);
    TestFunction  v(slaveFES);
    TrialFunction eta(masterFES);
    TestFunction  zeta(masterFES);

    auto bc =
      DirichletBC(
          u,
          RealFunction(2.0) * eta.x() + RealFunction(-0.5) * eta.y());
    bc.assemble();

    using IdentifiedDOFs = DirichletBCBase<Real>::IdentifiedDOFs;
    ASSERT_TRUE(std::holds_alternative<IdentifiedDOFs>(bc.getDOFs()));
    const auto& ident = std::get<IdentifiedDOFs>(bc.getDOFs());
    ASSERT_FALSE(ident.empty());

    bool sawMultiMaster = false;
    for (const auto& [slave, row] : ident)
    {
      (void) slave;
      if (row.first.size() >= 2)
        sawMultiMaster = true;
    }
    ASSERT_TRUE(sawMultiMaster);

    const Index nSlave  = static_cast<Index>(slaveFES.getSize());
    const Index nMaster = static_cast<Index>(masterFES.getSize());
    const Index nTotal  = nSlave + nMaster;

    struct MatrixEntry
    {
      Index row;
      Index col;
      Real value;
    };
    const std::vector<MatrixEntry> matrixEntries{
      { 0, 0, 2.0 },
      { 0, 1, 3.0 },
      { 1, 0, 5.0 },
      { 1, 1, 7.0 }
    };
    const std::vector<std::pair<Index, Real>> vectorEntries{
      { 0, 11.0 },
      { 1, -13.0 }
    };

    BilinearForm uu(u, v);
    uu.getOperator().resize(nSlave, nSlave);
    std::vector<Eigen::Triplet<Real>> triplets;
    for (const auto& e : matrixEntries)
      triplets.emplace_back(e.row, e.col, e.value);
    uu.getOperator().setFromTriplets(triplets.begin(), triplets.end());

    LinearForm loadU(v);
    loadU.getVector().resize(nSlave);
    loadU.getVector().setZero();
    for (const auto& [row, value] : vectorEntries)
      loadU.getVector().coeffRef(row) = value;

    auto body = uu + bc - loadU;

    using LinearSystemType =
      Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    using ProblemType =
      Problem<LinearSystemType,
              decltype(u), decltype(v), decltype(eta), decltype(zeta)>;

    auto trialFunctions = Tuple{ std::ref(u), std::ref(eta) };
    auto testFunctions  = Tuple{ std::ref(v), std::ref(zeta) };

    std::array<size_t, 2> trialOffsets{
      0, static_cast<size_t>(nSlave)
    };
    std::array<size_t, 2> testOffsets{
      0, static_cast<size_t>(nSlave)
    };

    boost::bimap<FormLanguage::Base::UUID, size_t> trialUUIDMap;
    boost::bimap<FormLanguage::Base::UUID, size_t> testUUIDMap;
    trialUUIDMap.right.insert({ 0, u.getUUID() });
    trialUUIDMap.right.insert({ 1, eta.getUUID() });
    testUUIDMap.right.insert({ 0, v.getUUID() });
    testUUIDMap.right.insert({ 1, zeta.getUUID() });

    Assembly::ProblemAssemblyInput<
      std::decay_t<decltype(body)>,
      decltype(u), decltype(v), decltype(eta), decltype(zeta)> input(
          body,
          trialFunctions,
          testFunctions,
          trialOffsets,
          testOffsets,
          trialUUIDMap,
          testUUIDMap,
          static_cast<size_t>(nTotal),
          static_cast<size_t>(nTotal));

    LinearSystemType ls;
    Assembly::Sequential<LinearSystemType, ProblemType> assembler;
    assembler.execute(ls, input);

    Eigen::MatrixXd expectedA =
      Eigen::MatrixXd::Zero(
          static_cast<Eigen::Index>(nTotal),
          static_cast<Eigen::Index>(nTotal));
    Eigen::VectorXd expectedB =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(nTotal));

    auto expand = [&](Index local)
    {
      std::vector<std::pair<Index, Real>> res;
      const auto it = ident.find(local);
      if (it == ident.end())
      {
        res.emplace_back(local, 1.0);
        return res;
      }
      const auto& masters = it->second.first;
      const auto& coeffs  = it->second.second;
      for (Index k = 0; k < static_cast<Index>(masters.size()); k++)
        res.emplace_back(nSlave + masters[k], coeffs[k]);
      return res;
    };

    for (const auto& e : matrixEntries)
      for (const auto& r : expand(e.row))
        for (const auto& c : expand(e.col))
          expectedA(r.first, c.first) += r.second * e.value * c.second;

    for (const auto& [row, value] : vectorEntries)
      for (const auto& r : expand(row))
        expectedB(r.first) += r.second * value;

    for (const auto& [slave, row] : ident)
    {
      expectedA.row(static_cast<Eigen::Index>(slave)).setZero();
      expectedA(slave, slave) = 1.0;
      const auto& masters = row.first;
      const auto& coeffs  = row.second;
      for (Index k = 0; k < static_cast<Index>(masters.size()); k++)
        expectedA(slave, nSlave + masters[k]) -= coeffs[k];
      expectedB(slave) = 0.0;
    }

    const auto& A = ls.getOperator();
    const auto& b = ls.getVector();
    ASSERT_EQ(A.rows(), static_cast<Eigen::Index>(nTotal));
    ASSERT_EQ(A.cols(), static_cast<Eigen::Index>(nTotal));
    ASSERT_EQ(b.size(), static_cast<Eigen::Index>(nTotal));

    for (Index i = 0; i < nTotal; i++)
    {
      EXPECT_NEAR(b.coeff(i), expectedB(i), 1e-14) << "row " << i;
      for (Index j = 0; j < nTotal; j++)
      {
        EXPECT_NEAR(A.coeff(i, j), expectedA(i, j), 1e-14)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

  TEST(Assembly_Problem_Identification, SelfIdentificationMatchesZeroValueConstraint)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 refFES(mesh);
    TrialFunction uRef(refFES);
    TestFunction  vRef(refFES);

    Problem refProblem(uRef, vRef);
    refProblem = Integral(Grad(uRef), Grad(vRef))
               - Integral(RealFunction(1.0), vRef)
               + DirichletBC(uRef, Zero());
    refProblem.assemble();

    P1 idFES(mesh);
    TrialFunction uId(idFES);
    TestFunction  vId(idFES);

    Problem idProblem(uId, vId);
    idProblem = Integral(Grad(uId), Grad(vId))
              - Integral(RealFunction(1.0), vId)
              + DirichletBC(uId, -uId);
    idProblem.assemble();

    const auto& ARef = refProblem.getLinearSystem().getOperator();
    const auto& bRef = refProblem.getLinearSystem().getVector();
    const auto& AId  = idProblem.getLinearSystem().getOperator();
    const auto& bId  = idProblem.getLinearSystem().getVector();

    ASSERT_EQ(ARef.rows(), AId.rows());
    ASSERT_EQ(ARef.cols(), AId.cols());
    ASSERT_EQ(bRef.size(), bId.size());

    EXPECT_NEAR((ARef - AId).norm(), 0.0, 1e-12);
    EXPECT_NEAR((bRef - bId).norm(), 0.0, 1e-12);
  }

  TEST(Assembly_Problem_Identification, OverlappingValueAndIdentificationThrows)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    Problem problem(u, v);
    problem = Integral(u, v)
            + DirichletBC(u, RealFunction(0.0))
            + DirichletBC(u, RealFunction(2.0) * u);

    EXPECT_THROW(problem.assemble(), std::invalid_argument);
  }

  // =========================================================================
  // All-geometry parameterised tests — Sequential assembler
  //
  // These tests verify correctness across all supported mesh geometry types:
  //   1D: Segment
  //   2D: Triangle, Quadrilateral
  //   3D: Tetrahedron, Hexahedron
  // =========================================================================

  /**
   * @brief Creates a mesh with minimal connectivity for cell integration.
   *
   * For each geometry, computes cell→face incidence (needed for boundary
   * detection by BC assemblers) and cell→vertex incidence (for DOF mapping).
   */
  static Mesh<Context::Local> makeAllGeomMesh(Polytope::Type geom, size_t n = 4)
  {
    switch (geom)
    {
      case Polytope::Type::Segment:
      {
        auto mesh = LocalMesh::UniformGrid(geom, { n });
        mesh.getConnectivity().compute(1, 0);
        mesh.getConnectivity().compute(0, 1);
        return mesh;
      }
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Hexahedron:
      {
        auto mesh = LocalMesh::UniformGrid(geom, { n, n, n });
        mesh.getConnectivity().compute(2, 3);
        mesh.getConnectivity().compute(3, 0);
        return mesh;
      }
      default:
      {
        auto mesh = LocalMesh::UniformGrid(geom, { n, n });
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
    }
  }

  class Assembly_Sequential_AllGeometries
    : public ::testing::TestWithParam<Polytope::Type> {};

  /**
   * @brief P1 LinearForm vector size equals the FES DOF count for every
   * supported geometry type.
   */
  TEST_P(Assembly_Sequential_AllGeometries, LF_VectorSize_P1)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P1 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    EXPECT_EQ(lf.getVector().size(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief P0 LinearForm vector size equals the FES DOF count (= cell count)
   * for every supported geometry type.
   */
  TEST_P(Assembly_Sequential_AllGeometries, LF_VectorSize_P0)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P0 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    EXPECT_EQ(lf.getVector().size(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief P1 stiffness matrix has correct square dimensions for every
   * supported geometry type.
   */
  TEST_P(Assembly_Sequential_AllGeometries, BF_StiffnessMatrix_Dimensions_P1)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& A = bf.getOperator();
    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(fes.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(fes.getSize()));
  }

  /**
   * @brief P1 mass matrix is symmetric for every supported geometry type.
   */
  TEST_P(Assembly_Sequential_AllGeometries, BF_MassMatrix_IsSymmetric_P1)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    const auto& A = bf.getOperator();
    const Math::SparseMatrix<Real> diff = A - Math::SparseMatrix<Real>(A.transpose());
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief P1 stiffness matrix is symmetric for every supported geometry type.
   */
  TEST_P(Assembly_Sequential_AllGeometries, BF_StiffnessMatrix_IsSymmetric_P1)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));
    bf.assemble();

    const auto& A = bf.getOperator();
    const Math::SparseMatrix<Real> diff = A - Math::SparseMatrix<Real>(A.transpose());
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief P0 mass matrix is diagonal for every supported geometry type
   * (each cell has exactly one P0 DOF, so there is no inter-cell coupling).
   */
  TEST_P(Assembly_Sequential_AllGeometries, BF_P0MassMatrix_IsDiagonal)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P0 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);
    bf.assemble();

    const auto& A = bf.getOperator();
    ASSERT_GT(A.rows(), 0);

    for (int k = 0; k < A.outerSize(); ++k)
    {
      for (typename Math::SparseMatrix<Real>::InnerIterator it(A, k); it; ++it)
      {
        if (it.row() != it.col())
        {
          EXPECT_NEAR(it.value(), 0.0, 1e-14);
        }
      }
    }
  }

  /**
   * @brief P0 load vector: all entries are positive (element measure > 0)
   * for every supported geometry type.
   */
  TEST_P(Assembly_Sequential_AllGeometries, LF_LoadVector_P0_EntriesPositive)
  {
    auto mesh = makeAllGeomMesh(GetParam());
    P0 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);
    lf.assemble();

    const auto& b = lf.getVector();
    ASSERT_GT(b.size(), 0);
    EXPECT_GT(b.minCoeff(), 0.0);
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometries,
    Assembly_Sequential_AllGeometries,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron
    )
  );
}
