/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Unit tests for Assembly::OpenMP.
 *
 * These tests exercise the OpenMP assembler directly — verifying that
 * LinearForm vectors and BilinearForm matrices assembled with OpenMP are
 * numerically identical to those assembled with Sequential for all
 * supported 2-D polygon types (Triangle, Quadrilateral).
 *
 * The file is only compiled when RODIN_USE_OPENMP is defined; the
 * CMakeLists.txt registers the target inside an if(RODIN_USE_OPENMP) block.
 */
#include <gtest/gtest.h>

#include <boost/bimap.hpp>

#include "Rodin/Configure.h"
#include "Rodin/Variational.h"
#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Assembly/OpenMP.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  // =========================================================================
  // Helper
  // =========================================================================

  static Mesh<Context::Local> makeMesh(Polytope::Type geom, size_t n = 4)
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

  void checkOpenMPSelfIdentificationMatchesZeroValueConstraint()
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 refFES(mesh);
    TrialFunction uRef(refFES);
    TestFunction  vRef(refFES);

    auto refBody =
      Integral(Grad(uRef), Grad(vRef))
      - Integral(RealFunction(1.0), vRef)
      + DirichletBC(uRef, Zero());

    P1 idFES(mesh);
    TrialFunction uId(idFES);
    TestFunction  vId(idFES);

    auto idBody =
      Integral(Grad(uId), Grad(vId))
      - Integral(RealFunction(1.0), vId)
      + DirichletBC(uId, -uId);

    Problem refProblem(uRef, vRef);
    refProblem = refBody;
    refProblem.assemble();

    Problem idProblem(uId, vId);
    idProblem = idBody;
    idProblem.assemble();

    const auto matrixDiff =
      refProblem.getLinearSystem().getOperator()
      - idProblem.getLinearSystem().getOperator();
    const auto vectorDiff =
      refProblem.getLinearSystem().getVector()
      - idProblem.getLinearSystem().getVector();

    EXPECT_NEAR(matrixDiff.norm(), 0.0, 1e-12);
    EXPECT_NEAR(vectorDiff.norm(), 0.0, 1e-12);
  }

  // =========================================================================
  // LinearForm — OpenMP matches Sequential across geometries
  // =========================================================================

  class Assembly_OpenMP_LinearForm : public ::testing::TestWithParam<Polytope::Type> {};

  /**
   * @brief Constant load f = 1: OpenMP result equals Sequential result.
   */
  TEST_P(Assembly_OpenMP_LinearForm, ConstantLoad_MatchesSequential_P1)
  {
    auto mesh = makeMesh(GetParam());
    P1 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);

    using LFType = decltype(lf);
    using VecType = Math::Vector<Real>;

    VecType seqResult;
    Assembly::Sequential<VecType, LFType> seqAsm;
    seqAsm.execute(seqResult, { fes, lf.getIntegrators() });

    VecType ompResult;
    Assembly::OpenMP<VecType, LFType> ompAsm;
    ompAsm.execute(ompResult, { fes, lf.getIntegrators() });

    ASSERT_EQ(seqResult.size(), ompResult.size());
    EXPECT_NEAR((seqResult - ompResult).norm(), 0.0, 1e-12);
  }

  /**
   * @brief P0 constant load: OpenMP result equals Sequential result.
   */
  TEST_P(Assembly_OpenMP_LinearForm, ConstantLoad_MatchesSequential_P0)
  {
    auto mesh = makeMesh(GetParam());
    P0 fes(mesh);
    TestFunction v(fes);

    LinearForm lf(v);
    lf = Integral(RealFunction(1.0), v);

    using LFType = decltype(lf);
    using VecType = Math::Vector<Real>;

    VecType seqResult;
    Assembly::Sequential<VecType, LFType> seqAsm;
    seqAsm.execute(seqResult, { fes, lf.getIntegrators() });

    VecType ompResult;
    Assembly::OpenMP<VecType, LFType> ompAsm;
    ompAsm.execute(ompResult, { fes, lf.getIntegrators() });

    ASSERT_EQ(seqResult.size(), ompResult.size());
    EXPECT_NEAR((seqResult - ompResult).norm(), 0.0, 1e-12);
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometries,
    Assembly_OpenMP_LinearForm,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron
    )
  );

  // =========================================================================
  // BilinearForm — OpenMP matches Sequential across geometries
  // =========================================================================

  class Assembly_OpenMP_BilinearForm : public ::testing::TestWithParam<Polytope::Type> {};

  /**
   * @brief P1 stiffness matrix: OpenMP result equals Sequential result.
   */
  TEST_P(Assembly_OpenMP_BilinearForm, P1StiffnessMatrix_MatchesSequential)
  {
    auto mesh = makeMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));

    using BFType = decltype(bf);
    using MatType = Math::SparseMatrix<Real>;

    MatType seqResult;
    Assembly::Sequential<MatType, BFType> seqAsm;
    seqAsm.execute(seqResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    MatType ompResult;
    Assembly::OpenMP<MatType, BFType> ompAsm;
    ompAsm.execute(ompResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    ASSERT_EQ(seqResult.rows(), ompResult.rows());
    ASSERT_EQ(seqResult.cols(), ompResult.cols());
    const MatType diff = seqResult - ompResult;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief P1 mass matrix: OpenMP result equals Sequential result.
   */
  TEST_P(Assembly_OpenMP_BilinearForm, P1MassMatrix_MatchesSequential)
  {
    auto mesh = makeMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);

    using BFType = decltype(bf);
    using MatType = Math::SparseMatrix<Real>;

    MatType seqResult;
    Assembly::Sequential<MatType, BFType> seqAsm;
    seqAsm.execute(seqResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    MatType ompResult;
    Assembly::OpenMP<MatType, BFType> ompAsm;
    ompAsm.execute(ompResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    ASSERT_EQ(seqResult.rows(), ompResult.rows());
    ASSERT_EQ(seqResult.cols(), ompResult.cols());
    const MatType diff = seqResult - ompResult;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief P0 mass matrix: OpenMP result equals Sequential result.
   */
  TEST_P(Assembly_OpenMP_BilinearForm, P0MassMatrix_MatchesSequential)
  {
    auto mesh = makeMesh(GetParam());
    P0 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(u, v);

    using BFType = decltype(bf);
    using MatType = Math::SparseMatrix<Real>;

    MatType seqResult;
    Assembly::Sequential<MatType, BFType> seqAsm;
    seqAsm.execute(seqResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    MatType ompResult;
    Assembly::OpenMP<MatType, BFType> ompAsm;
    ompAsm.execute(ompResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    ASSERT_EQ(seqResult.rows(), ompResult.rows());
    ASSERT_EQ(seqResult.cols(), ompResult.cols());
    const MatType diff = seqResult - ompResult;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief Non-square mixed BF (P0 test × P1 trial): OpenMP equals Sequential.
   */
  TEST_P(Assembly_OpenMP_BilinearForm, NonSquare_P0test_P1trial_MatchesSequential)
  {
    auto mesh = makeMesh(GetParam());
    P1 p1h(mesh);
    P0 p0h(mesh);
    TrialFunction u(p1h);
    TestFunction  q(p0h);

    BilinearForm bf(u, q);
    bf = Integral(u, q);

    using BFType = decltype(bf);
    using MatType = Math::SparseMatrix<Real>;

    MatType seqResult;
    Assembly::Sequential<MatType, BFType> seqAsm;
    seqAsm.execute(seqResult, {
        p1h, p0h, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    MatType ompResult;
    Assembly::OpenMP<MatType, BFType> ompAsm;
    ompAsm.execute(ompResult, {
        p1h, p0h, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    ASSERT_EQ(seqResult.rows(), ompResult.rows());
    ASSERT_EQ(seqResult.cols(), ompResult.cols());
    const MatType diff = seqResult - ompResult;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  /**
   * @brief OpenMP with explicit thread count 1 gives same result as Sequential.
   */
  TEST_P(Assembly_OpenMP_BilinearForm, SingleThread_EqualsSequential_P1Stiffness)
  {
    auto mesh = makeMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    BilinearForm bf(u, v);
    bf = Integral(Grad(u), Grad(v));

    using BFType = decltype(bf);
    using MatType = Math::SparseMatrix<Real>;

    MatType seqResult;
    Assembly::Sequential<MatType, BFType> seqAsm;
    seqAsm.execute(seqResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    MatType ompResult;
    Assembly::OpenMP<MatType, BFType> ompAsm(1);
    ompAsm.execute(ompResult, {
        fes, fes, bf.getLocalIntegrators(), bf.getGlobalIntegrators() });

    ASSERT_EQ(seqResult.rows(), ompResult.rows());
    const MatType diff = seqResult - ompResult;
    EXPECT_NEAR(diff.norm(), 0.0, 1e-12);
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometries,
    Assembly_OpenMP_BilinearForm,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron
    )
  );

  // =========================================================================
  // Problem (multi-variable) — OpenMP matches Sequential
  // =========================================================================

  /**
   * @brief Single-variable Problem: Default assembler (OpenMP when available)
   * produces same operator as Sequential for all 2D geometries.
   */
  class Assembly_OpenMP_Problem : public ::testing::TestWithParam<Polytope::Type> {};

  TEST_P(Assembly_OpenMP_Problem, SingleVar_P1_DimensionsCorrect)
  {
    auto mesh = makeMesh(GetParam());
    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    Problem problem(u, v);
    problem = Integral(Grad(u), Grad(v))
            - Integral(RealFunction(1.0), v);
    problem.assemble();

    const auto& A = problem.getLinearSystem().getOperator();
    const auto& b = problem.getLinearSystem().getVector();

    EXPECT_EQ(A.rows(), static_cast<Eigen::Index>(fes.getSize()));
    EXPECT_EQ(A.cols(), static_cast<Eigen::Index>(fes.getSize()));
    EXPECT_EQ(b.size(), static_cast<Eigen::Index>(fes.getSize()));
  }

  TEST_P(Assembly_OpenMP_Problem, MultiVar_P0P1_DimensionsCorrect)
  {
    auto mesh = makeMesh(GetParam());

    P0 p0h(mesh);
    P1 p1h(mesh);

    TrialFunction u(p1h);
    TrialFunction p(p0h);
    TestFunction  v(p1h);
    TestFunction  q(p0h);

    Problem mixed(u, v, p, q);
    mixed = Integral(u, v) + Integral(p, q);
    mixed.assemble();

    const auto& A = mixed.getLinearSystem().getOperator();
    const Eigen::Index expected =
      static_cast<Eigen::Index>(p1h.getSize() + p0h.getSize());

    EXPECT_EQ(A.rows(), expected);
    EXPECT_EQ(A.cols(), expected);
  }

  TEST(Assembly_OpenMP_Problem, IdentificationVectorMasterMatchesSequential)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
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

    BilinearForm uu(u, v);
    uu.getOperator().resize(nSlave, nSlave);
    std::vector<Eigen::Triplet<Real>> triplets{
      Eigen::Triplet<Real>(0, 0, 2.0),
      Eigen::Triplet<Real>(0, 1, 3.0),
      Eigen::Triplet<Real>(1, 0, 5.0),
      Eigen::Triplet<Real>(1, 1, 7.0)
    };
    uu.getOperator().setFromTriplets(triplets.begin(), triplets.end());

    LinearForm loadU(v);
    loadU.getVector().resize(nSlave);
    loadU.getVector().setZero();
    loadU.getVector().coeffRef(0) = 11.0;
    loadU.getVector().coeffRef(1) = -13.0;

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

    LinearSystemType seqLS;
    LinearSystemType ompLS;
    Assembly::Sequential<LinearSystemType, ProblemType> seqAsm;
    Assembly::OpenMP<LinearSystemType, ProblemType> ompAsm;
    seqAsm.execute(seqLS, input);
    ompAsm.execute(ompLS, input);

    const auto& ASeq = seqLS.getOperator();
    const auto& AOmp = ompLS.getOperator();
    const auto& bSeq = seqLS.getVector();
    const auto& bOmp = ompLS.getVector();

    ASSERT_EQ(ASeq.rows(), static_cast<Eigen::Index>(nTotal));
    ASSERT_EQ(AOmp.rows(), ASeq.rows());
    ASSERT_EQ(AOmp.cols(), ASeq.cols());
    ASSERT_EQ(bOmp.size(), bSeq.size());

    for (Index i = 0; i < nTotal; i++)
    {
      EXPECT_NEAR(bOmp.coeff(i), bSeq.coeff(i), 1e-14) << "row " << i;
      for (Index j = 0; j < nTotal; j++)
      {
        EXPECT_NEAR(AOmp.coeff(i, j), ASeq.coeff(i, j), 1e-14)
          << "entry (" << i << ", " << j << ")";
      }
    }
  }

  TEST(Assembly_OpenMP_Problem, SelfIdentificationMatchesZeroValueConstraint)
  {
    checkOpenMPSelfIdentificationMatchesZeroValueConstraint();
  }

  INSTANTIATE_TEST_SUITE_P(
    AllGeometries,
    Assembly_OpenMP_Problem,
    ::testing::Values(
      Polytope::Type::Segment,
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron
    )
  );
}
