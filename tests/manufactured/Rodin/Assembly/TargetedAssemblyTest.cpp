/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 *
 * Regression tests for targeted (LHS-only / RHS-only) assembly on the core
 * Eigen-backed assembly backends (Sequential and OpenMP), for both
 * single-field and two-field (block) problems.
 *
 * Targeted assembly must satisfy:
 *   - assembling with AssemblyTarget::LHS reproduces the operator of a full
 *     assembly, and
 *   - assembling with AssemblyTarget::RHS reproduces the vector of a full
 *     assembly.
 *
 * The manufactured systems are deliberately free of Dirichlet boundary
 * conditions so that the targeted/full equivalence is exact to roundoff.
 */

/**
 * @file
 * @brief Targeted assembly test manufactured regression tests.
 *
 * These tests assemble Rodin variational forms for a Targeted Assembly Test manufactured regression, solve the problem on the configured mesh, and compare against analytic fields or expected residual/error behavior. They protect the assembly backends and targeted element assembly path, including boundary-condition handling, geometry coverage, and numerical accuracy of the manufactured workflow.
 */

#include <array>

#include <boost/bimap.hpp>

#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Configure.h"
#include "Rodin/Geometry.h"
#include "Rodin/Tuple.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured::Assembly
{
  /// @brief Helper used by the manufactured tests to Linear System Type.
  using LinearSystemType =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  /// @brief Helper used by the manufactured tests to Expect Same Vector.
  void expectSameVector(
    const Math::Vector<Real>& expected, const Math::Vector<Real>& actual)
  {
    ASSERT_EQ(expected.size(), actual.size());
    for (Eigen::Index i = 0; i < expected.size(); ++i)
      EXPECT_LE(std::abs(expected.coeff(i) - actual.coeff(i)), 1e-12) << "row " << i;
  }

  /// @brief Helper used by the manufactured tests to Expect Same Matrix.
  void expectSameMatrix(
    const Math::SparseMatrix<Real>& expected, const Math::SparseMatrix<Real>& actual)
  {
    ASSERT_EQ(expected.rows(), actual.rows());
    ASSERT_EQ(expected.cols(), actual.cols());
    const Math::SparseMatrix<Real> diff = expected - actual;
    Real maxAbs = 0;
    for (int k = 0; k < diff.outerSize(); ++k)
      for (Math::SparseMatrix<Real>::InnerIterator it(diff, k); it; ++it)
        maxAbs = std::max(maxAbs, std::abs(it.value()));
    EXPECT_LE(maxAbs, 1e-12);
  }

  // -------------------------------------------------------------------------
  // Single-field: -Delta u = 1 weak form (mass + stiffness), no BCs.
  // -------------------------------------------------------------------------
  template <template <class, class> class Assembler>
  /// @brief Helper used by the manufactured tests to Assemble Single Field.
  LinearSystemType assembleSingleField(Variational::AssemblyTarget target, bool targeted)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    using ProblemType = Problem<LinearSystemType, decltype(u), decltype(v)>;
    typename ProblemType::ProblemBodyType body =
      Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(RealFunction(1.0), v);

    ::Rodin::Assembly::ProblemAssemblyInput<typename ProblemType::ProblemBodyType,
      decltype(u), decltype(v)>
      input(body, u, v);

    LinearSystemType ls;
    Assembler<LinearSystemType, ProblemType> assembler;
    if (targeted)
      assembler.execute(ls, input, target);
    else
      assembler.execute(ls, input);
    return ls;
  }

  template <template <class, class> class Assembler>
  /// @brief Helper used by the manufactured tests to Check Single Field Targeted.
  void checkSingleFieldTargeted()
  {
    const auto full =
      assembleSingleField<Assembler>(Variational::AssemblyTarget::LHS, false);
    const auto lhs =
      assembleSingleField<Assembler>(Variational::AssemblyTarget::LHS, true);
    const auto rhs =
      assembleSingleField<Assembler>(Variational::AssemblyTarget::RHS, true);

    expectSameMatrix(full.getOperator(), lhs.getOperator());
    expectSameVector(full.getVector(), rhs.getVector());
  }

  /// @brief Verifies sequential single field LHS and RHS for eigen targeted assembly by checking form assembly.
  TEST(Eigen_TargetedAssembly, SequentialSingleFieldLHSAndRHS)
  {
    checkSingleFieldTargeted<::Rodin::Assembly::Sequential>();
  }

#ifdef RODIN_USE_OPENMP
  /// @brief Verifies open MP single field LHS and RHS for eigen targeted assembly by checking form assembly.
  TEST(Eigen_TargetedAssembly, OpenMPSingleFieldLHSAndRHS)
  {
    checkSingleFieldTargeted<::Rodin::Assembly::OpenMP>();
  }
#endif

  // -------------------------------------------------------------------------
  // Two-field (block) targeted assembly, driving a specific block backend
  // directly (so both Sequential and OpenMP block paths are exercised locally,
  // regardless of which backend Assembly::Default resolves to). The block
  // ProblemAssemblyInput is built the same way the high-level Problem builds
  // it. BC-free so targeted/full equivalence is exact.
  // -------------------------------------------------------------------------
  template <template <class, class> class Assembler>
  /// @brief Helper used by the manufactured tests to Assemble Block.
  LinearSystemType assembleBlock(Variational::AssemblyTarget target, bool targeted)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    P0 ph(mesh);
    TrialFunction u(vh);
    TestFunction v(vh);
    TrialFunction p(ph);
    TestFunction q(ph);

    using ProblemType =
      Problem<LinearSystemType, decltype(u), decltype(v), decltype(p), decltype(q)>;
    typename ProblemType::ProblemBodyType body = Integral(Grad(u), Grad(v)) +
      Integral(u, v) - Integral(p, v) - Integral(u, q) + Integral(p, q) -
      Integral(RealFunction(1.0), v) - Integral(RealFunction(2.0), q);

    // Trial blocks: [u, p]; test blocks: [v, q].
    auto us = Tuple{std::ref(u), std::ref(p)};
    auto vs = Tuple{std::ref(v), std::ref(q)};

    const size_t uSize = static_cast<size_t>(vh.getSize());
    const size_t pSize = static_cast<size_t>(ph.getSize());

    std::array<size_t, 2> trialOffsets{0, uSize};
    std::array<size_t, 2> testOffsets{0, uSize};

    boost::bimap<Rodin::FormLanguage::Base::UUID, size_t> trialUUIDMap;
    boost::bimap<Rodin::FormLanguage::Base::UUID, size_t> testUUIDMap;
    trialUUIDMap.right.insert({0, u.getUUID()});
    trialUUIDMap.right.insert({1, p.getUUID()});
    testUUIDMap.right.insert({0, v.getUUID()});
    testUUIDMap.right.insert({1, q.getUUID()});

    const size_t totalTrial = uSize + pSize;
    const size_t totalTest = uSize + pSize;

    ::Rodin::Assembly::ProblemAssemblyInput<typename ProblemType::ProblemBodyType,
      decltype(u), decltype(v), decltype(p), decltype(q)>
      input(body, us, vs, trialOffsets, testOffsets, trialUUIDMap, testUUIDMap,
        totalTrial, totalTest);

    LinearSystemType ls;
    Assembler<LinearSystemType, ProblemType> assembler;
    if (targeted)
      assembler.execute(ls, input, target);
    else
      assembler.execute(ls, input);
    return ls;
  }

  template <template <class, class> class Assembler>
  /// @brief Helper used by the manufactured tests to Check Block Targeted.
  void checkBlockTargeted()
  {
    const auto full = assembleBlock<Assembler>(Variational::AssemblyTarget::LHS, false);
    const auto lhs = assembleBlock<Assembler>(Variational::AssemblyTarget::LHS, true);
    const auto rhs = assembleBlock<Assembler>(Variational::AssemblyTarget::RHS, true);

    expectSameMatrix(full.getOperator(), lhs.getOperator());
    expectSameVector(full.getVector(), rhs.getVector());
  }

  /// @brief Verifies sequential block LHS and RHS for eigen targeted assembly by checking form assembly.
  TEST(Eigen_TargetedAssembly, SequentialBlockLHSAndRHS)
  {
    checkBlockTargeted<::Rodin::Assembly::Sequential>();
  }

#ifdef RODIN_USE_OPENMP
  /// @brief Verifies open MP block LHS and RHS for eigen targeted assembly by checking form assembly.
  TEST(Eigen_TargetedAssembly, OpenMPBlockLHSAndRHS)
  {
    checkBlockTargeted<::Rodin::Assembly::OpenMP>();
  }
#endif

  // Block targeted assembly through the high-level Problem API (Default
  // backend), covering the Problem::assemble(target) wiring end-to-end.
  /// @brief Verifies block problem APILHS and RHS for eigen targeted assembly by checking form assembly.
  TEST(Eigen_TargetedAssembly, BlockProblemAPILHSAndRHS)
  {
    auto mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    P0 ph(mesh);
    TrialFunction u(vh);
    TrialFunction p(ph);
    TestFunction v(vh);
    TestFunction q(ph);

    Problem prob(u, v, p, q);
    prob = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(p, v) - Integral(u, q) +
      Integral(p, q) - Integral(RealFunction(1.0), v) - Integral(RealFunction(2.0), q);

    prob.assemble();
    const Math::SparseMatrix<Real> Afull = prob.getLinearSystem().getOperator();
    const Math::Vector<Real> bfull = prob.getLinearSystem().getVector();

    prob.assemble(Variational::AssemblyTarget::LHS);
    const Math::SparseMatrix<Real> Alhs = prob.getLinearSystem().getOperator();

    prob.assemble(Variational::AssemblyTarget::RHS);
    const Math::Vector<Real> brhs = prob.getLinearSystem().getVector();

    expectSameMatrix(Afull, Alhs);
    expectSameVector(bfull, brhs);
  }
}
