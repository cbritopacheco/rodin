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
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Configure.h"
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Manufactured::Assembly
{
  using LinearSystemType =
    Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;

  void expectSameVector(const Math::Vector<Real>& expected,
                        const Math::Vector<Real>& actual)
  {
    ASSERT_EQ(expected.size(), actual.size());
    for (Eigen::Index i = 0; i < expected.size(); ++i)
      EXPECT_LE(std::abs(expected.coeff(i) - actual.coeff(i)), 1e-12)
        << "row " << i;
  }

  void expectSameMatrix(const Math::SparseMatrix<Real>& expected,
                        const Math::SparseMatrix<Real>& actual)
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
  LinearSystemType assembleSingleField(
      Variational::AssemblyTarget target, bool targeted)
  {
    auto mesh =
      Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 fes(mesh);
    TrialFunction u(fes);
    TestFunction  v(fes);

    using ProblemType = Problem<LinearSystemType, decltype(u), decltype(v)>;
    typename ProblemType::ProblemBodyType body =
      Integral(Grad(u), Grad(v)) + Integral(u, v)
        - Integral(RealFunction(1.0), v);

    ::Rodin::Assembly::ProblemAssemblyInput<
      typename ProblemType::ProblemBodyType,
      decltype(u), decltype(v)> input(body, u, v);

    LinearSystemType ls;
    Assembler<LinearSystemType, ProblemType> assembler;
    if (targeted)
      assembler.execute(ls, input, target);
    else
      assembler.execute(ls, input);
    return ls;
  }

  template <template <class, class> class Assembler>
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

  TEST(Eigen_TargetedAssembly, SequentialSingleFieldLHSAndRHS)
  {
    checkSingleFieldTargeted<::Rodin::Assembly::Sequential>();
  }

#ifdef RODIN_USE_OPENMP
  TEST(Eigen_TargetedAssembly, OpenMPSingleFieldLHSAndRHS)
  {
    checkSingleFieldTargeted<::Rodin::Assembly::OpenMP>();
  }
#endif

  // -------------------------------------------------------------------------
  // Two-field (block) targeted assembly through the high-level Problem API.
  // The Problem routes to the Default backend (Sequential or OpenMP depending
  // on the build), exercising the block backend + Problem::assemble(target)
  // wiring end-to-end. BC-free so the targeted/full equivalence is exact.
  // -------------------------------------------------------------------------
  TEST(Eigen_TargetedAssembly, BlockProblemLHSAndRHS)
  {
    auto mesh =
      Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    P1 vh(mesh);
    P0 ph(mesh);
    TrialFunction u(vh);
    TrialFunction p(ph);
    TestFunction  v(vh);
    TestFunction  q(ph);

    Problem prob(u, v, p, q);
    prob = Integral(Grad(u), Grad(v)) + Integral(u, v)
         - Integral(p, v) - Integral(u, q) + Integral(p, q)
         - Integral(RealFunction(1.0), v) - Integral(RealFunction(2.0), q);

    prob.assemble();
    const Math::SparseMatrix<Real> Afull = prob.getLinearSystem().getOperator();
    const Math::Vector<Real>       bfull = prob.getLinearSystem().getVector();

    prob.assemble(Variational::AssemblyTarget::LHS);
    const Math::SparseMatrix<Real> Alhs = prob.getLinearSystem().getOperator();

    prob.assemble(Variational::AssemblyTarget::RHS);
    const Math::Vector<Real>       brhs = prob.getLinearSystem().getVector();

    expectSameMatrix(Afull, Alhs);
    expectSameVector(bfull, brhs);
  }
}
