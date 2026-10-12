/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  template <class Scalar = Real, class Operator = Math::SparseMatrix<Scalar>>
  void checkNonlocalBlocks(
    bool restrictTrial = false, bool restrictTest = false, bool missingAttribute = false)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.setAttribute({1, 0}, 7);
    mesh.setAttribute({1, 1}, 8);
    P1<Scalar> fes(mesh);
    TrialFunction u(fes);
    TrialFunction p(fes);
    TestFunction v(fes);
    TestFunction q(fes);
    using System = Math::LinearSystem<Operator, Math::Vector<Scalar>>;
    Problem<System, decltype(u), decltype(v), decltype(p), decltype(q)> problem(
      u, v, p, q);
    const Scalar factor = [] {
      if constexpr (std::is_same_v<Scalar, Complex>)
        return Scalar(1, -2);
      else
        return Scalar(1);
    }();
    const Scalar coefficients[2][2] = {
      {Real(2) * factor, Real(-3) * factor}, {Real(0.5) * factor, Real(4) * factor}};
    auto kernel = [](Scalar value) {
      return [value](const Point&, const Point&) { return value; };
    };
    const auto k00 = kernel(coefficients[0][0]);
    const auto k01 = kernel(coefficients[0][1]);
    const auto k10 = kernel(coefficients[1][0]);
    const auto k11 = kernel(coefficients[1][1]);
    typename decltype(problem)::ProblemBodyType body;
    body.getGlobalBFIs().add(Integral(Potential(k00, u), v));
    body.getGlobalBFIs().add(Integral(Potential(k01, p), v));
    body.getGlobalBFIs().add(Integral(Potential(k10, u), q));
    body.getGlobalBFIs().add(Integral(Potential(k11, p), q));
    for (auto& integral : body.getGlobalBFIs())
    {
      if (restrictTrial)
        integral.setTrialAttributes({missingAttribute ? 999u : 7u});
      if (restrictTest)
        integral.setTestAttributes({missingAttribute ? 999u : 8u});
    }
    problem = body;

    // The two unit-length cells give exact nodal integrals (1/2, 1, 1/2).
    // A constant kernel separates the double integral into their outer product.
    const Real trialMoments[] = {
      0.5, restrictTrial ? 0.5 : 1.0, restrictTrial ? 0.0 : 0.5};
    const Real testMoments[] = {restrictTest ? 0.0 : 0.5, restrictTest ? 0.5 : 1.0, 0.5};
    auto check = [&](auto& instance, size_t stage) {
      instance.assemble();
      const auto& matrix = instance.getLinearSystem().getOperator();
      ASSERT_EQ(matrix.rows(), 6);
      ASSERT_EQ(matrix.cols(), 6);
      for (size_t row = 0; row < 6; ++row)
      {
        for (size_t col = 0; col < 6; ++col)
        {
          const Scalar expected = missingAttribute ? Scalar(0)
                                                   : coefficients[row / 3][col / 3] *
              testMoments[row % 3] * trialMoments[col % 3];
          EXPECT_NEAR(std::abs(matrix.coeff(row, col) - expected), 0.0, 1e-12)
            << "row=" << row << ", col=" << col << ", assembly=" << stage;
        }
      }
    };
    for (size_t stage = 0; stage < 2; ++stage)
      check(problem, stage);
    auto copied = problem;
    check(copied, 2);
    auto moved = std::move(copied);
    check(moved, 3);
    auto assigned = problem;
    assigned = problem;
    check(assigned, 4);
    auto moveSource = problem;
    auto moveAssigned = problem;
    moveAssigned = std::move(moveSource);
    check(moveAssigned, 5);
    std::unique_ptr<decltype(problem)> clone(problem.copy());
    check(*clone, 6);
  }

  /// @brief Checks every sparse block against exact double integrals.
  TEST(NonlocalBlock, SparseReal)
  {
    checkNonlocalBlocks();
  }

  /// @brief Trial-region restriction changes the column moments of every block.
  TEST(NonlocalBlock, RestrictedTrial)
  {
    checkNonlocalBlocks(true, false);
  }

  /// @brief Test-region restriction changes the row moments of every block.
  TEST(NonlocalBlock, RestrictedTest)
  {
    checkNonlocalBlocks(false, true);
  }

  /// @brief Independent region restrictions produce the selected cross-cell coupling.
  TEST(NonlocalBlock, RestrictedBoth)
  {
    checkNonlocalBlocks(true, true);
  }

  /// @brief Nonmatching attributes produce zero entries without retaining old assembly data.
  TEST(NonlocalBlock, MissingAttributes)
  {
    checkNonlocalBlocks(true, true, true);
  }

  /// @brief Dense real mixed assembly uses every nonlocal block and both independent regions.
  TEST(NonlocalBlock, DenseReal)
  {
    checkNonlocalBlocks<Real, Math::Matrix<Real>>();
    checkNonlocalBlocks<Real, Math::Matrix<Real>>(true, true);
    checkNonlocalBlocks<Real, Math::Matrix<Real>>(true, true, true);
  }

  /// @brief Dense complex mixed assembly preserves the kernel's nonreal coefficient.
  TEST(NonlocalBlock, DenseComplex)
  {
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>>();
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>>(true, true);
    checkNonlocalBlocks<Complex, Math::Matrix<Complex>>(true, true, true);
  }

  /// @brief Sparse complex mixed assembly agrees with exact separated double integrals.
  TEST(NonlocalBlock, SparseComplex)
  {
    checkNonlocalBlocks<Complex>();
    checkNonlocalBlocks<Complex>(true, true);
    checkNonlocalBlocks<Complex>(true, true, true);
  }

}
