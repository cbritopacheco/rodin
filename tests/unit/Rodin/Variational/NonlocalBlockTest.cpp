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
  void checkNonlocalBlocks(
    bool restrictTrial = false, bool restrictTest = false, bool missingAttribute = false)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.setAttribute({1, 0}, 7);
    mesh.setAttribute({1, 1}, 8);
    P1<Real> fes(mesh);
    TrialFunction u(fes);
    TrialFunction p(fes);
    TestFunction v(fes);
    TestFunction q(fes);
    Problem problem(u, v, p, q);
    const Real coefficients[2][2] = {{2, -3}, {0.5, 4}};
    auto kernel = [](Real value) {
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
    for (size_t assembly = 0; assembly < 2; ++assembly)
    {
      problem.assemble();
      const auto& matrix = problem.getLinearSystem().getOperator();
      ASSERT_EQ(matrix.rows(), 6);
      ASSERT_EQ(matrix.cols(), 6);
      for (size_t row = 0; row < 6; ++row)
      {
        for (size_t col = 0; col < 6; ++col)
        {
          const Real expected = missingAttribute ? Real(0)
                                                 : coefficients[row / 3][col / 3] *
              testMoments[row % 3] * trialMoments[col % 3];
          EXPECT_NEAR(std::abs(matrix.coeff(row, col) - expected), 0.0, 1e-12)
            << "row=" << row << ", col=" << col << ", assembly=" << assembly;
        }
      }
    }
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

}
