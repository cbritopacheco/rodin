#include <gtest/gtest.h>

#include "Rodin/Math.h"

using namespace Rodin;

class UmbrellaIncludeTest : public ::testing::Test
{};

TEST_F(UmbrellaIncludeTest, UmbrellaExposesExpectedMathAPIs)
{
  Math::Vector<Real> v(2);
  v << 1.0, 2.0;
  Math::Matrix<Real> m(2, 2);
  m.setIdentity();
  Math::SparseMatrix<Real> s(2, 2);
  s.setIdentity();

  Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>> lsDense;
  Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>> lsSparse;
  lsDense.getOperator().resize(2, 2);
  lsSparse.getOperator().resize(2, 2);

  Math::Rad angle(Math::Constants::pi());
  EXPECT_NEAR(static_cast<Real>(angle), Math::Constants::pi(), 1e-12);

  Math::RungeKutta::RK2 rk2;
  Math::RungeKutta::RK4 rk4;
  (void) rk2;
  (void) rk4;

  Math::RootFinding::NewtonRaphson<Real> solver;
  auto f = [](Real x)
  {
    return std::make_pair(x * x - 2.0, 2.0 * x);
  };
  const auto root = solver.solve(f, 1.5, 1.0, 2.0);
  ASSERT_TRUE(root.has_value());
}
