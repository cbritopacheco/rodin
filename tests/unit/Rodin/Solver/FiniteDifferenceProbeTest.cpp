/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <cmath>

#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Solver/FiniteDifferenceProbe.h"
#include "Rodin/Solver/LinearSolver.h"
#include "Rodin/Variational/Problem.h"

using namespace Rodin;

namespace
{
  using DenseLinearSystem = Math::LinearSystem<Math::Matrix<Real>, Math::Vector<Real>>;

  class TwoVarProblem final : public Variational::ProblemBase<DenseLinearSystem>
  {
    public:
      using Parent = Variational::ProblemBase<DenseLinearSystem>;
      using ProblemBodyType = typename Parent::ProblemBodyType;

      explicit TwoVarProblem(Math::Vector<Real>& state)
        : m_state(state)
      {
        m_system.getOperator().resize(2, 2);
        m_system.getVector().resize(2);
        m_system.getSolution().resize(2);
      }

      Parent& operator=(const ProblemBodyType&) override
      {
        return *this;
      }

      void solve(Solver::LinearSolverBase<DenseLinearSystem>&) override {}

      TwoVarProblem& assemble() override
      {
        const Real x = m_state.get()(0);
        const Real y = m_state.get()(1);

        // Residuals:
        // R0 = x^2 + 3y - 1
        // R1 = sin(x) + y^3
        // Tangent in Newton form J d = -R => RHS = -R
        m_system.getVector()(0) = -(x * x + 3.0 * y - 1.0);
        m_system.getVector()(1) = -(std::sin(x) + y * y * y);

        m_system.getOperator()(0, 0) = 2.0 * x;
        m_system.getOperator()(0, 1) = 3.0;
        m_system.getOperator()(1, 0) = std::cos(x);
        m_system.getOperator()(1, 1) = 3.0 * y * y;

        return *this;
      }

      DenseLinearSystem& getLinearSystem() override { return m_system; }
      const DenseLinearSystem& getLinearSystem() const override { return m_system; }

      TwoVarProblem* copy() const noexcept override
      {
        return new TwoVarProblem(*this);
      }

    private:
      std::reference_wrapper<Math::Vector<Real>> m_state;
      DenseLinearSystem m_system;
  };
}

TEST(FiniteDifferenceProbeTest, RecoversResidualJacobian)
{
  Math::Vector<Real> x(2);
  x << 0.4, -0.2;

  auto residual = [](const Math::Vector<Real>& state, Math::Vector<Real>& out)
  {
    out.resize(2);
    out(0) = state(0) * state(0) + 3.0 * state(1) - 1.0;
    out(1) = std::sin(state(0)) + state(1) * state(1) * state(1);
  };

  const auto Jfd = Solver::FiniteDifferenceProbe::jacobian(x, residual);

  Math::Matrix<Real> Janalytic(2, 2);
  Janalytic << 2.0 * x(0), 3.0,
               std::cos(x(0)), 3.0 * x(1) * x(1);

  EXPECT_LT((Jfd - Janalytic).norm(), 1e-7);
}

TEST(FiniteDifferenceProbeTest, ProbesTangentialProblemInNewtonForm)
{
  Math::Vector<Real> state(2);
  state << 0.3, 0.2;

  TwoVarProblem pb(state);
  pb.assemble();
  const auto Janalytic = pb.getLinearSystem().getOperator();

  const auto Jfd = Solver::FiniteDifferenceProbe::tangentialJacobian(
    pb,
    state,
    [&](const Math::Vector<Real>& s) { state = s; });

  const auto rel = Solver::FiniteDifferenceProbe::relativeError(Janalytic, Jfd);
  EXPECT_LT(rel, 1e-7);
}
