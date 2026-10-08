/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_HINGEPROBLEM_H
#define RODIN_ADAPTATION_WNGIR_HINGEPROBLEM_H

#include <chrono>
#include <type_traits>
#include <utility>
#include "Rodin/Assembly.h"
#include "Rodin/Variational/Problem.h"

namespace Rodin::Adaptation::WNGIR
{
  /**
   * @brief Native variational tangent problem for the affine quadratic hinges.
   *
   * @par Architecture
   * The assigned forms assemble H and b such that the absolute linearized
   * minimizer satisfies H w = b. After symmetry and homogeneous boundary
   * elimination, assembly converts the load to b - H v, so NewtonSolver solves
   * the correction equation H delta = b - H v. The current increment is a
   * referenced GridFunction; its lifetime must include the Newton solve.
   */
  template <class TrialFunction, class TestFunction>
  class HingeProblem final
    : public std::decay_t<decltype(Variational::Problem(
        std::declval<TrialFunction&>(), std::declval<TestFunction&>()))>
  {
    public:
      using Parent = std::decay_t<decltype(Variational::Problem(
        std::declval<TrialFunction&>(), std::declval<TestFunction&>()))>;
      using Displacement = typename Parent::SolutionType;
      using Parent::operator=;

      HingeProblem(TrialFunction& trial, TestFunction& test)
        : Parent(trial, test)
      {}

      HingeProblem& setState(const Displacement& state)
      {
        m_state = &state;
        return *this;
      }

      HingeProblem& setBoundaryDOFs(const IndexMap<Real>& dofs)
      {
        m_boundaryDOFs = &dofs;
        return *this;
      }

      HingeProblem& setCentering(const Math::Matrix<Real>& couplings)
      {
        m_centering = couplings;
        return *this;
      }

      HingeProblem& assemble() override
      {
        assert(m_state);
        const auto start = std::chrono::steady_clock::now();
        Parent::assemble();
        auto& system = this->getLinearSystem();
        auto& matrix = system.getOperator();
        matrix =
          (Real(0.5) * (matrix + Math::SparseMatrix<Real>(matrix.transpose()))).eval();
        if (m_boundaryDOFs && !m_boundaryDOFs->empty())
          system.eliminate(*m_boundaryDOFs);
        matrix.makeCompressed();
        system.getVector() -= matrix * m_state->getData();
        if (m_centering.size() != 0)
          system.getVector() += m_centering * (m_centering.transpose() * m_state->getData());
        m_assemblySeconds =
          std::chrono::duration<Real>(std::chrono::steady_clock::now() - start).count();
        return *this;
      }

      Real getAssemblySeconds() const noexcept
      {
        return m_assemblySeconds;
      }

      HingeProblem* copy() const noexcept override
      {
        return new HingeProblem(*this);
      }

    private:
      const Displacement* m_state = nullptr;
      const IndexMap<Real>* m_boundaryDOFs = nullptr;
      Math::Matrix<Real> m_centering;
      Real m_assemblySeconds = 0;
  };
}
#endif
