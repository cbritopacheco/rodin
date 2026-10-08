/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_HINGEPROBLEM_H
#define RODIN_ADAPTATION_SWIFT_HINGEPROBLEM_H

#include <chrono>
#include <type_traits>
#include <utility>
#include "Rodin/Assembly.h"
#include "Rodin/Variational/Problem.h"

namespace Rodin::Adaptation::SWIFT
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
      /// @brief Native variational tangent problem type.
      using Parent = std::decay_t<decltype(Variational::Problem(
        std::declval<TrialFunction&>(), std::declval<TestFunction&>()))>;
      /// @brief Solution field holding the current increment.
      using Displacement = typename Parent::SolutionType;
      using Parent::operator=;

      /**
       * @brief Constructs a tangent problem on the supplied fields.
       * @param trial Increment trial function.
       * @param test Increment test function.
       */
      HingeProblem(TrialFunction& trial, TestFunction& test)
        : Parent(trial, test)
      {}

      /**
       * @brief Binds the current Newton iterate without taking ownership.
       * @param state Increment field whose lifetime includes the Newton solve.
       * @returns This problem.
       */
      HingeProblem& setState(const Displacement& state)
      {
        m_state = &state;
        return *this;
      }

      /**
       * @brief Binds the homogeneous increment boundary conditions.
       * @param dofs Boundary degree-of-freedom map that must outlive assembly.
       * @returns This problem.
       */
      HingeProblem& setBoundaryDOFs(const IndexMap<Real>& dofs)
      {
        m_boundaryDOFs = &dofs;
        return *this;
      }

      /**
       * @brief Copies the global mean-strain couplings.
       * @param couplings Low-rank factor of the centered metric subtraction.
       * @returns This problem.
       */
      HingeProblem& setCentering(const Math::Matrix<Real>& couplings)
      {
        m_centering = couplings;
        return *this;
      }

      /**
       * @brief Assembles the tangent and converts its load to a Newton residual.
       * @returns This assembled problem.
       */
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
          system.getVector() +=
            m_centering * (m_centering.transpose() * m_state->getData());
        m_assemblySeconds =
          std::chrono::duration<Real>(std::chrono::steady_clock::now() - start).count();
        return *this;
      }

      /**
       * @brief Returns the duration of the most recent assembly.
       * @returns Elapsed wall-clock seconds.
       */
      Real getAssemblySeconds() const noexcept
      {
        return m_assemblySeconds;
      }

      /**
       * @brief Clones the problem while retaining its referenced state.
       * @returns A newly allocated copy owned by the caller.
       */
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
