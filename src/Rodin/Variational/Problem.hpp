/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_HPP
#define RODIN_VARIATIONAL_PROBLEM_HPP

#include <chrono>

#include "Rodin/Utility.h"

#include "Assembly/Native.h"

#include "GridFunction.h"
#include "DirichletBC.h"

#include "Problem.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
  // ------------------------------------------------------------------------
  // ---- Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>
  // ------------------------------------------------------------------------

  template <class TrialFES, class TestFES>
  constexpr
  Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>
  ::Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
     :  m_trialFunction(u),
        m_testFunction(v),
        m_linearForm(v),
        m_bilinearForm(u, v),
        m_assembled(false)
  {}

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>
  ::operator=(const ProblemBody& rhs)
  {
    for (auto& bfi : rhs.getBFIs())
      getBilinearForm().add(bfi);

    for (auto& lfi : rhs.getLFIs())
      getLinearForm().add(UnaryMinus(lfi)); // Negate every linear form

    m_dbcs = rhs.getDBCs();

     return *this;
  }

  template <class TrialFES, class TestFES>
  void
  Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>::assemble()
  {
    auto& trial = getTrialFunction();
    auto& trialFES = trial.getFiniteElementSpace();

    auto& test = getTestFunction();
    auto& testFES = test.getFiniteElementSpace();

    // Emplace data
    trial.emplace();

    // Assemble both sides
    auto& lhs = getBilinearForm();
    lhs.assemble();
    Math::SparseMatrix& stiffness = lhs.getOperator();

    auto& rhs = getLinearForm();
    rhs.assemble();
    Math::Vector& mass = rhs.getVector();

    // Impose Dirichlet boundary conditions
    if (trialFES == testFES)
    {
      for (auto& dbc : m_dbcs)
      {
         dbc.assemble();
         const auto& dofs = dbc.getDOFs();

         // Move essential degrees of freedom in the LHS to the RHS
         for (const auto& kv : dofs)
         {
            const Index& global = kv.first;
            const auto& dof = kv.second;
            for (Math::SparseMatrix::InnerIterator it(stiffness, global); it; ++it)
               mass.coeffRef(it.row()) -= it.value() * dof;
         }

         // Impose essential degrees of freedom on both sides
         for (const auto& kv : dofs)
         {
           const Index& global = kv.first;
           const auto& dof = kv.second;

            // Set essential DOF in the LHS
            mass.coeffRef(global) = dof;

            // Diagonalize the DOF in the RHS
            stiffness.row(global) *= 0;
            for (Math::SparseMatrix::InnerIterator it(stiffness, global); it; ++it)
            {
              assert(it.row() >= 0);
              it.valueRef() = Scalar(static_cast<size_t>(it.row()) == global);
            }
         }
      }
    }
    else
    {
      assert(false); // Not handled yet
    }

    m_assembled = true;
  }

  template <class TrialFES, class TestFES>
  void
  Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>
  ::solve(const Solver::SolverBase<OperatorType, VectorType>& solver)
  {
     // Assemble the system
     if (!m_assembled)
        assemble();

     // Solve the system AX = B
     solver.solve(getStiffnessOperator(), m_guess, getMassVector());

     // Recover solution
     getTrialFunction().getSolution().setWeights(std::move(m_guess));
  }
}

#endif
