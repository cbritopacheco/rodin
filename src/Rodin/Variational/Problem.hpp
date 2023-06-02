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
   ::operator=(ProblemBody&& rhs)
   {
      Parent::operator=(std::move(rhs));

      for (auto& bfi : getProblemBody().getBFIs())
         getBilinearForm().add(bfi);

      for (auto& lfi : getProblemBody().getLFIs())
         getLinearForm().add(UnaryMinus(lfi)); // Negate every linear form

      return *this;
   }

   template <class TrialFES, class TestFES>
   void
   Problem<TrialFES, TestFES, Context::Serial, Math::SparseMatrix, Math::Vector>::assemble()
   {
      // Assemble both sides
      getLinearForm().assemble();
      m_mass = std::move(getLinearForm().getVector());

      getBilinearForm().assemble();
      m_stiffness = std::move(getBilinearForm().getOperator());

      // Emplace data
      getTrialFunction().emplace();

      // Project values onto the essential boundary and compute essential dofs
      IndexSet trialEssentialDOFs;
      for (const auto& dbc : getProblemBody().getDBCs())
      {
         dbc.project();
         trialEssentialDOFs.insert(dbc.getDOFs());
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
      assert(false);
   }
}

#endif
