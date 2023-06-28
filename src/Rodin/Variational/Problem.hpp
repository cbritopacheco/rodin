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
      m_bilinearForm.add(bfi);

    for (auto& lfi : rhs.getLFIs())
      m_linearForm.add(UnaryMinus(lfi)); // Negate every linear form

    m_dbcs = rhs.getDBCs();

     return *this;
  }

  namespace Internal
  {
    class IteratorTripletWrapper
    {
      public:
        IteratorTripletWrapper(UnorderedMap<std::pair<Index, Index>, Scalar>::iterator it) : m_it(it) {}
        inline bool operator==(const IteratorTripletWrapper& other) const { return m_it == other.m_it; }
        inline bool operator!=(const IteratorTripletWrapper& other) const { return m_it != other.m_it; }
        inline IteratorTripletWrapper& operator++() { ++m_it; return *this; }
        inline const Eigen::Triplet<Scalar>& operator*() const
        { m_value = Eigen::Triplet<Scalar>(m_it->first.first, m_it->first.second, m_it->second); return m_value; }
        inline const Eigen::Triplet<Scalar>* operator->() const { m_value = Eigen::Triplet<Scalar>(m_it->first.first, m_it->first.second, m_it->second); return &m_value; }
      private:
        UnorderedMap<std::pair<Index, Index>, Scalar>::iterator m_it;
        mutable Eigen::Triplet<Scalar> m_value;
    };
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
    m_linearForm.assemble();
    m_mass = std::move(m_linearForm.getVector());

    m_bilinearForm.assemble();
    m_stiffness = std::move(m_bilinearForm.getOperator());

    // Impose Dirichlet boundary conditions
    if (trialFES == testFES)
    {
      for (auto& dbc : m_dbcs)
      {
        dbc.assemble();
        const auto& dofs = dbc.getDOFs();

        // Impose essential degrees of freedom on RHS
        for (const auto& [global, dof] : dofs)
        {
          m_mass.coeffRef(global) = dof;

          Math::SparseMatrix::StorageIndex* outerPtr = m_stiffness.outerIndexPtr();
          Math::SparseMatrix::StorageIndex* innerPtr = m_stiffness.innerIndexPtr();
          Scalar* valuePtr = m_stiffness.valuePtr();
          for (auto i = outerPtr[global]; i < outerPtr[global + 1]; ++i)
          {
            assert(innerPtr[i] >= 0);
            const Index row = innerPtr[i];
            valuePtr[i] = Scalar(row == global);
            if (row != global)
            {
              for (auto k = outerPtr[row]; 1; k++)
              {
                if (static_cast<Index>(innerPtr[k]) == global)
                {
                   valuePtr[k] = 0.0;
                   break;
                }
              }
            }
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
     solver.solve(m_stiffness, m_guess, m_mass);

     // Recover solution
     getTrialFunction().getSolution().setWeights(std::move(m_guess));
  }
}

#endif
