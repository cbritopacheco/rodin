/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DENSEPROBLEM_H
#define RODIN_VARIATIONAL_DENSEPROBLEM_H

#include <set>
#include <variant>
#include <functional>

#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Solver/Solver.h"

#include "ForwardDecls.h"

#include "Problem.h"
#include "ProblemBody.h"
#include "LinearForm.h"
#include "BilinearForm.h"
#include "TrialFunction.h"
#include "TestFunction.h"


namespace Rodin::Variational
{
  template <class TrialFES, class TestFES>
  DenseProblem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> DenseProblem<TrialFES, TestFES, Math::Matrix<Scalar>, Math::Vector<Scalar>>;

  /**
   * @defgroup DenseProblemSpecializations DenseProblem Template Specializations
   * @brief Template specializations of the DenseProblem class.
   * @see DenseProblem
   */

  /**
   * @ingroup DenseProblemSpecializations
   * @brief General class to assemble linear systems with `Math::Matrix<Scalar>`
   * and `Math::Vector<Scalar>` types in a serial context.
   */
  template <class TrialFES, class TestFES>
  class DenseProblem<
    TrialFES, TestFES,
    Math::Matrix<typename FormLanguage::Traits<TrialFES>::NumberType>,
    Math::Vector<typename FormLanguage::Traits<TestFES>::NumberType>>
    : public ProblemBase<
        Math::Matrix<typename FormLanguage::Traits<TrialFES>::NumberType>,
        Math::Vector<typename FormLanguage::Traits<TestFES>::NumberType>>
  {
    public:
      using TrialFESNumberType =
        typename FormLanguage::Traits<TrialFES>::NumberType;

      using TestFESNumberType =
        typename FormLanguage::Traits<TestFES>::NumberType;

      using OperatorNumberType =
        decltype(std::declval<TrialFESNumberType>() * std::declval<TestFESNumberType>());

      using VectorNumberType = TestFESNumberType;

      using ContextType = Context::Sequential;

      using OperatorType = Math::Matrix<Scalar>;

      using VectorType = Math::Vector<Scalar>;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<TestFESNumberType>;

      using Parent = ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>>;

      /**
       * @brief Constructs an empty DenseProblem involving the trial function @f$ u @f$
       * and the test function @f$ v @f$.
       *
       * @param[in,out] u Trial function
       * @param[in,out] v %Test function
       */
      constexpr
      DenseProblem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
        : m_trialFunction(u),
          m_testFunction(v),
          m_linearForm(v),
          m_bilinearForm(u, v),
          m_assembled(false)
      {}

      /**
       * @brief Deleted copy constructor.
       */
      DenseProblem(const DenseProblem& other) = delete;

      /**
       * @brief Deleted copy assignment operator.
       */
      void operator=(const DenseProblem& other) = delete;

      constexpr
      TrialFunction<TrialFES>& getTrialFunction()
      {
        return m_trialFunction;
      }

      constexpr
      TestFunction<TestFES>& getTestFunction()
      {
        return m_testFunction;
      }

      constexpr
      const TrialFunction<TrialFES>& getTrialFunction() const
      {
        return m_trialFunction.get();
      }

      constexpr
      const TestFunction<TestFES>& getTestFunction() const
      {
        return m_testFunction.get();
      }

      constexpr
      const PeriodicBoundary& getPeriodicBoundary() const
      {
        return m_pbcs;
      }

      constexpr
      const EssentialBoundary& getEssentialBoundary() const
      {
        return m_dbcs;
      }

      DenseProblem& imposePeriodicBCs()
      {
        const auto& trial = getTrialFunction();
        const auto& trialFES = trial.getFiniteElementSpace();

        const auto& test = getTestFunction();
        const auto& testFES = test.getFiniteElementSpace();

        if (trialFES == testFES)
        {
          for (auto& pbc : m_pbcs)
          {
            pbc.assemble();
            const auto& dofs = pbc.getDOFs();

            std::deque<Index> q;
            IndexSet dependents;
            dependents.reserve(dofs.size());
            for (const auto& [k, v] : dofs)
              dependents.insert(v.first.begin(), v.first.end());

            for (auto it = dofs.begin(); it != dofs.end(); ++it)
            {
              const Index k = it->first;
              if (!dependents.contains(k))
                q.push_front(k);
            }

            // Perform breadth-first traversal
            while (q.size() > 0)
            {
              const Index parent = q.back();
              assert(m_stiffness.rows() >= 0);
              assert(m_stiffness.cols() >= 0);
              assert(parent < static_cast<size_t>(m_stiffness.rows()));
              assert(parent < static_cast<size_t>(m_stiffness.cols()));
              q.pop_back();

              auto find = dofs.find(parent);
              if (find == dofs.end())
                continue;

              const auto& [children, coeffs] = find->second;
              assert(children.size() > 0);

              for (const auto& child : children)
                q.push_front(child);

              assert(children.size() == coeffs.size());
              const size_t count = children.size();

              // Eliminate the parent column, adding it to the child columns
              for (size_t i = 0; i < count; i++)
              {
                const Scalar coeff = coeffs.coeff(i);
                const Index child = children.coeff(i);
                m_stiffness.col(child) += coeff * m_stiffness.col(parent);
              }

              m_stiffness.col(parent).setZero();

              // Eliminate the parent row, adding it to the child rows
              IndexMap<Scalar> parentLookup;
              std::vector<IndexMap<Scalar>> childrenLookup(children.size());
              for (size_t col = 0; col < static_cast<size_t>(m_stiffness.cols()); col++)
              {
                bool parentFound = false;
                size_t childrenFound = 0;
                for (Math::Matrix<Scalar>::InnerIterator it(m_stiffness, col); it; ++it)
                {
                  if (parentFound && childrenFound == count)
                  {
                    break;
                  }
                  else
                  {
                    const Index row = it.row();
                    if (row == parent)
                    {
                      parentLookup[col] = it.value();
                      m_stiffness.coeffRef(it.row(), it.col()) = 0;
                      parentFound = true;
                    }
                    else
                    {
                      for (size_t i = 0; i < count; i++)
                      {
                        const Index child = children.coeff(i);
                        if (row == child)
                        {
                          childrenLookup[i][col] = it.value();
                          childrenFound += 1;
                        }
                      }
                    }
                  }
                }
              }

              for (const auto& [col, value] : parentLookup)
              {
                for (size_t i = 0; i < count; i++)
                {
                  const Scalar coeff = coeffs.coeff(i);
                  childrenLookup[i][col] += coeff * value;
                }
              }

              for (size_t i = 0; i < count; i++)
              {
                const Index child = children.coeff(i);
                for (const auto& [col, value] : childrenLookup[i])
                  m_stiffness.coeffRef(child, col) = value;
              }

              // Eliminate the parent entry, adding it to the child entries
              for (size_t i = 0; i < count; i++)
              {
                const Scalar coeff = coeffs.coeff(i);
                const Index child = children.coeff(i);
                m_mass.coeffRef(child) += coeff * m_mass.coeff(parent);
              }
              m_mass.coeffRef(parent) = 0;
            }

            for (const auto& [parent, node] : dofs)
            {
              m_stiffness.coeffRef(parent, parent) = 1.0;
              const auto& [children, coeffs] = node;
              assert(children.size() >= 0);
              for (size_t i = 0; i < static_cast<size_t>(children.size()); i++)
              {
                const Scalar coeff = coeffs.coeff(i);
                const Index child = children.coeff(i);
                m_stiffness.coeffRef(parent, child) = -coeff;
              }
            }
          }
        }
        else
        {
          assert(false); // Not handled yet
        }

        return *this;
      }

      DenseProblem& imposeDirichletBCs()
      {
        const auto& trial = getTrialFunction();
        const auto& trialFES = trial.getFiniteElementSpace();
        const auto& test = getTestFunction();
        const auto& testFES = test.getFiniteElementSpace();
        for (auto& dbc : m_dbcs)
        {
          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          // Move essential degrees of freedom in the LHS to the RHS
          for (const auto& kv : dofs)
          {
            const Index& global = kv.first;
            const auto& dof = kv.second;
            m_mass -= dof * m_stiffness.col(global);
          }
          for (const auto& [global, dof] : dofs)
          {
            // Impose essential degrees of freedom on RHS
            m_mass.coeffRef(global) = dof;

            // Impose essential degrees of freedom on LHS
            m_stiffness.col(global).setZero();
            m_stiffness.row(global).setZero();
            m_stiffness.coeffRef(global, global) = 1;
          }
        }
        return *this;
      }

      DenseProblem& assemble() override
      {
        auto& trial = getTrialFunction();

        // Emplace data
        trial.emplace();

        // Assemble both sides
        m_linearForm.assemble();
        m_mass = std::move(m_linearForm.getVector());

        m_bilinearForm.assemble();
        m_stiffness = std::move(m_bilinearForm.getOperator());

        for (auto& bf : m_bfs)
        {
          bf.assemble();
          m_stiffness += bf.getOperator();
        }

        // Impose Dirichlet boundary conditions
        imposeDirichletBCs();

        // Impose periodic boundary conditions
        imposePeriodicBCs();

        m_assembled = true;

        return *this;
      }

      void solve(Solver::SolverBase<OperatorType, VectorType>& solver) override
      {
         // Assemble the system
         if (!m_assembled)
            assemble();

         // Solve the system AX = B
         solver.solve(m_stiffness, m_guess, m_mass);

         // Recover solution
         getTrialFunction().getSolution().setWeights(std::move(m_guess));
      }

      DenseProblem& operator=(const ProblemBody<OperatorType, VectorType>& rhs) override
      {
        for (auto& bfi : rhs.getLocalBFIs())
          m_bilinearForm.add(bfi);
        for (auto& bfi : rhs.getGlobalBFIs())
          m_bilinearForm.add(bfi);

        for (auto& lfi : rhs.getLFIs())
          m_linearForm.add(UnaryMinus(lfi)); // Negate every linear form

        m_bfs = rhs.getBFs();

        m_dbcs = rhs.getDBCs();
        m_pbcs = rhs.getPBCs();

        return *this;
      }

      virtual VectorType& getMassVector() override
      {
        return m_mass;
      }

      virtual const VectorType& getMassVector() const override
      {
        return m_mass;
      }

      virtual OperatorType& getStiffnessOperator() override
      {
        return m_stiffness;
      }

      virtual const OperatorType& getStiffnessOperator() const override
      {
        return m_stiffness;
      }

      virtual DenseProblem* copy() const noexcept override
      {
        assert(false);
        return nullptr;
      }

    private:
      std::reference_wrapper<TrialFunction<TrialFES>> m_trialFunction;
      std::reference_wrapper<TestFunction<TestFES>>   m_testFunction;

      LinearForm<TestFES, VectorType> m_linearForm;
      BilinearForm<TrialFES, TestFES, OperatorType> m_bilinearForm;

      FormLanguage::List<BilinearFormBase<OperatorType>> m_bfs;

      EssentialBoundary m_dbcs;
      PeriodicBoundary  m_pbcs;

      bool m_assembled;
      VectorType      m_mass;
      VectorType      m_guess;
      OperatorType    m_stiffness;
  };
}

#include "DenseProblem.hpp"

#endif

