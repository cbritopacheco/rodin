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

#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Assembly/Input.h"

#include "Exceptions/TrialFunctionMismatchException.h"
#include "Exceptions/TestFunctionMismatchException.h"

#include "GridFunction.h"
#include "DirichletBC.h"

#include "Problem.h"
#include "UnaryMinus.h"

namespace Rodin::Variational
{
  // ------------------------------------------------------------------------
  // ---- Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  // ------------------------------------------------------------------------

  template <class TrialFES, class TestFES>
  constexpr
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
     :  m_trialFunction(u),
        m_testFunction(v),
        m_linearForm(v),
        m_bilinearForm(u, v),
        m_assembled(false)
  {}

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator=(const ProblemBody<Math::SparseMatrix, Math::Vector>& rhs)
  {
    for (auto& bfi : rhs.getLocalBFIs())
    {
      if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
        TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
      if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
        TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
      m_bilinearForm.add(bfi);
    }

    for (auto& bfi : rhs.getGlobalBFIs())
    {
      if (bfi.getTrialFunction().getUUID() != getTrialFunction().getUUID())
        TrialFunctionMismatchException(bfi.getTrialFunction()) << Alert::Raise;
      if (bfi.getTestFunction().getUUID() != getTestFunction().getUUID())
        TestFunctionMismatchException(bfi.getTestFunction()) << Alert::Raise;
      m_bilinearForm.add(bfi);
    }

    for (auto& lfi : rhs.getLFIs())
    {
      if (lfi.getTestFunction().getUUID() != getTestFunction().getUUID())
      {

        std::cout << getTestFunction().getUUID() << std::endl;
        TestFunctionMismatchException(lfi.getTestFunction()) << Alert::Raise;
      }
      m_linearForm.add(UnaryMinus(lfi)); // Negate every linear form
    }

    m_bfs = rhs.getBFs();

    m_dbcs = rhs.getDBCs();
    m_pbcs = rhs.getPBCs();

    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator+=(const LocalBilinearFormIntegratorBase& rhs)
  {
    if (rhs.getTrialFunction().getUUID() != getTrialFunction().getUUID())
      TrialFunctionMismatchException(rhs.getTrialFunction()) << Alert::Raise;
    if (rhs.getTestFunction().getUUID() != getTestFunction().getUUID())
      TestFunctionMismatchException(rhs.getTestFunction()) << Alert::Raise;
    m_bilinearForm.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator-=(const LocalBilinearFormIntegratorBase& rhs)
  {
    if (rhs.getTrialFunction().getUUID() != getTrialFunction().getUUID())
      TrialFunctionMismatchException(rhs.getTrialFunction()) << Alert::Raise;
    if (rhs.getTestFunction().getUUID() != getTestFunction().getUUID())
      TestFunctionMismatchException(rhs.getTestFunction()) << Alert::Raise;
    m_bilinearForm.add(UnaryMinus(rhs));
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator+=(const LinearFormIntegratorBase& rhs)
  {
    if (rhs.getTestFunction().getUUID() != getTestFunction().getUUID())
      TestFunctionMismatchException(rhs.getTestFunction()) << Alert::Raise;
    m_linearForm.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator-=(const LinearFormIntegratorBase& rhs)
  {
    if (rhs.getTestFunction().getUUID() != getTestFunction().getUUID())
      TestFunctionMismatchException(rhs.getTestFunction()) << Alert::Raise;
    m_linearForm.add(UnaryMinus(rhs));
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator+=(const DirichletBCBase& rhs)
  {
    m_dbcs.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::operator+=(const PeriodicBCBase& rhs)
  {
    m_pbcs.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>::imposeDirichletBCs()
  {
    const auto& trial = getTrialFunction();
    const auto& trialFES = trial.getFiniteElementSpace();
    const auto& test = getTestFunction();
    const auto& testFES = test.getFiniteElementSpace();
    Scalar* const valuePtr = m_stiffness.valuePtr();
    Math::SparseMatrix::StorageIndex* const outerPtr = m_stiffness.outerIndexPtr();
    Math::SparseMatrix::StorageIndex* const innerPtr = m_stiffness.innerIndexPtr();
    for (auto& dbc : m_dbcs)
    {
      dbc.assemble();
      const auto& dofs = dbc.getDOFs();
      // Move essential degrees of freedom in the LHS to the RHS
      for (const auto& kv : dofs)
      {
        const Index& global = kv.first;
        const auto& dof = kv.second;
        for (Math::SparseMatrix::InnerIterator it(m_stiffness, global); it; ++it)
           m_mass.coeffRef(it.row()) -= it.value() * dof;
      }
      for (const auto& [global, dof] : dofs)
      {
        // Impose essential degrees of freedom on RHS
        m_mass.coeffRef(global) = dof;

        // Impose essential degrees of freedom on LHS
        for (auto i = outerPtr[global]; i < outerPtr[global + 1]; ++i)
        {
          assert(innerPtr[i] >= 0);
          // Assumes CCS format
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
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>::imposePeriodicBCs()
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

          // Assumes CCS format
          for (Math::SparseMatrix::InnerIterator it(m_stiffness, parent); it; ++it)
            it.valueRef() = 0;

          // Eliminate the parent row, adding it to the child rows
          IndexMap<Scalar> parentLookup;
          std::vector<IndexMap<Scalar>> childrenLookup(children.size());
          for (size_t col = 0; col < static_cast<size_t>(m_stiffness.cols()); col++)
          {
            bool parentFound = false;
            size_t childrenFound = 0;
            for (Math::SparseMatrix::InnerIterator it(m_stiffness, col); it; ++it)
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
                  it.valueRef() = 0;
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

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>::assemble()
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

  template <class TrialFES, class TestFES>
  void
  Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::solve(Solver::SolverBase<OperatorType, VectorType>& solver)
  {
     // Assemble the system
     if (!m_assembled)
        assemble();

     // Solve the system AX = B
     solver.solve(m_stiffness, m_guess, m_mass);

     // Recover solution
     getTrialFunction().getSolution().setWeights(std::move(m_guess));
  }

  // ------------------------------------------------------------------------
  // ---- Problem<Tuple<U1, U2, Us...>, Context::Sequential, Math::SparseMatrix, Math::Vector>
  // ------------------------------------------------------------------------

  template <class U1, class U2, class ... Us>
  Problem<Tuple<U1, U2, Us...>, Context::Sequential, Math::SparseMatrix, Math::Vector>
  ::Problem(U1& u1, U2& u2, Us&... us)
    : m_us(
        Tuple{std::ref(u1), std::ref(u2), std::ref(us)...}
        .template filter<IsTrialFunctionReferenceWrapper>()),
      m_vs(
        Tuple{std::ref(u1), std::ref(u2), std::ref(us)...}
        .template filter<IsTestFunctionReferenceWrapper>()),
      m_lft(m_vs.map(
            [](const auto& v)
            { return LinearFormType<
                typename std::decay_t<
                typename Utility::UnwrapRefDecay<decltype(v)>::Type>::FES>(v.get());
            })),
      m_bft(m_us.product(
            [](const auto& u, const auto& v) { return Pair(u, v); }, m_vs)
                .map(
                  [](const auto& uv)
                  { return BilinearFormType<
                      typename std::decay_t<
                      typename Utility::UnwrapRefDecay<decltype(uv.first())>::Type>::FES,
                      typename std::decay_t<
                      typename Utility::UnwrapRefDecay<decltype(uv.second())>::Type>::FES>(
                          uv.first().get(), uv.second().get());
                  }))
  {
    m_assembly.reset(new SequentialAssembly);
  }

  template <class U1, class U2, class ... Us>
  Problem<Tuple<U1, U2, Us...>, Context::Sequential, Math::SparseMatrix, Math::Vector>&
  Problem<Tuple<U1, U2, Us...>, Context::Sequential, Math::SparseMatrix, Math::Vector>::assemble()
  {
    m_us.apply([](auto& u) { u.get().emplace(); });
    auto t =
      m_bft.map(
          [](auto& bf)
          {
            auto& u = bf.getTrialFunction();
            auto& v = bf.getTestFunction();
            return Assembly::BilinearFormAssemblyInput(
                u.getFiniteElementSpace(), v.getFiniteElementSpace(),
                bf.getLocalIntegrators(), bf.getGlobalIntegrators());
          });

    m_assembled = true;

    return *this;
  }
}

#endif
