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
  // ---- Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  // ------------------------------------------------------------------------

  template <class TrialFES, class TestFES>
  constexpr
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
     :  m_trialFunction(u),
        m_testFunction(v),
        m_linearForm(v),
        m_bilinearForm(u, v),
        m_assembled(false)
  {}

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator=(const ProblemBody<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>& rhs)
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

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator+=(const LocalBilinearFormIntegratorBase& rhs)
  {
    m_bilinearForm.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator-=(const LocalBilinearFormIntegratorBase& rhs)
  {
    m_bilinearForm.add(UnaryMinus(rhs));
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator+=(const LinearFormIntegratorBase& rhs)
  {
    m_linearForm.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator-=(const LinearFormIntegratorBase& rhs)
  {
    m_linearForm.add(UnaryMinus(rhs));
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator+=(const DirichletBCBase& rhs)
  {
    m_dbcs.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator+=(const PeriodicBCBase& rhs)
  {
    m_pbcs.add(rhs);
    return *this;
  }

  template <class TrialFES, class TestFES>
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>::imposePeriodicBCs()
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
          for (Math::SparseMatrix<Scalar>::InnerIterator it(m_stiffness, parent); it; ++it)
            it.valueRef() = 0;

          // Eliminate the parent row, adding it to the child rows
          IndexMap<Scalar> parentLookup;
          std::vector<IndexMap<Scalar>> childrenLookup(children.size());
          for (size_t col = 0; col < static_cast<size_t>(m_stiffness.cols()); col++)
          {
            bool parentFound = false;
            size_t childrenFound = 0;
            for (Math::SparseMatrix<Scalar>::InnerIterator it(m_stiffness, col); it; ++it)
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
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>::assemble()
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
    const auto& trialFES = trial.getFiniteElementSpace();
    const auto& test = getTestFunction();
    const auto& testFES = test.getFiniteElementSpace();
    for (auto& dbc : m_dbcs)
    {
      dbc.assemble();
      const auto& dofs = dbc.getDOFs();
      if (dbc.isComponent())
      {
        assert(false);
      }
      else
      {
        Math::Kernels::eliminate(m_stiffness, m_mass, dofs);
      }
    }

    // Impose periodic boundary conditions
    imposePeriodicBCs();

    m_assembled = true;

    return *this;
  }

  template <class TrialFES, class TestFES>
  void
  Problem<TrialFES, TestFES, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
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
  // ---- Problem<Tuple<U1, U2, Us...>, Context::Sequential, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  // ------------------------------------------------------------------------

  template <class U1, class U2, class ... Us>
  Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
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
    m_bfa.reset(new BilinearFormTupleSequentialAssembly);
    m_lfa.reset(new LinearFormTupleSequentialAssembly);
    m_us.iapply([&](size_t i, const auto& u)
        { m_trialUUIDMap.right.insert({ i, u.get().getUUID() }); });
    m_vs.iapply([&](size_t i, const auto& v)
        { m_testUUIDMap.right.insert({ i, v.get().getUUID() }); });
  }

  template <class U1, class U2, class ... Us>
  Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>::assemble()
  {
    m_us.apply([](auto& u) { u.get().emplace(); });

    auto bt =
      m_bft.map(
          [](auto& bf)
          {
            auto& u = bf.getTrialFunction();
            auto& v = bf.getTestFunction();
            return Assembly::BilinearFormAssemblyInput(
                u.getFiniteElementSpace(), v.getFiniteElementSpace(),
                bf.getLocalIntegrators(), bf.getGlobalIntegrators());
          });

    auto lt =
      m_lft.map(
          [](auto& lf)
          {
            auto& v = lf.getTestFunction();
            return Assembly::LinearFormAssemblyInput(
                v.getFiniteElementSpace(), lf.getIntegrators());
          });

    // Compute trial offsets
    {
      std::array<size_t, TrialFunctionTuple::Size> sz;
      m_us.map(
            [](const auto& u) { return u.get().getFiniteElementSpace().getSize(); })
          .iapply(
            [&](const Index i, size_t s) { sz[i] = s; });
      m_trialOffsets[0] = 0;
      for (size_t i = 0; i < TrialFunctionTuple::Size - 1; i++)
        m_trialOffsets[i + 1] = sz[i] + m_trialOffsets[i];
    }

    // Compute test offsets
    {
      std::array<size_t, TestFunctionTuple::Size> sz;
      m_vs.map(
            [](const auto& u) { return u.get().getFiniteElementSpace().getSize(); })
          .iapply(
            [&](const Index i, size_t s) { sz[i] = s; });
      m_testOffsets[0] = 0;
      for (size_t i = 0; i < TrialFunctionTuple::Size - 1; i++)
        m_testOffsets[i + 1] = sz[i] + m_testOffsets[i];
    }

    size_t rows =
      m_vs.reduce(
        [](const auto& a, const auto& b)
        {
          return a.get().getFiniteElementSpace().getSize() + b.get().getFiniteElementSpace().getSize();
        });

    size_t cols =
      m_us.reduce(
        [](const auto& a, const auto& b)
        {
          return a.get().getFiniteElementSpace().getSize() + b.get().getFiniteElementSpace().getSize();
        });

    // Compute block offsets to build the triplets
    std::array<Pair<size_t, size_t>, decltype(bt)::Size> boffsets;
    std::array<size_t, decltype(lt)::Size> loffsets;

    m_bft.iapply(
        [&](const Index i, const auto& bf)
        {
          auto ui = m_trialUUIDMap.left.find(bf.getTrialFunction().getUUID());
          auto vi = m_testUUIDMap.left.find(bf.getTestFunction().getUUID());
          if (ui != m_trialUUIDMap.left.end() && vi != m_testUUIDMap.left.end())
            boffsets[i] = Pair{ m_trialOffsets[ui->second], m_testOffsets[vi->second] };
        });

    m_lft.iapply(
        [&](const Index i, const auto& lf)
        {
          auto vi = m_testUUIDMap.left.find(lf.getTestFunction().getUUID());
          if (vi != m_testUUIDMap.left.end())
            loffsets[i] = m_testOffsets[vi->second];
        });

    // Assemble stiffness operator
    m_stiffness = m_bfa->execute(
        Assembly::BilinearFormTupleAssemblyInput(rows, cols, boffsets, bt));

    // Assemble mass vector
    m_mass = m_lfa->execute(
        Assembly::LinearFormTupleAssemblyInput(rows, loffsets, lt));

    // Impose Dirichlet boundary conditions
    m_us.apply(
        [&](const auto& u)
        {
          auto ui = m_trialUUIDMap.left.find(u.get().getUUID());
          size_t offset = m_trialOffsets[ui->second];
          for (auto& dbc : m_dbcs)
          {
            if (dbc.getOperand().getUUID() == u.get().getUUID())
            {
              dbc.assemble();
              const auto& dofs = dbc.getDOFs();
              Math::Kernels::eliminate(m_stiffness, m_mass, dofs, offset);
            }
          }
        });

    m_assembled = true;

    return *this;
  }

  template <class U1, class U2, class ... Us>
  Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>&
  Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::operator=(const ProblemBody<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>& rhs)
  {
    for (auto& bfi : rhs.getLocalBFIs())
    {
      m_bft.apply(
          [&](auto& bf)
          {
            if (bfi.getTrialFunction().getUUID() == bf.getTrialFunction().getUUID() &&
                bfi.getTestFunction().getUUID() == bf.getTestFunction().getUUID())
            {
              bf.add(bfi);
            }
          });
    }

    for (auto& bfi : rhs.getGlobalBFIs())
    {
      m_bft.apply(
          [&](auto& bf)
          {
            if (bfi.getTrialFunction().getUUID() == bf.getTrialFunction().getUUID() &&
                bfi.getTestFunction().getUUID() == bf.getTestFunction().getUUID())
            {
              bf.add(bfi);
            }
          });
    }

    for (auto& lfi : rhs.getLFIs())
    {
      m_lft.apply(
          [&](auto& lf)
          {
            if (lfi.getTestFunction().getUUID() == lf.getTestFunction().getUUID())
            {
              lf.add(UnaryMinus(lfi));
            }
          });
    }

    m_dbcs = rhs.getDBCs();

    return *this;
  }

  template <class U1, class U2, class ... Us>
  void
  Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>
  ::solve(Solver::SolverBase<OperatorType, VectorType>& solver)
  {
     // Assemble the system
     if (!m_assembled)
        assemble();

     // Solve the system AX = B
     solver.solve(m_stiffness, m_guess, m_mass);

     // Recover solutions
     m_us.iapply(
         [&](size_t i, auto& u)
         {
          const size_t n = u.get().getFiniteElementSpace().getSize();
          u.get().getSolution().setWeights(m_guess.segment(m_trialOffsets[i], n));
         });
  }
}

#endif
