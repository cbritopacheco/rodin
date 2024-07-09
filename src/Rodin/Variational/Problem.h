/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_H
#define RODIN_VARIATIONAL_PROBLEM_H

#include <set>
#include <variant>
#include <functional>
#include <boost/mp11.hpp>

#include "Rodin/Pair.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/Solver/Solver.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/Math/BlockSparseMatrix.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Tuple.h"
#include "Rodin/Utility/Extract.h"
#include "Rodin/Utility/Product.h"
#include "Rodin/Utility/Wrap.h"

#include "ForwardDecls.h"

#include "ProblemBody.h"
#include "LinearForm.h"
#include "BilinearForm.h"
#include "TrialFunction.h"
#include "TestFunction.h"


namespace Rodin::Variational
{
  /**
   * @defgroup ProblemSpecializations Problem Template Specializations
   * @brief Template specializations of the Problem class.
   * @see Problem
   */

  /**
   * @brief Abstract base class for variational problems.
   */
  template <class Operator, class Vector, class Scalar>
  class ProblemBase : public FormLanguage::Base
  {
    public:
      using VectorType = Vector;

      using OperatorType = Operator;

      using VectorScalarType = typename FormLanguage::Traits<Vector>::ScalarType;

      using OperatorScalarType = typename FormLanguage::Traits<Operator>::ScalarType;

      using ScalarType = Scalar;

      ProblemBase() = default;

      ProblemBase(ProblemBase&& other) = default;

      ProblemBase(const ProblemBase& other) = default;

      virtual ProblemBase& operator=(const ProblemBody<OperatorType, VectorType, ScalarType>& rhs) = 0;

      virtual void solve(Solver::SolverBase<OperatorType, VectorType>& solver) = 0;

      /**
       * @brief Assembles the underlying linear system to solve.
       */
      virtual ProblemBase& assemble() = 0;

      /**
       * @returns Reference to the stiffness operator.
       *
       * This must be called only after assemble() has been called.
       */
      virtual OperatorType& getStiffnessOperator() = 0;

      /**
       * @returns Constant reference to the stiffness operator.
       *
       * This must be called only after assemble() has been called.
       */
      virtual const OperatorType& getStiffnessOperator() const = 0;

      /**
       * @returns Reference to the mass vector.
       *
       * This must be called only after assemble() has been called.
       */
      virtual VectorType& getMassVector() = 0;

      /**
       * @returns Constant reference to the mass vector.
       *
       * This must be called only after assemble() has been called.
       */
      virtual const VectorType& getMassVector() const = 0;

      virtual ProblemBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup ProblemSpecializations
   * @brief General class to assemble linear systems with `Math::SparseMatrix`
   * and `Math::Vector` types in a sequential context.
   */
  template <class TrialFES, class TestFES>
  class Problem<
    TrialFES, TestFES,
    Math::SparseMatrix<
      typename FormLanguage::Mult<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TrialFES>::ScalarType>
      ::Type>,
    Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    : public ProblemBase<
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TrialFES>::ScalarType>
          ::Type>,
        Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>,
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TrialFES>::ScalarType>
          ::Type>
  {
    public:
      using TrialFESScalarType =
        typename FormLanguage::Traits<TrialFES>::ScalarType;

      using TestFESScalarType =
        typename FormLanguage::Traits<TestFES>::ScalarType;

      using OperatorScalarType =
        typename FormLanguage::Mult<TrialFESScalarType, TestFESScalarType>::Type;

      using VectorScalarType = TestFESScalarType;

      using ScalarType = OperatorScalarType;

      using ContextType = Context::Sequential;

      using OperatorType = Math::SparseMatrix<OperatorScalarType>;

      using VectorType = Math::Vector<TestFESScalarType>;

      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<TestFESScalarType>;

      using Parent = ProblemBase<OperatorType, VectorType, ScalarType>;

      /**
       * @brief Constructs an empty problem involving the trial function @f$ u @f$
       * and the test function @f$ v @f$.
       *
       * @param[in,out] u Trial function
       * @param[in,out] v %Test function
       */
      constexpr
      Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v)
         :  m_trialFunction(u),
            m_testFunction(v),
            m_linearForm(v),
            m_bilinearForm(u, v),
            m_assembled(false)
      {}

      /**
       * @brief Deleted copy constructor.
       */
      Problem(const Problem& other) = delete;

      /**
       * @brief Deleted copy assignment operator.
       */
      void operator=(const Problem& other) = delete;

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

      Problem& imposePeriodicBCs()
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
                const OperatorScalarType coeff = coeffs.coeff(i);
                const Index child = children.coeff(i);
                m_stiffness.col(child) += coeff * m_stiffness.col(parent);
              }

              // Assumes CCS format
              for (typename OperatorType::InnerIterator it(m_stiffness, parent); it; ++it)
                it.valueRef() = 0;

              // Eliminate the parent row, adding it to the child rows
              IndexMap<OperatorScalarType> parentLookup;
              std::vector<IndexMap<OperatorScalarType>> childrenLookup(children.size());
              for (size_t col = 0; col < static_cast<size_t>(m_stiffness.cols()); col++)
              {
                bool parentFound = false;
                size_t childrenFound = 0;
                for (typename OperatorType::InnerIterator it(m_stiffness, col); it; ++it)
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
                  const OperatorScalarType coeff = coeffs.coeff(i);
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
                const OperatorScalarType coeff = coeffs.coeff(i);
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
                const OperatorScalarType coeff = coeffs.coeff(i);
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

      const LinearForm<TestFES, VectorType>& getLinearForm() const
      {
        return m_linearForm;
      }

      Problem& assemble() override
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

      Problem& operator=(const ProblemBody<OperatorType, VectorType, ScalarType>& rhs) override
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

      virtual Problem* copy() const noexcept override
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

      bool            m_assembled;
      VectorType      m_mass;
      VectorType      m_guess;
      OperatorType    m_stiffness;
  };

  template <class TrialFES, class TestFES>
  Problem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> Problem<TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
        Math::Vector<
          typename FormLanguage::Traits<TestFES>::ScalarType>>;

  template <class U1, class U2, class ... Us>
  class Problem<
      Tuple<U1, U2, Us...>, Math::SparseMatrix<Real>, Math::Vector<Real>>
    : public ProblemBase<Math::SparseMatrix<Real>, Math::Vector<Real>, Real>
  {

    template <class T>
    struct IsTrialOrTestFunction
    {
      static constexpr bool Value = IsTrialFunction<T>::Value || IsTestFunction<T>::Value;
    };

    static_assert(Utility::ParameterPack<U1, U2, Us...>::template All<IsTrialOrTestFunction>::Value);

    public:
      using ScalarType = Real;

      using ContextType = Context::Sequential;

      using OperatorType = Math::SparseMatrix<ScalarType>;

      using VectorType = Math::Vector<ScalarType>;

      using Parent = ProblemBase<OperatorType, VectorType, Real>;

    private:
      template <class T>
      struct GetFES;

      template <class T>
      struct GetFES<std::reference_wrapper<T>>
      {
        using Type = typename FormLanguage::Traits<T>::FES;
      };

      template <class T>
      struct IsTrialFunctionReferenceWrapper
      {
        static constexpr Boolean Value = false;
      };

      template <class T>
      struct IsTrialFunctionReferenceWrapper<std::reference_wrapper<T>>
      {
        static constexpr Boolean Value = IsTrialFunction<T>::Value;
      };

      template <class T>
      struct IsTestFunctionReferenceWrapper
      {
        static constexpr Boolean Value = false;
      };

      template <class T>
      struct IsTestFunctionReferenceWrapper<std::reference_wrapper<T>>
      {
        static constexpr Boolean Value = IsTestFunction<T>::Value;
      };

      using TrialFunctionTuple =
        decltype(std::declval<
          Tuple<
            std::reference_wrapper<U1>,
            std::reference_wrapper<U2>,
            std::reference_wrapper<Us>...>>()
            .template filter<IsTrialFunctionReferenceWrapper>());

      using TestFunctionTuple =
        decltype(std::declval<
          Tuple<
            std::reference_wrapper<U1>,
            std::reference_wrapper<U2>,
            std::reference_wrapper<Us>...>>()
            .template filter<IsTestFunctionReferenceWrapper>());

      using TrialFESTuple = typename Utility::Extract<TrialFunctionTuple>::template Type<GetFES>;

      using TestFESTuple = typename Utility::Extract<TestFunctionTuple>::template Type<GetFES>;

      template <class TrialFES, class TestFES>
      using BilinearFormType = BilinearForm<TrialFES, TestFES, OperatorType>;

      template <class TestFES>
      using LinearFormType = LinearForm<TestFES, VectorType>;

      using BilinearFormTuple =
        typename Utility::Product<TrialFESTuple, TestFESTuple>::template Type<BilinearFormType>;

      using LinearFormTuple =
        typename Utility::Wrap<TestFESTuple>::template Type<LinearFormType>;

    public:
      using BilinearFormTupleSequentialAssembly =
        Assembly::Sequential<OperatorType, BilinearFormTuple>;

      using LinearFormTupleSequentialAssembly =
        Assembly::Sequential<VectorType, LinearFormTuple>;

      Problem(U1& u1, U2& u2, Us&... us)
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

      Problem& assemble() override
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

      void solve(Solver::SolverBase<OperatorType, VectorType>& solver) override
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

      Problem& operator=(const ProblemBody<OperatorType, VectorType, Real>& rhs) override
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

      inline
      const auto& getTrialOffsets() const
      {
        return m_trialOffsets;
      }

      inline
      const auto& getTestOffsets() const
      {
        return m_testOffsets;
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

      virtual Problem* copy() const noexcept override
      {
        assert(false);
        return nullptr;
      }

    private:
      TrialFunctionTuple m_us;
      TestFunctionTuple  m_vs;

      LinearFormTuple   m_lft;
      BilinearFormTuple m_bft;

      EssentialBoundary m_dbcs;

      bool            m_assembled;
      VectorType      m_mass;
      VectorType      m_guess;
      OperatorType    m_stiffness;

      std::array<size_t, TrialFunctionTuple::Size> m_trialOffsets;
      std::array<size_t, TestFunctionTuple::Size> m_testOffsets;

      boost::bimap<FormLanguage::Base::UUID, size_t> m_trialUUIDMap;
      boost::bimap<FormLanguage::Base::UUID, size_t> m_testUUIDMap;

      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearFormTuple>> m_bfa;
      std::unique_ptr<Assembly::AssemblyBase<VectorType, LinearFormTuple>> m_lfa;
  };

  template <class U1, class U2, class ... Us>
  Problem(U1& u1, U2& u2, Us&... us)
    -> Problem<Tuple<U1, U2, Us...>, Math::SparseMatrix<Real>, Math::Vector<Real>>;
}

#include "Problem.hpp"

#endif
