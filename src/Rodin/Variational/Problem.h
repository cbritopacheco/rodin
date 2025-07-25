/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_H
#define RODIN_VARIATIONAL_PROBLEM_H

#include <variant>
#include <functional>
#include <boost/mp11.hpp>

#include "Rodin/Pair.h"
#include "Rodin/Tuple.h"
#include "Rodin/Alert.h"
#include "Rodin/Geometry.h"
#include "Rodin/Solver/Solver.h"
#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/SparseMatrix.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Utility/Extract.h"
#include "Rodin/Utility/Product.h"
#include "Rodin/Utility/Wrap.h"
#include "Rodin/Assembly/Input.h"

#include "ForwardDecls.h"

#include "ProblemBody.h"
#include "LinearForm.h"
#include "BilinearForm.h"

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
      using OperatorType = Operator;

      using VectorType = Vector;

      using ScalarType = Scalar;

      ProblemBase() = default;

      ProblemBase(ProblemBase&& other) = default;

      ProblemBase(const ProblemBase& other) = default;

      virtual ProblemBase& operator=(
          const ProblemBody<OperatorType, VectorType, ScalarType>& rhs) = 0;

      virtual void solve(
          Solver::SolverBase<OperatorType, VectorType, ScalarType>& solver) = 0;

      virtual Math::LinearSystem<OperatorType, VectorType>& getLinearSystem() = 0;

      virtual const Math::LinearSystem<OperatorType, VectorType>& getLinearSystem() const = 0;

      /**
       * @brief Assembles the underlying linear system to solve.
       */
      virtual ProblemBase& assemble() = 0;

      virtual ProblemBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup ProblemSpecializations
   * @brief General class to assemble linear systems with `Operator`
   * and `Vector` generic types in a sequential context.
   */
  template <class LinearSystem, class U, class V>
  class ProblemUV
    : public ProblemBase<
        typename FormLanguage::Traits<LinearSystem>::OperatorType,
        typename FormLanguage::Traits<LinearSystem>::VectorType,
        typename FormLanguage::Traits<LinearSystem>::ScalarType>
  {
    public:
      using TrialFunctionType =
        U;

      using TestFunctionType =
        V;

      using LinearSystemType =
        LinearSystem;

      using SolutionType =
        typename FormLanguage::Traits<TrialFunctionType>::SolutionType;

      using TrialFESType =
        typename FormLanguage::Traits<U>::FESType;

      using TestFESType =
        typename FormLanguage::Traits<V>::FESType;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      using ProblemBodyType = ProblemBody<OperatorType, VectorType, ScalarType>;

      using TrialFESScalarType =
        typename FormLanguage::Traits<TrialFESType>::ScalarType;

      using TestFESScalarType =
        typename FormLanguage::Traits<TestFESType>::ScalarType;

      using LinearFormIntegratorBaseType =
        LinearFormIntegratorBase<TestFESScalarType>;

      using Parent = ProblemBase<OperatorType, VectorType, ScalarType>;

      constexpr
      ProblemUV(U& u, V& v)
        : m_assembled(false),
          m_trialFunction(u), m_testFunction(v)
      {}

      ProblemUV(const ProblemUV& other)
        : Parent(other),
          m_assembled(other.m_assembled),
          m_trialFunction(other.m_trialFunction.get()),
          m_testFunction(other.m_testFunction.get()),
          m_pb(other.m_pb)
      {}

      ProblemUV(ProblemUV&& other)
        : Parent(std::move(other)),
          m_assembled(std::exchange(other.m_assembled, false)),
          m_trialFunction(std::move(other.m_trialFunction)),
          m_testFunction(std::move(other.m_testFunction)),
          m_pb(std::move(other.m_pb))
      {}

      ProblemUV& operator=(const ProblemUV& other)
      {
        if (this != &other)
        {
          m_assembled = other.m_assembled;
          m_trialFunction = other.m_trialFunction.get();
          m_testFunction = other.m_testFunction.get();
          m_pb = other.m_pb;
        }
        return *this;
      }

      ProblemUV& operator=(ProblemUV&& other) noexcept
      {
        if (this != &other)
        {
          m_assembled = std::exchange(other.m_assembled, false);
          m_trialFunction = std::move(other.m_trialFunction);
          m_testFunction = std::move(other.m_testFunction);
          m_pb = std::move(other.m_pb);
        }
        return *this;
      }

      constexpr
      TrialFunctionType& getTrialFunction()
      {
        return m_trialFunction;
      }

      constexpr
      TestFunctionType& getTestFunction()
      {
        return m_testFunction;
      }

      constexpr
      const TrialFunctionType& getTrialFunction() const
      {
        return m_trialFunction.get();
      }

      constexpr
      const TestFunctionType& getTestFunction() const
      {
        return m_testFunction.get();
      }

      ProblemUV& assemble() override
      {
        using LinearFormType = LinearForm<TestFESType, VectorType>;
        using BilinearFormType =
          BilinearForm<SolutionType, TrialFESType, TestFESType, OperatorType>;

        auto& pb = m_pb;
        auto& u = this->getTrialFunction();
        auto& v = this->getTestFunction();
        auto& axb = this->getLinearSystem();
        auto& mass = axb.getVector();
        auto& stiffness = axb.getOperator();
        auto& bfs = pb.getBFs();
        auto& dbcs = pb.getDBCs();
        auto& pbcs = pb.getPBCs();

        LinearFormType lf(v);
        for (auto& lfi : pb.getLFIs())
          lf.add(UnaryMinus(lfi)); // Negate every linear form integrator
        lf.assemble();
        mass = std::move(lf.getVector());

        BilinearFormType bf(u, v);
        for (auto& bfi : pb.getLocalBFIs())
          bf.add(bfi);
        for (auto& bfi : pb.getGlobalBFIs())
          bf.add(bfi);
        for (auto& _bf : bfs)
        {
          _bf.assemble();
          stiffness += _bf.getOperator();
        }
        bf.assemble();
        stiffness = std::move(bf.getOperator());

        // Impose Dirichlet boundary conditions
        auto& trial = getTrialFunction();
        const auto& trialFES = trial.getFiniteElementSpace();
        const auto& test = getTestFunction();
        const auto& testFES = test.getFiniteElementSpace();
        for (auto& dbc : dbcs)
        {
          dbc.assemble();
          const auto& dofs = dbc.getDOFs();
          if (dbc.isComponent())
          {
            assert(false);
          }
          else
          {
            axb.eliminate(dofs);
          }
        }

        // Impose periodic boundary conditions
        if (trialFES == testFES)
        {
          for (auto& pbc : pbcs)
          {
            pbc.assemble();
            const auto& dofs = pbc.getDOFs();

            if (pbc.isComponent())
            {
              assert(false);
            }
            else
            {
              axb.merge(dofs);
            }
          }
        }
        else
        {
          assert(false); // Not handled yet
        }

        m_assembled = true;

        return *this;
      }

      void solve(Solver::SolverBase<OperatorType, VectorType, ScalarType>& solver) override
      {
         // Assemble the system
         if (!m_assembled)
            assemble();

         // Solve the system AX = B
         auto& axb = this->getLinearSystem();
         auto& a = axb.getOperator();
         auto& x = axb.getSolution();
         auto& b = axb.getVector();
         solver.solve(a, x, b);

         // Recover solution
         getTrialFunction().getSolution().setData(x);
      }

      ProblemUV& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      virtual ProblemUV* copy() const noexcept override = 0;

    private:
      Boolean m_assembled;
      std::reference_wrapper<TrialFunctionType> m_trialFunction;
      std::reference_wrapper<TestFunctionType>  m_testFunction;
      ProblemBodyType m_pb;
  };

  template <class LinearSystem, class U, class V>
  class Problem<LinearSystem, U, V> : public ProblemUV<LinearSystem, U, V>
  {
    public:
      using LinearSystemType = LinearSystem;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      using ProblemBodyType = ProblemBody<OperatorType, VectorType, ScalarType>;

      using Parent = ProblemUV<LinearSystem, U, V>;

      constexpr
      Problem(U& u, V& v)
        : Parent(u, v)
      {}

      constexpr
      Problem(const Problem& other)
        : Parent(other),
          m_axb(other.m_axb)
      {}

      constexpr
      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_axb(std::move(other.m_axb))
      {}

      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_axb = other.m_axb;
        }
        return *this;
      }

      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_axb = std::move(other.m_axb);
        }
        return *this;
      }

      Problem& operator=(const ProblemBodyType& rhs) override
      {
        Parent::operator=(rhs);
        return *this;
      }

      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      Problem* copy() const noexcept override
      {
        return new Problem(*this);
      }

    private:
      LinearSystemType m_axb;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class U, class V>
  Problem(U& u, V& v)
    -> Problem<
          Math::LinearSystem<
            Math::SparseMatrix<
              typename FormLanguage::Mult<
                typename FormLanguage::Traits<typename FormLanguage::Traits<U>::FESType>::ScalarType,
                typename FormLanguage::Traits<typename FormLanguage::Traits<V>::FESType>::ScalarType>::Type>,
            Math::Vector<
              typename FormLanguage::Traits<typename FormLanguage::Traits<V>::FESType>::ScalarType>>,
          U, V>;

  template <class LinearSystem, class U1, class U2, class ... Us>
  class ProblemUs
    : public ProblemBase<
        typename FormLanguage::Traits<LinearSystem>::OperatorType,
        typename FormLanguage::Traits<LinearSystem>::VectorType,
        typename FormLanguage::Traits<LinearSystem>::ScalarType>
  {
    template <class T>
    struct IsTrialOrTestFunction
    {
      static constexpr Boolean Value = IsTrialFunction<T>::Value || IsTestFunction<T>::Value;
    };

    static_assert(
        Utility::ParameterPack<U1, U2, Us...>::template All<IsTrialOrTestFunction>::Value);

    public:
      using LinearSystemType = LinearSystem;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using ProblemBodyType =
        ProblemBody<OperatorType, VectorType, ScalarType>;

      using Parent = ProblemBase<OperatorType, VectorType, ScalarType>;

    private:
      template <class T>
      struct GetFES;

      template <class T>
      struct GetFES<std::reference_wrapper<T>>
      {
        using Type = typename FormLanguage::Traits<T>::FESType;
      };

      template <class T>
      struct GetSolution;

      template <class T>
      struct GetSolution<std::reference_wrapper<T>>
      {
        using Type = typename FormLanguage::Traits<T>::SolutionType;
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

      using SolutionTuple = typename Utility::Extract<TrialFunctionTuple>::template Type<GetSolution>;

      using TrialFESTuple = typename Utility::Extract<TrialFunctionTuple>::template Type<GetFES>;

      using TestFESTuple = typename Utility::Extract<TestFunctionTuple>::template Type<GetFES>;

      template<class U, class V>
      using BilinearFormType =
        BilinearForm<
          typename FormLanguage::Traits<
            std::decay_t<typename Utility::UnwrapRefDecay<U>::Type>
          >::SolutionType,
          typename FormLanguage::Traits<
            std::decay_t<typename Utility::UnwrapRefDecay<U>::Type>
          >::FESType,
          typename FormLanguage::Traits<
            std::decay_t<typename Utility::UnwrapRefDecay<V>::Type>
          >::FESType,
          OperatorType
        >;

      template <class TestFES>
      using LinearFormType = LinearForm<TestFES, VectorType>;

      using BilinearFormTuple =
        typename Utility::Product<TrialFunctionTuple, TestFunctionTuple>::template Type<BilinearFormType>;

      using LinearFormTuple =
        typename Utility::Wrap<TestFESTuple>::template Type<LinearFormType>;

      using BilinearFormTupleSequentialAssembly =
        Assembly::Sequential<OperatorType, BilinearFormTuple>;

      using LinearFormTupleSequentialAssembly =
        Assembly::Sequential<VectorType, LinearFormTuple>;

    public:
      ProblemUs(U1& u1, U2& u2, Us&... us)
        : m_assembled(false),
          m_us(
            Tuple{std::ref(u1), std::ref(u2), std::ref(us)...}
            .template filter<IsTrialFunctionReferenceWrapper>()),
          m_vs(
            Tuple{std::ref(u1), std::ref(u2), std::ref(us)...}
            .template filter<IsTestFunctionReferenceWrapper>()),
          m_lft(m_vs.map(
                [](const auto& v)
                { return LinearFormType<
                    typename std::decay_t<
                    typename Utility::UnwrapRefDecay<decltype(v)>::Type>::FESType>(v.get());
                })),
          m_bft(m_us.product([](const auto& u, const auto& v) { return Pair(u, v); }, m_vs)
                    .map([](const auto& uv)
                         { return BilinearFormType<
                             decltype(uv.first()), decltype(uv.second())>(
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

      ProblemUs(const ProblemUs& other)
        : Parent(other),
          m_assembled(other.m_assembled),
          m_us(other.m_us),
          m_vs(other.m_vs),
          m_lft(other.m_lft),
          m_bft(other.m_bft),
          m_trialOffsets(other.m_trialOffsets),
          m_testOffsets(other.m_testOffsets),
          m_trialUUIDMap(other.m_trialUUIDMap),
          m_testUUIDMap(other.m_testUUIDMap),
          m_bfa(other.m_bfa->copy()),
          m_lfa(other.m_lfa->copy())
      {}

      ProblemUs(ProblemUs&& other) noexcept
        : Parent(std::move(other)),
          m_assembled(std::exchange(other.m_assembled, false)),
          m_us(std::move(other.m_us)),
          m_vs(std::move(other.m_vs)),
          m_lft(std::move(other.m_lft)),
          m_bft(std::move(other.m_bft)),
          m_trialOffsets(std::move(other.m_trialOffsets)),
          m_testOffsets(std::move(other.m_testOffsets)),
          m_trialUUIDMap(std::move(other.m_trialUUIDMap)),
          m_testUUIDMap(std::move(other.m_testUUIDMap)),
          m_bfa(std::move(other.m_bfa)),
          m_lfa(std::move(other.m_lfa))
      {}

      ProblemUs& operator=(const ProblemUs& other)
      {
        if (this != &other)
        {
          m_assembled = other.m_assembled;
          m_us = other.m_us;
          m_vs = other.m_vs;
          m_lft = other.m_lft;
          m_bft = other.m_bft;
          m_trialOffsets = other.m_trialOffsets;
          m_testOffsets = other.m_testOffsets;
          m_trialUUIDMap = other.m_trialUUIDMap;
          m_testUUIDMap = other.m_testUUIDMap;
          m_bfa.reset(other.m_bfa->copy());
          m_lfa.reset(other.m_lfa->copy());
        }
        return *this;
      }

      ProblemUs& operator=(ProblemUs&& other)
      {
        if (this != &other)
        {
          m_assembled = std::exchange(other.m_assembled, false);
          m_us = std::move(other.m_us);
          m_vs = std::move(other.m_vs);
          m_lft = std::move(other.m_lft);
          m_bft = std::move(other.m_bft);
          m_trialOffsets = std::move(other.m_trialOffsets);
          m_testOffsets = std::move(other.m_testOffsets);
          m_trialUUIDMap = std::move(other.m_trialUUIDMap);
          m_testUUIDMap = std::move(other.m_testUUIDMap);
          m_bfa.reset(std::move(other.m_bfa));
          m_lfa.reset(std::move(other.m_lfa));
        }
        return *this;
      }

      ProblemUs& assemble() override
      {
        auto& axb = getLinearSystem();

        m_lft.apply([&](auto& lf) { lf.clear(); });
        m_bft.apply([&](auto& bf) { bf.clear(); });

        for (auto& bfi : m_pb.getLocalBFIs())
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

        for (auto& bfi : m_pb.getGlobalBFIs())
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

        for (auto& lfi : m_pb.getLFIs())
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

        auto lt =
          m_lft.map(
              [](auto& lf)
              {
                auto& v = lf.getTestFunction();
                return Assembly::LinearFormAssemblyInput(
                    v.getFiniteElementSpace(), lf.getIntegrators());
              });

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

        m_lft.iapply(
            [&](const Index i, const auto& lf)
            {
              auto vi = m_testUUIDMap.left.find(lf.getTestFunction().getUUID());
              if (vi != m_testUUIDMap.left.end())
                loffsets[i] = m_testOffsets[vi->second];
            });

        m_bft.iapply(
            [&](const Index i, const auto& bf)
            {
              auto ui = m_trialUUIDMap.left.find(bf.getTrialFunction().getUUID());
              auto vi = m_testUUIDMap.left.find(bf.getTestFunction().getUUID());
              if (ui != m_trialUUIDMap.left.end() && vi != m_testUUIDMap.left.end())
                boffsets[i] = Pair{ m_trialOffsets[ui->second], m_testOffsets[vi->second] };
            });

        // Assemble stiffness operator
        m_bfa->execute(axb.getOperator(),
          Assembly::BilinearFormTupleAssemblyInput(rows, cols, boffsets, bt));

        // Assemble mass vector
        m_lfa->execute(axb.getVector(),
          Assembly::LinearFormTupleAssemblyInput(rows, loffsets, lt));

        // Impose Dirichlet boundary conditions
        m_us.apply(
            [&](const auto& u)
            {
              const auto ui = m_trialUUIDMap.left.find(u.get().getUUID());
              size_t offset = m_trialOffsets[ui->second];
              for (auto& dbc : m_pb.getDBCs())
              {
                if (dbc.getOperand().getUUID() == u.get().getUUID())
                {
                  dbc.assemble();
                  const auto& dofs = dbc.getDOFs();
                  axb.eliminate(dofs, offset);
                }
              }
            });

        m_assembled = true;

        return *this;
      }

      void solve(Solver::SolverBase<OperatorType, VectorType, ScalarType>& solver) override
      {
        auto& axb = getLinearSystem();

        // Assemble the system
        if (!m_assembled)
           assemble();

        // Solve the system AX = B
        solver.solve(axb.getOperator(), axb.getSolution(), axb.getVector());

        // Recover solutions
        m_us.iapply(
            [&](size_t i, auto& u)
            {
              u.get().getSolution().setData(axb.getSolution(), m_trialOffsets[i]);
            });
      }

      ProblemUs& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
        return *this;
      }

      const auto& getTrialOffsets() const
      {
        return m_trialOffsets;
      }

      const auto& getTestOffsets() const
      {
        return m_testOffsets;
      }

      virtual LinearSystemType& getLinearSystem() override = 0;

      virtual const LinearSystemType& getLinearSystem() const override = 0;

      virtual ProblemUs* copy() const noexcept override = 0;

    private:
      Boolean m_assembled;

      TrialFunctionTuple  m_us;
      TestFunctionTuple   m_vs;

      LinearFormTuple     m_lft;
      BilinearFormTuple   m_bft;

      ProblemBody<OperatorType, VectorType, ScalarType> m_pb;

      std::array<size_t, TrialFunctionTuple::Size>  m_trialOffsets;
      std::array<size_t, TestFunctionTuple::Size>   m_testOffsets;

      boost::bimap<FormLanguage::Base::UUID, size_t> m_trialUUIDMap;
      boost::bimap<FormLanguage::Base::UUID, size_t> m_testUUIDMap;

      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearFormTuple>>  m_bfa;
      std::unique_ptr<Assembly::AssemblyBase<VectorType, LinearFormTuple>>      m_lfa;
  };

  template <class LinearSystem, class U1, class U2, class ... Us>
  class Problem<LinearSystem, U1, U2, Us...> : public ProblemUs<LinearSystem, U1, U2, Us...>
  {
    public:
      using LinearSystemType = LinearSystem;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystemType>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystemType>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystemType>::ScalarType;

      using ProblemBodyType =
        ProblemBody<OperatorType, VectorType, ScalarType>;

      using Parent = ProblemUs<LinearSystem, U1, U2, Us...>;

      Problem(U1& u1, U2& u2, Us&... us)
        : Parent(u1, u2, us...)
      {}

      Problem(const Problem& other)
        : Parent(other),
          m_axb(other.m_axb)
      {}

      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_axb(std::move(other.m_axb))
      {}

      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_axb = other.m_axb;
        }
        return *this;
      }

      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_axb = std::move(other.m_axb);
        }
        return *this;
      }

      Problem& operator=(const ProblemBodyType& rhs) override
      {
        Parent::operator=(rhs);
        return *this;
      }

      LinearSystemType& getLinearSystem() override
      {
        return m_axb;
      }

      const LinearSystemType& getLinearSystem() const override
      {
        return m_axb;
      }

      Problem* copy() const noexcept override
      {
        return new Problem(*this);
      }

    private:
      LinearSystemType m_axb;
  };

  template <class U1, class U2, class ... Us>
  Problem(U1& u1, U2& u2, Us&... us)
    -> Problem<
        Math::LinearSystem<
          Math::SparseMatrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<typename FormLanguage::Traits<U1>::FESType>::ScalarType,
              typename FormLanguage::Traits<typename FormLanguage::Traits<U2>::FESType>::ScalarType>::Type>,
          Math::Vector<
            typename FormLanguage::Traits<typename FormLanguage::Traits<U2>::FESType>::ScalarType>>,
        U1, U2, Us...>;
}

#endif
