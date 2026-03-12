/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Problem.h
 * @brief Variational problem classes for finite element systems.
 *
 * This file defines the Problem classes which encapsulate complete finite element
 * problems including bilinear forms, linear forms, boundary conditions, and
 * solution methods. Problems represent the discrete system @f$ Au = b @f$ arising
 * from weak formulations.
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_H
#define RODIN_VARIATIONAL_PROBLEM_H

#include <functional>
#include <boost/mp11.hpp>

#include "Rodin/Pair.h"

#include "Rodin/Math/ForwardDecls.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/LinearSystem.h"
#include "Rodin/Math/SparseMatrix.h"

#include "Rodin/Utility/Extract.h"
#include "Rodin/Utility/Product.h"
#include "Rodin/Utility/Wrap.h"

#include "Rodin/Assembly/ForwardDecls.h"
#include "Rodin/Assembly/Input.h"
#include "Rodin/Assembly/Sequential.h"

#include "Rodin/Solver/LinearSolver.h"

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ProblemBody.h"

namespace Rodin::Variational
{
  /**
   * @defgroup RodinVariational Variational Module
   * @brief Variational formulations and finite element problem definitions.
   *
   * The Variational module provides the infrastructure for defining and solving
   * variational problems in finite element analysis. It supports the formulation
   * of linear and nonlinear problems, boundary conditions, and various finite
   * element spaces.
   *
   * ## Mathematical Foundation
   * The module implements the standard finite element methodology:
   * - **Bilinear Forms**: @f$ a(u,v) : V \times V \to \mathbb{R} @f$
   * - **Linear Forms**: @f$ l(v) : V \to \mathbb{R} @f$ 
   * - **Weak Formulation**: Find @f$ u \in V @f$ such that @f$ a(u,v) = l(v) @f$ for all @f$ v \in V @f$
   * - **Discrete System**: @f$ Au = b @f$ where @f$ A_{ij} = a(\phi_j, \psi_i) @f$ and @f$ b_i = l(\psi_i) @f$
   */

  /**
   * @defgroup ProblemSpecializations Problem Template Specializations
   * @brief Template specializations of the Problem class.
   * @see Problem
   */

  /**
   * @ingroup RodinVariational
   * @brief Abstract base class for variational problems.
   *
   * This class provides the foundation for defining and solving variational
   * problems in finite element analysis. It encapsulates the discrete linear
   * system @f$ Au = b @f$ arising from the weak formulation of PDEs and provides
   * a unified interface for assembly and solution operations.
   *
   * @tparam LinearSystem Type of linear system (e.g., sparse or dense)
   *
   * ## Usage Example
   * ```cpp
   * // Define a problem and solve
   * auto problem = Problem(trialSpace, testSpace);
   * problem = bilinearForm == linearForm;
   * problem.assemble();
   * problem.solve(solver);
   * ```
   */
  template <class LinearSystem>
  class ProblemBase : public FormLanguage::Base
  {
    public:
      /// @brief Matrix operator type from linear system traits
      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      /// @brief Vector type from linear system traits
      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      /// @brief Scalar type from linear system traits  
      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      /// @brief Problem body type containing bilinear and linear forms
      using ProblemBodyType =
        ProblemBody<OperatorType, VectorType, ScalarType>;

      /// @brief Default constructor
      ProblemBase() = default;

      /// @brief Move constructor
      ProblemBase(ProblemBase&& other) = default;

      /// @brief Copy constructor
      ProblemBase(const ProblemBase& other) = default;

      /**
       * @brief Assigns a problem body to this problem.
       *
       * @param rhs Problem body containing bilinear and linear forms
       * @return Reference to this problem instance
       *
       * This operator allows setting up the variational problem by assigning
       * a problem body that contains the bilinear form @f$ a(u,v) @f$ and
       * linear form @f$ l(v) @f$ defining the weak formulation.
       */
      virtual ProblemBase& operator=(const ProblemBodyType& rhs) = 0;

      /**
       * @brief Solves the assembled linear system.
       *
       * @param solver Solver instance to use for solving the linear system
       *
       * Solves the discrete system @f$ Au = b @f$ using the provided solver.
       * The problem must be assembled before calling this method.
       */
      virtual void solve(Solver::LinearSolverBase<LinearSystem>& solver) = 0;

      /**
       * @brief Assembles the underlying linear system to solve.
       *
       * @return Reference to this problem instance for method chaining
       *
       * Assembles the bilinear form into the system matrix @f$ A @f$ and the
       * linear form into the right-hand side vector @f$ b @f$, creating the
       * discrete linear system @f$ Au = b @f$.
       */
      virtual ProblemBase& assemble() = 0;

      /**
       * @brief Gets the assembled linear system.
       * @return Reference to the linear system @f$ Au = b @f$
       *
       * Provides access to the assembled discrete system consisting of the
       * system matrix @f$ A @f$ and right-hand side vector @f$ b @f$.
       */
      virtual LinearSystem& getLinearSystem() = 0;

      /**
       * @brief Gets the assembled linear system (const version).
       * @return Const reference to the linear system @f$ Au = b @f$
       */
      virtual const LinearSystem& getLinearSystem() const = 0;

      /**
       * @brief Creates a polymorphic copy of this problem.
       * @return Pointer to a new copy of this problem
       */
      virtual ProblemBase* copy() const noexcept override = 0;
  };

  /**
   * @ingroup ProblemSpecializations
   * @brief General class to assemble linear systems with `Operator`
   * and `Vector` generic types in a sequential context.
   */
  template <class LinearSystem, class U, class V>
  class ProblemUVBase : public ProblemBase<LinearSystem>
  {
    public:
      using TrialFunctionType =
        U;

      using TestFunctionType =
        V;

      using LinearSystemType =
        LinearSystem;

      using SolverBaseType =
        Solver::LinearSolverBase<LinearSystem>;

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

      using TrialFESScalarType =
        typename FormLanguage::Traits<TrialFESType>::ScalarType;

      using TestFESScalarType =
        typename FormLanguage::Traits<TestFESType>::ScalarType;

      using LinearFormIntegratorBaseType =
        LinearFormIntegratorBase<TestFESScalarType>;

      using ProblemBodyType =
        ProblemBody<OperatorType, VectorType, ScalarType>;

      using Parent =
        ProblemBase<LinearSystemType>;

      constexpr
      ProblemUVBase(U& u, V& v)
        : m_trialFunction(u), m_testFunction(v)
      {}

      ProblemUVBase(const ProblemUVBase& other)
        : Parent(other),
          m_trialFunction(other.m_trialFunction),
          m_testFunction(other.m_testFunction)
      {}

      ProblemUVBase(ProblemUVBase&& other)
        : Parent(std::move(other)),
          m_trialFunction(std::move(other.m_trialFunction)),
          m_testFunction(std::move(other.m_testFunction))
      {}

      ProblemUVBase& operator=(const ProblemUVBase& other)
      {
        if (this != &other)
        {
          m_trialFunction = other.m_trialFunction;
          m_testFunction = other.m_testFunction;
        }
        return *this;
      }

      ProblemUVBase& operator=(ProblemUVBase&& other) noexcept
      {
        if (this != &other)
        {
          m_trialFunction = std::move(other.m_trialFunction);
          m_testFunction = std::move(other.m_testFunction);
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

      virtual ProblemUVBase& assemble() override = 0;

      virtual void solve(SolverBaseType& solver) override = 0;

      virtual ProblemUVBase& operator=(const ProblemBodyType& rhs) override = 0;

      virtual ProblemUVBase* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<TrialFunctionType> m_trialFunction;
      std::reference_wrapper<TestFunctionType> m_testFunction;
  };

  template <class LinearSystem, class TrialFunction, class TestFunction>
  class Problem<LinearSystem, TrialFunction, TestFunction>
    : public ProblemUVBase<LinearSystem, TrialFunction, TestFunction>
  {
    using TrialFESType =
      typename FormLanguage::Traits<TrialFunction>::FESType;

    using TestFESType =
      typename FormLanguage::Traits<TestFunction>::FESType;

    using TrialFESMeshType =
      typename FormLanguage::Traits<TrialFESType>::MeshType;

    using TestFESMeshType =
      typename FormLanguage::Traits<TestFESType>::MeshType;

    using TrialFESContextType =
      typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

    using TestFESContextType =
      typename FormLanguage::Traits<TestFESMeshType>::ContextType;

    public:
      using LinearSystemType =
        LinearSystem;

      using SolverBaseType =
        Solver::LinearSolverBase<LinearSystemType>;

      using OperatorType =
        typename FormLanguage::Traits<LinearSystem>::OperatorType;

      using VectorType =
        typename FormLanguage::Traits<LinearSystem>::VectorType;

      using ScalarType =
        typename FormLanguage::Traits<LinearSystem>::ScalarType;

      using ProblemBodyType =
        ProblemBody<OperatorType, VectorType, ScalarType>;

      using TrialFESMeshContextType =
        typename FormLanguage::Traits<TrialFESMeshType>::ContextType;

      using TestFESMeshContextType =
        typename FormLanguage::Traits<TestFESMeshType>::ContextType;

      using AssemblyType =
        typename Assembly::Default<TrialFESMeshContextType, TestFESMeshContextType>
          ::template Type<LinearSystem, Problem>;

      using Parent =
        ProblemUVBase<LinearSystem, TrialFunction, TestFunction>;

      constexpr
      Problem(TrialFunction& u, TestFunction& v)
        : Parent(u, v)
      {}

      constexpr
      Problem(const Problem& other)
        : Parent(other),
          m_assembled(other.m_assembled),
          m_pb(other.m_pb),
          m_axb(other.m_axb)
      {}

      constexpr
      Problem(Problem&& other) noexcept
        : Parent(std::move(other)),
          m_assembled(other.m_assembled),
          m_pb(std::move(other.m_pb)),
          m_axb(std::move(other.m_axb)),
          m_assembly(std::move(other.m_assembly))
      {}

      Problem& operator=(const Problem& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_assembled = other.m_assembled;
          m_pb = other.m_pb;
          m_axb = other.m_axb;
          m_assembly.reset(other.m_assembly->copy());
        }
        return *this;
      }

      Problem& operator=(Problem&& other) noexcept
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_assembled = std::exchange(other.m_assembled, false);
          m_pb = std::move(other.m_pb);
          m_axb = std::move(other.m_axb);
          m_assembly = std::move(other.m_assembly);
        }
        return *this;
      }

      Problem& assemble() override
      {
        m_assembly.execute(m_axb, { m_pb, this->getTrialFunction(), this->getTestFunction() });
        m_assembled = true;
        return *this;
      }

      void solve(SolverBaseType& solver) override
      {
         auto& axb = this->getLinearSystem();
         if (!m_assembled)
            this->assemble();
         solver.solve(axb);
         this->getTrialFunction().getSolution().setData(axb.getSolution());
      }

      Problem& operator=(const ProblemBodyType& rhs) override
      {
        m_pb = rhs;
        m_assembled = false;
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
      Boolean m_assembled;
      ProblemBodyType m_pb;
      LinearSystemType m_axb;
      AssemblyType m_assembly;
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

  template <class LinearSystem, class U1, class U2, class U3, class ... Us>
  class ProblemUsBase : public ProblemBase<LinearSystem>
  {
    template <class T>
    struct IsTrialOrTestFunction
    {
      static constexpr Boolean Value = IsTrialFunction<T>::Value || IsTestFunction<T>::Value;
    };

    static_assert(
        Utility::ParameterPack<U1, U2, U3, Us...>::template All<IsTrialOrTestFunction>::Value);

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

      using Parent = ProblemBase<LinearSystemType>;

      using AssemblyInput =
        Assembly::ProblemAssemblyInput<ProblemBodyType, U1, U2, U3, Us...>;

      using U1FESType = typename FormLanguage::Traits<U1>::FESType;
      using U2FESType = typename FormLanguage::Traits<U2>::FESType;
      using U3FESType = typename FormLanguage::Traits<U3>::FESType;

      using U1FESMeshType = typename FormLanguage::Traits<U1FESType>::MeshType;
      using U2FESMeshType = typename FormLanguage::Traits<U2FESType>::MeshType;
      using U3FESMeshType = typename FormLanguage::Traits<U3FESType>::MeshType;

      using U1FESMeshContextType = typename FormLanguage::Traits<U1FESMeshType>::ContextType;
      using U2FESMeshContextType = typename FormLanguage::Traits<U2FESMeshType>::ContextType;
      using U3FESMeshContextType = typename FormLanguage::Traits<U3FESMeshType>::ContextType;

      using AssemblyType =
        typename Assembly::Default<U1FESMeshContextType, U2FESMeshContextType>
          ::template Type<LinearSystemType, Problem<LinearSystemType, U1, U2, U3, Us...>>;

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
            std::reference_wrapper<U3>,
            std::reference_wrapper<Us>...>>()
            .template filter<IsTrialFunctionReferenceWrapper>());

      using TestFunctionTuple =
        decltype(std::declval<
          Tuple<
            std::reference_wrapper<U1>,
            std::reference_wrapper<U2>,
            std::reference_wrapper<U3>,
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
      ProblemUsBase(U1& u1, U2& u2, U3& u3, Us&... us)
        : m_assembled(false),
          m_us(
            Tuple{std::ref(u1), std::ref(u2), std::ref(u3), std::ref(us)...}
            .template filter<IsTrialFunctionReferenceWrapper>()),
          m_vs(
            Tuple{std::ref(u1), std::ref(u2), std::ref(u3), std::ref(us)...}
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

      ProblemUsBase(const ProblemUsBase& other)
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
          m_lfa(other.m_lfa->copy()),
          m_assembly(other.m_assembly)
      {}

      ProblemUsBase(ProblemUsBase&& other) noexcept
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
          m_lfa(std::move(other.m_lfa)),
          m_assembly(std::move(other.m_assembly))
      {}

      ProblemUsBase& operator=(const ProblemUsBase& other)
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
          m_assembly = other.m_assembly;
        }
        return *this;
      }

      ProblemUsBase& operator=(ProblemUsBase&& other)
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
          m_assembly = std::move(other.m_assembly);
        }
        return *this;
      }

      virtual ProblemUsBase& assemble() override
      {
        auto& axb = getLinearSystem();

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
          for (size_t i = 0; i < TestFunctionTuple::Size - 1; i++)
            m_testOffsets[i + 1] = sz[i] + m_testOffsets[i];
        }

        size_t rows =
          m_vs
            .map([](const auto& v)
            {
              return static_cast<size_t>(v.get().getFiniteElementSpace().getSize());
            })
            .reduce([](size_t a, size_t b) { return a + b; });

        size_t cols =
          m_us
            .map([](const auto& u)
            {
              return static_cast<size_t>(u.get().getFiniteElementSpace().getSize());
            })
            .reduce([](size_t a, size_t b) { return a + b; });

        AssemblyInput input(
            m_pb, m_us, m_vs, m_trialOffsets, m_testOffsets,
            m_trialUUIDMap, m_testUUIDMap, cols, rows);
        m_assembly.execute(axb, input);

        m_assembled = true;

        return *this;
      }

      void solve(Solver::LinearSolverBase<LinearSystemType>& solver) override
      {
        auto& axb = getLinearSystem();
        if (!m_assembled)
           this->assemble();
        solver.solve(axb);
        m_us.iapply(
            [&](size_t i, auto& u)
            {
              u.get().getSolution().setData(axb.getSolution(), m_trialOffsets[i]);
            });
      }

      ProblemUsBase& operator=(const ProblemBodyType& rhs) override
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

      virtual ProblemUsBase* copy() const noexcept override = 0;

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
      AssemblyType m_assembly;
  };

  template <class LinearSystem, class U1, class U2, class U3, class ... Us>
  class Problem<LinearSystem, U1, U2, U3, Us...>
    : public ProblemUsBase<LinearSystem, U1, U2, U3, Us...>
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

      using Parent = ProblemUsBase<LinearSystem, U1, U2, U3, Us...>;

      Problem(U1& u1, U2& u2, U3& u3, Us&... us)
        : Parent(u1, u2, u3, us...)
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

  template <class U1, class U2, class U3, class ... Us>
  Problem(U1& u1, U2& u2, U3& u3, Us&... us)
    -> Problem<
        Math::LinearSystem<
          Math::SparseMatrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<typename FormLanguage::Traits<U1>::FESType>::ScalarType,
              typename FormLanguage::Traits<typename FormLanguage::Traits<U2>::FESType>::ScalarType>::Type>,
          Math::Vector<
            typename FormLanguage::Traits<typename FormLanguage::Traits<U2>::FESType>::ScalarType>>,
        U1, U2, U3, Us...>;
}

#endif
