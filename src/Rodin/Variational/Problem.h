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
  template <class OperatorType, class VectorType>
  class ProblemBase : public FormLanguage::Base
  {
    public:
      ProblemBase() = default;

      ProblemBase(ProblemBase&& other) = default;

      ProblemBase(const ProblemBase& other) = default;

      virtual ProblemBase& operator=(const ProblemBody<OperatorType, VectorType>& rhs) = 0;

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

  template <>
  class Problem<> {};

  /**
   * @ingroup ProblemSpecializations
   * @brief General class to assemble linear systems with `Math::SparseMatrix`
   * and `Math::Vector` types in a serial context.
   */
  template <class TrialFES, class TestFES>
  class Problem<TrialFES, TestFES, Context::Sequential, Math::SparseMatrix, Math::Vector>
    : public ProblemBase<Math::SparseMatrix, Math::Vector>
  {
      static_assert(std::is_same_v<typename TrialFES::Context, Context::Sequential>);
      static_assert(std::is_same_v<typename TestFES::Context, Context::Sequential>);

    public:
      using Context = Context::Sequential;
      using OperatorType = Math::SparseMatrix;
      using VectorType = Math::Vector;
      using Parent = ProblemBase<Math::SparseMatrix, Math::Vector>;

      /**
       * @brief Constructs an empty problem involving the trial function @f$ u @f$
       * and the test function @f$ v @f$.
       *
       * @param[in,out] u Trial function
       * @param[in,out] v %Test function
       */
      explicit
      constexpr
      Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v);

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

      Problem& imposePeriodicBCs();

      Problem& imposeDirichletBCs();

      Problem& operator+=(const LocalBilinearFormIntegratorBase& lfi);

      Problem& operator-=(const LocalBilinearFormIntegratorBase& lfi);

      Problem& operator+=(const LinearFormIntegratorBase& lfi);

      Problem& operator-=(const LinearFormIntegratorBase& lfi);

      Problem& operator+=(const DirichletBCBase& dbc);

      Problem& operator+=(const PeriodicBCBase& pbc);

      Problem& assemble() override;

      void solve(Solver::SolverBase<OperatorType, VectorType>& solver) override;

      Problem& operator=(const ProblemBody<OperatorType, VectorType>& rhs) override;

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
      BilinearForm<TrialFES, TestFES, Math::SparseMatrix> m_bilinearForm;

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
    -> Problem<TrialFES, TestFES, typename TrialFES::Context, Math::SparseMatrix, Math::Vector>;

  template <class U1, class U2, class ... Us>
  class Problem<
      Tuple<U1, U2, Us...>,
      Context::Sequential, Math::SparseMatrix, Math::Vector>
    : public ProblemBase<Math::SparseMatrix, Math::Vector>
  {

    template <class T>
    struct IsTrialOrTestFunction
    {
      static constexpr bool Value = IsTrialFunction<T>::Value || IsTestFunction<T>::Value;
    };

    static_assert(Utility::ParameterPack<U1, U2, Us...>::template All<IsTrialOrTestFunction>::Value);

    public:
      using Context = Context::Sequential;
      using OperatorType = Math::SparseMatrix;
      using VectorType = Math::Vector;
      using Parent = ProblemBase<Math::SparseMatrix, Math::Vector>;

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
      using SequentialAssembly = Assembly::Sequential<OperatorType, BilinearFormTuple>;

      Problem(U1& u1, U2& u2, Us&... us);

      Problem& assemble() override;

      void solve(Solver::SolverBase<OperatorType, VectorType>& solver) override
      {}

      Problem& operator=(const ProblemBody<OperatorType, VectorType>& rhs) override
      {
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
      TrialFunctionTuple m_us;
      TestFunctionTuple  m_vs;

      LinearFormTuple   m_lft;
      BilinearFormTuple m_bft;

      bool            m_assembled;
      VectorType      m_mass;
      VectorType      m_guess;
      OperatorType    m_stiffness;

      std::unique_ptr<Assembly::AssemblyBase<OperatorType, BilinearFormTuple>> m_assembly;
  };

  template <class U1, class U2, class ... Us>
  Problem(U1& u1, U2& u2, Us&... us)
    -> Problem<
        Tuple<U1, U2, Us...>,
        Context::Sequential, Math::SparseMatrix, Math::Vector>;
}

#include "Problem.hpp"

#endif
