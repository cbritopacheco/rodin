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
    -> DenseProblem<TrialFES, TestFES, typename TrialFES::Context, Math::Matrix<Scalar>, Math::Vector<Scalar>>;

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
  class DenseProblem<TrialFES, TestFES, Context::Sequential, Math::Matrix<Scalar>, Math::Vector<Scalar>>
    : public ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>>
  {
      static_assert(std::is_same_v<typename TrialFES::Context, Context::Sequential>);
      static_assert(std::is_same_v<typename TestFES::Context, Context::Sequential>);

    public:
      using Context = Context::Sequential;
      using OperatorType = Math::Matrix<Scalar>;
      using VectorType = Math::Vector<Scalar>;
      using Parent = ProblemBase<Math::Matrix<Scalar>, Math::Vector<Scalar>>;

      /**
       * @brief Constructs an empty DenseProblem involving the trial function @f$ u @f$
       * and the test function @f$ v @f$.
       *
       * @param[in,out] u Trial function
       * @param[in,out] v %Test function
       */
      explicit
      constexpr
      DenseProblem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v);

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

      DenseProblem& imposePeriodicBCs();

      DenseProblem& imposeDirichletBCs();

      DenseProblem& operator+=(const LocalBilinearFormIntegratorBase& lfi);

      DenseProblem& operator-=(const LocalBilinearFormIntegratorBase& lfi);

      DenseProblem& operator+=(const LinearFormIntegratorBase& lfi);

      DenseProblem& operator-=(const LinearFormIntegratorBase& lfi);

      DenseProblem& operator+=(const DirichletBCBase& dbc);

      DenseProblem& operator+=(const PeriodicBCBase& pbc);

      DenseProblem& assemble() override;

      void solve(Solver::SolverBase<OperatorType, VectorType>& solver) override;

      DenseProblem& operator=(const ProblemBody<OperatorType, VectorType>& rhs) override;

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
      BilinearForm<TrialFES, TestFES, Math::Matrix<Scalar>> m_bilinearForm;

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

