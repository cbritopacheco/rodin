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

#include <mfem.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Mesh.h"

#include "FormLanguage/ForwardDecls.h"

#include "TrialFunction.h"
#include "TestFunction.h"
#include "ForwardDecls.h"
#include "BilinearForm.h"
#include "LinearForm.h"


namespace Rodin::Variational
{
   /**
    * @brief Base class for Variational::Problem.
    */
   class ProblemBase
   {
      public:
         /**
          * @brief Gets a reference to the solution.
          * @returns Reference to the solution.
          */
         virtual GridFunctionBase& getSolution() = 0;

         virtual std::set<int>& getEssentialBoundary() = 0;

         virtual BilinearFormBase& getBilinearForm() = 0;

         virtual LinearFormBase& getLinearForm() = 0;

         /**
          * @brief Builds the problem using the FormLanguage::ProblemBody
          * object.
          *
          * @param[in] rhs Problem body constructed using expressions from the
          * Variational::FormLanguage.
          */
         virtual ProblemBase& operator=(const FormLanguage::ProblemBody& rhs) = 0;

         /**
          * @brief Assembles the underlying linear system to solve.
          *
          * @note This is typically the first thing that is done when calling
          * the @ref Solver::solve() method.
          */
         virtual void assemble() = 0;

         virtual void update() = 0;
   };

   /**
    * @brief Represents a variational problem to be solved.
    *
    * The problem may be specified via the overloaded operator
    * @ref Variational::Problem::operator=(const FormLanguage::ProblemBody&).
    *
    * The problem may then be solved by utilizing any derived instance Solver
    * class in the Rodin::Solver namespace.
    *
    * @note The underlying linear system is only assembled until the
    * @ref assemble() method is called. It is usually the responsibility of
    * the derived @ref Solver object to assemble the problem when solving it
    * via the @ref solve(Problem&) method.
    */
   template <class TrialFEC, class TestFEC>
   class Problem : public ProblemBase
   {
      public:
         /**
          * @brief Deleted copy constructor.
          */
         Problem(const Problem& other) = delete;

         /**
          * @brief Deleted copy assignment operator.
          */
         void operator=(const Problem& other) = delete;

         /**
          * @brief Constructs an empty problem involving the trial function @f$ u @f$
          * and the test function @f$ v @f$.
          *
          * @param[in,out] u Trial function @f$ u @f$ belonging to a suitable
          * finite element space.
          */
         Problem(TrialFunction<TrialFEC>& u, TestFunction<TestFEC>& v);

         void assemble() override;
         void update() override;

         GridFunction<TrialFEC>& getSolution() override;
         LinearForm<TrialFEC>& getLinearForm() override;
         BilinearForm<TrialFEC>& getBilinearForm() override;
         std::set<int>& getEssentialBoundary() override;

         Problem& operator=(const FormLanguage::ProblemBody& rhs) override;
      private:
         TrialFunction<TrialFEC>&   m_solution;
         LinearForm<TrialFEC>       m_linearForm;
         BilinearForm<TrialFEC>     m_bilinearForm;
         std::set<int> m_essBdr;
         std::unique_ptr<FormLanguage::ProblemBody> m_pb;
   };
}

#include "Problem.hpp"

#endif
