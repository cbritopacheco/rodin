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

#include "ForwardDecls.h"

#include "ProblemBody.h"
#include "LinearForm.h"
#include "BilinearForm.h"
#include "TrialFunction.h"
#include "TestFunction.h"
#include "EssentialBoundary.h"


namespace Rodin::Variational
{
   /**
    * @brief Base class for Variational::Problem.
    */
   class ProblemBase
   {
      public:
         virtual EssentialBoundary& getEssentialBoundary() = 0;

         /**
          * @brief Assembles the underlying linear system to solve.
          */
         virtual void assemble() = 0;

         /**
          * @brief Updates the ProblemBase instance after a refinement in the mesh
          */
         virtual ProblemBase& update() = 0;

         /**
          * @brief After solving the problem using a linear solver, this method
          * must be called to recover the solution.
          *
          * After this call, the solution(s) will be contained in the
          * respective TrialFunction object(s) and can be obtained via
          * TrialFunction::getGridFunction().
          */
         virtual void recoverSolution() = 0;

         /**
          * @returns Reference to the mfem::Operator representing the stiffness
          * matrix.
          *
          * This must be called only after assemble() has been called.
          */
         virtual mfem::Operator& getStiffnessMatrix() = 0;

         /**
          * @returns Constant reference to the mfem::Operator representing the stiffness
          * matrix, i.e. the LHS of the weak formulation.
          *
          * This must be called only after assemble() has been called.
          */
         virtual const mfem::Operator& getStiffnessMatrix() const = 0;

         /**
          * @returns Reference to the mfem::Vector representing the mass
          * vector, i.e. the RHS of the weak formulation.
          *
          * This must be called only after assemble() has been called.
          */
         virtual mfem::Vector& getMassVector() = 0;

         /**
          * @returns Constant reference to the mfem::Vector representing the mass
          * vector, i.e. the RHS of the weak formulation.
          *
          * This must be called only after assemble() has been called.
          */
         virtual const mfem::Vector& getMassVector() const = 0;

         virtual mfem::Vector& getInitialGuess() = 0;

         virtual const mfem::Vector& getInitialGuess() const = 0;
   };

   /**
    * @brief Represents a variational problem to be solved.
    *
    * The problem may then be solved by utilizing any derived instance Solver
    * class in the Rodin::Solver namespace.
    *
    */
   template <class TrialFES, class TestFES, class OperatorType>
   class Problem<TrialFES, TestFES, OperatorType> : public ProblemBase
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
          * @param[in,out] u Trial function @f$ u @f$
          * @param[in,out] v Test function @f$ v @f$
          */
         explicit
         Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v, OperatorType&& = OperatorType());

         explicit
         Problem(TrialFunction<TrialFES>& u, TestFunction<TestFES>& v, OperatorType&);

         Problem& operator=(ProblemBody&& rhs);

         Problem& update() override;

         void assemble() override;

         void recoverSolution() override;

         EssentialBoundary& getEssentialBoundary() override;


         OperatorType& getStiffnessMatrix() override
         {
            return *m_stiffnessOp.As<OperatorType>();
         }

         const OperatorType& getStiffnessMatrix() const override
         {
            return *m_stiffnessOp.As<OperatorType>();
         }

         mfem::Vector& getMassVector() override
         {
            return m_massVector;
         }

         const mfem::Vector& getMassVector() const override
         {
            return m_massVector;
         }

         mfem::Vector& getInitialGuess() override
         {
            return m_guess;
         }

         const mfem::Vector& getInitialGuess() const override
         {
            return m_guess;
         }

      private:
         ProblemBody m_pb;

         LinearForm<TestFES>                m_linearForm;
         BilinearForm<TrialFES, TestFES>    m_bilinearForm;

         mfem::OperatorHandle    m_stiffnessOp;
         mfem::Vector            m_massVector;
         mfem::Vector            m_guess;

         const std::map<
            boost::uuids::uuid,
            std::reference_wrapper<TrialFunction<TrialFES>>> m_trialFunctions;
         const std::map<
            boost::uuids::uuid,
            std::reference_wrapper<TestFunction<TestFES>>> m_testFunctions;

         mfem::Array<int> m_essTrueDofList;
   };
   template <class TrialFES, class TestFES>
   Problem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
      -> Problem<TrialFES, TestFES, mfem::SparseMatrix>;

   template <class TrialFES, class TestFES, class OperatorType>
   Problem(TrialFunction<TrialFES>&, TestFunction<TestFES>&, OperatorType&&)
      -> Problem<TrialFES, TestFES, OperatorType>;

   template <class TrialFES, class TestFES, class OperatorType>
   Problem(TrialFunction<TrialFES>&, TestFunction<TestFES>&, OperatorType&)
      -> Problem<TrialFES, TestFES, OperatorType>;
}

#include "Problem.hpp"

#endif
