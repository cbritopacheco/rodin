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
          * @brief Builds the problem using the FormLanguage::ProblemBody
          * object.
          *
          * @param[in] rhs Problem body constructed using expressions from the
          * Variational::FormLanguage.
          */
         virtual ProblemBase& operator=(const ProblemBody& rhs) = 0;

         /**
          * @brief Assembles the underlying linear system to solve.
          */
         virtual void assemble() = 0;

         /**
          * @brief Updates the ProblemBase instance after a refinement in the mesh
          */
         virtual ProblemBase& update() = 0;

         virtual void recoverSolution() = 0;

         virtual mfem::Operator& getStiffnessMatrix() = 0;

         virtual const mfem::Operator& getStiffnessMatrix() const = 0;

         virtual mfem::Vector& getMassVector() = 0;

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
   template <class TrialFEC, class TestFEC, class OperatorType>
   class Problem<TrialFEC, TestFEC, OperatorType, Traits::Serial>
      : public ProblemBase
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
         Problem(
               TrialFunction<TrialFEC, Traits::Serial>& u,
               TestFunction<TestFEC, Traits::Serial>& v,
               OperatorType* = new OperatorType);

         Problem& update() override;

         void assemble() override;

         void recoverSolution() override;

         EssentialBoundary& getEssentialBoundary() override;

         Problem& operator=(const ProblemBody& rhs) override;

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
         LinearForm<TrialFEC, Traits::Serial>                m_linearForm;
         BilinearForm<TrialFEC, TestFEC, Traits::Serial>     m_bilinearForm;

         mfem::OperatorHandle    m_stiffnessOp;
         mfem::Vector            m_massVector;
         mfem::Vector            m_guess;

         const std::map<
            boost::uuids::uuid,
            std::reference_wrapper<TrialFunction<TrialFEC, Traits::Serial>>> m_trialFunctions;
         const std::map<
            boost::uuids::uuid,
            std::reference_wrapper<TestFunction<TestFEC, Traits::Serial>>> m_testFunctions;
         std::unique_ptr<ProblemBody> m_pb;

         mfem::Array<int> m_essTrueDofList;
   };
   template <class TrialFEC, class TestFEC>
   Problem(TrialFunction<TrialFEC, Traits::Serial>&, TestFunction<TestFEC, Traits::Serial>&)
      -> Problem<TrialFEC, TestFEC, mfem::SparseMatrix, Traits::Serial>;

   template <class TrialFEC, class TestFEC, class OperatorType, class Trait>
   Problem(
         TrialFunction<TrialFEC, Trait>&,
         TestFunction<TestFEC, Trait>&,
         OperatorType*)
      -> Problem<TrialFEC, TestFEC, OperatorType, Trait>;
}

#include "Problem.hpp"

#endif
