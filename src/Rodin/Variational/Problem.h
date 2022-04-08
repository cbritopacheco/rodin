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
         virtual void getLinearSystem(mfem::SparseMatrix& A, mfem::Vector& B, mfem::Vector& X) = 0;
         virtual void getSolution(mfem::Vector& X) = 0;

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
         virtual void update() = 0;
   };

   /**
    * @brief Represents a variational problem to be solved.
    *
    * The problem may be specified via the overloaded operator
    * @ref Variational::Problem::operator=(const ProblemBody&).
    *
    * The problem may then be solved by utilizing any derived instance Solver
    * class in the Rodin::Solver namespace.
    *
    */
   template <class TrialFEC, class TestFEC>
   class Problem<TrialFEC, TestFEC, Traits::Serial>
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
          * @param[in,out] u Trial function @f$ u @f$ belonging to a suitable
          * finite element space.
          */
         Problem(
               TrialFunction<TrialFEC, Traits::Serial>& u,
               TestFunction<TestFEC, Traits::Serial>& v);

         void assemble() override;
         void update() override;


         void getSolution(mfem::Vector& X) override;

         EssentialBoundary& getEssentialBoundary() override;

         void getLinearSystem(mfem::SparseMatrix& A, mfem::Vector& B, mfem::Vector& X) override;

         Problem& operator=(const ProblemBody& rhs) override;

      private:
         LinearForm<TrialFEC, Traits::Serial>                m_linearForm;
         BilinearForm<TrialFEC, TestFEC, Traits::Serial>     m_bilinearForm;

         const std::map<
            boost::uuids::uuid,
            std::reference_wrapper<TrialFunction<TrialFEC, Traits::Serial>>> m_trialFunctions;
         const std::map<
            boost::uuids::uuid,
            std::reference_wrapper<TestFunction<TestFEC, Traits::Serial>>> m_testFunctions;
         std::unique_ptr<ProblemBody> m_pb;

         mfem::Array<int> m_essTrueDofList;
   };
   template <class TrialFEC, class TestFEC, class Trait>
   Problem(TrialFunction<TrialFEC, Trait>&, TestFunction<TestFEC, Trait>&)
      -> Problem<TrialFEC, TestFEC, Trait>;
}

#include "Problem.hpp"

#endif
