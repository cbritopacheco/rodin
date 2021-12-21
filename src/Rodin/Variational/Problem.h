/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_PROBLEM_H
#define RODIN_VARIATIONAL_PROBLEM_H

#include <functional>

#include <mfem.hpp>

#include "Rodin/Alert.h"
#include "Rodin/Mesh.h"

#include "FormLanguage/ForwardDecls.h"

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

         virtual mfem::Array<int>& getEssentialBoundary() = 0;

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

         virtual void assemble()
         {
            getLinearForm().assemble();
            getBilinearForm().assemble();
         }
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
    * @see FormLanguage::ProblemBody
    */
   template <class FEC>
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
         Problem(GridFunction<FEC>& u);

         Problem& operator=(const FormLanguage::ProblemBody& rhs) override;

         GridFunction<FEC>& getSolution() override;
         LinearForm<FEC>& getLinearForm() override;
         BilinearForm<FEC>& getBilinearForm() override;
         mfem::Array<int>& getEssentialBoundary() override;

      private:
         std::reference_wrapper<GridFunction<FEC>>          m_solution;

         std::vector<mfem::Array<int>> m_bdr;

         BilinearForm<FEC>    m_bilinearForm;
         LinearForm<FEC>      m_linearForm;

         mfem::Array<int> m_essBdr;

         std::unique_ptr<FormLanguage::ProblemBody> m_pb;
   };
}

#include "Problem.hpp"

#endif
