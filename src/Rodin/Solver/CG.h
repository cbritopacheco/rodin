/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_CG_H
#define RODIN_SOLVER_CG_H

#include <optional>
#include <functional>

#include <mfem.hpp>

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
   CG() -> CG<mfem::SparseMatrix, mfem::Vector>;

   /**
    * @defgroup CGSpecializations CG Template Specializations
    * @brief Template specializations of the CG class.
    * @see CG
    */

   /**
    * @ingroup CGSpecializations
    * @brief Conjugate Gradient for use with `mfem::Operator` and
    * `mfem::Vector`.
    */
   template <>
   class CG<mfem::Operator, mfem::Vector> : public SolverBase<mfem::Operator, mfem::Vector>
   {
      public:
         using OperatorType = mfem::Operator;
         using VectorType = mfem::Vector;

         /**
          * @brief Constructs the CG object with default parameters.
          */
         CG()
            : m_maxIterations(200),
              m_printIterations(false),
              m_rtol(1e-12),
              m_atol(0)
         {}

         ~CG() = default;

         /**
          * @brief Sets whether some information will be printed at each
          * iteration.
          * @param[in] printIterations If set to true, will print the
          * iterations. Otherwise, no information will be printed to the
          * screen.
          * @returns Reference to self (for method chaining)
          */
         CG& printIterations(bool printIterations)
         {
            m_printIterations = printIterations;
            return *this;
         }

         /**
          * @brief Sets the maximum amount of iterations the solver will
          * perform.
          * @param[in] maxIterations Maximum amount of iterations
          * @returns Reference to self (for method chaining)
          */
         CG& setMaxIterations(int maxIterations)
         {
            m_maxIterations = maxIterations;
            return *this;
         }

         /**
          * @brief Sets the relative tolerance of the solver.
          * @param[in] rtol Relative tolerance
          * @returns Reference to self (for method chaining)
          */
         CG& setRelativeTolerance(double rtol)
         {
            m_rtol = rtol;
            return *this;
         }

         /**
          * @brief Sets the absolute tolerance of the solver.
          * @param[in] atol Absolute tolerance
          * @returns Reference to self (for method chaining)
          */
         CG& setAbsoluteTolerance(double atol)
         {
            m_atol = atol;
            return *this;
         }

         virtual
         void solve(OperatorType& A, VectorType& X, VectorType& B)
         const override
         {
            mfem::CGSolver pcg;
            pcg.SetPrintLevel(static_cast<int>(m_printIterations));
            pcg.SetMaxIter(m_maxIterations);
            pcg.SetRelTol(sqrt(m_rtol));
            pcg.SetAbsTol(sqrt(m_atol));
            pcg.SetOperator(A);
            pcg.Mult(B, X);
         }

      private:
         int  m_maxIterations;
         bool m_printIterations;
         double m_rtol, m_atol;
   };

   /**
    * @ingroup CGSpecializations
    * @brief Conjugate Gradient for use with `mfem::SparseMatrix` and
    * `mfem::Vector`.
    */
   template <>
   class CG<mfem::SparseMatrix, mfem::Vector>
      : public SolverBase<mfem::SparseMatrix, mfem::Vector>
   {
      public:
         using OperatorType = mfem::SparseMatrix;
         using VectorType = mfem::Vector;

         /**
          * @brief Constructs the CG object with default parameters.
          */
         CG()
            : m_maxIterations(200),
              m_printIterations(false),
              m_rtol(1e-12),
              m_atol(0)
         {}

         ~CG() = default;

         /**
          * @brief Sets whether some information will be printed at each
          * iteration.
          * @param[in] printIterations If set to true, will print the
          * iterations. Otherwise, no information will be printed to the
          * screen.
          * @returns Reference to self (for method chaining)
          */
         CG& printIterations(bool printIterations)
         {
            m_printIterations = printIterations;
            return *this;
         }

         /**
          * @brief Sets the maximum amount of iterations the solver will
          * perform.
          * @param[in] maxIterations Maximum amount of iterations
          * @returns Reference to self (for method chaining)
          */
         CG& setMaxIterations(int maxIterations)
         {
            m_maxIterations = maxIterations;
            return *this;
         }

         /**
          * @brief Sets the relative tolerance of the solver.
          * @param[in] rtol Relative tolerance
          * @returns Reference to self (for method chaining)
          */
         CG& setRelativeTolerance(double rtol)
         {
            m_rtol = rtol;
            return *this;
         }

         /**
          * @brief Sets the absolute tolerance of the solver.
          * @param[in] atol Absolute tolerance
          * @returns Reference to self (for method chaining)
          */
         CG& setAbsoluteTolerance(double atol)
         {
            m_atol = atol;
            return *this;
         }

         CG& setPreconditioner(mfem::Solver& smoother)
         {
            m_smoother.emplace(std::ref(smoother));
            return *this;
         }

         virtual
         void solve(OperatorType& A, VectorType& x, VectorType& b)
         const override
         {
            mfem::CGSolver pcg;
            pcg.SetPrintLevel(static_cast<int>(m_printIterations));
            pcg.SetMaxIter(m_maxIterations);
            pcg.SetRelTol(sqrt(m_rtol));
            pcg.SetAbsTol(sqrt(m_atol));
            if (m_smoother)
               pcg.SetPreconditioner(*m_smoother);
            pcg.SetOperator(A);
            pcg.Mult(b, x);
         }

      private:
         int  m_maxIterations;
         bool m_printIterations;
         double m_rtol, m_atol;
         std::optional<std::reference_wrapper<mfem::Solver>> m_smoother;
   };

}

#endif

