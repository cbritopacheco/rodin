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

#include "Rodin/Variational/Problem.h"

#include "Solver.h"

namespace Rodin::Solver
{
   /**
    * @brief Preconditioned Conjugate Gradient
    */
   class CG : public Solver
   {
      public:
         /**
          * @brief Constructs the CG object with default parameters.
          */
         CG() = default;

         ~CG() = default;

         /**
          * @brief Sets whether some information will be printed at each
          * iteration.
          * @param[in] printIterations If set to true, will print the
          * iterations. Otherwise, no information will be printed to the
          * screen.
          * @returns Reference to self (for method chaining)
          */
         CG& printIterations(bool printIterations);

         /**
          * @brief Sets the maximum amount of iterations the solver will
          * perform.
          * @param[in] maxIterations Maximum amount of iterations
          * @returns Reference to self (for method chaining)
          */
         CG& setMaxIterations(int maxIterations);

         /**
          * @brief Sets the relative tolerance of the solver.
          * @param[in] rtol Relative tolerance
          * @returns Reference to self (for method chaining)
          */
         CG& setRelativeTolerance(double rtol);

         /**
          * @brief Sets the absolute tolerance of the solver.
          * @param[in] atol Absolute tolerance
          * @returns Reference to self (for method chaining)
          */
         CG& setAbsoluteTolerance(double atol);

         void solve(Variational::ProblemBase& problem) const override;

      private:
         bool m_printIterations;
         int  m_maxIterations;
         double m_rtol, m_atol;
   };
}

#endif

