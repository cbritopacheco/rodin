/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_SOLVER_H
#define RODIN_SOLVER_SOLVER_H

#include "Rodin/Variational/Problem.h"

namespace Rodin::Solver
{
   /**
    * @brief Base class for solving variational problems represented by a
    * Variational::Problem instance.
    */
   class Solver
   {
      public:
         /**
          * @brief Default virtual destructor.
          */
         virtual ~Solver() = default;

         /**
          * Solves the specified Variational::Problem.
          *
          * @param[in,out] problem Variational problem to solve.
          */
         virtual void solve(Variational::ProblemBase& problem) const = 0;
   };
}

#endif
