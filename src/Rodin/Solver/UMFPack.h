/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_UMFPACK_H
#define RODIN_SOLVER_UMFPACK_H

#include <optional>
#include <functional>

#include <mfem.hpp>

#include "Rodin/Variational/Problem.h"

#include "Solver.h"

namespace Rodin::Solver
{
   /**
    * @brief UMFPack
    */
   class UMFPack : public Solver
   {
      public:
         /**
          * @brief Constructs the UMFPack object with default parameters.
          */
         UMFPack(bool useLongInts = false)
            : m_umfpack(useLongInts)
         {}

         ~UMFPack() = default;

         UMFPack& useLongInts(bool v = true);

         void solve(Variational::ProblemBase& problem) override;

      private:
         bool m_useLongInts;
         mfem::UMFPackSolver m_umfpack;
   };
}

#endif


