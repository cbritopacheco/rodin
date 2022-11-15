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

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
   /**
    * @brief UMFPack
    */
   template <>
   class UMFPack<mfem::SparseMatrix, mfem::Vector>
      : public SolverBase<mfem::SparseMatrix, mfem::Vector>
   {
      public:
         using OperatorType = mfem::SparseMatrix;
         using VectorType = mfem::Vector;

         /**
          * @brief Constructs the UMFPack object with default parameters.
          */
         UMFPack()
            : m_useLongInts(false)
         {}

         ~UMFPack() = default;

         UMFPack& useLongInts(bool v = true)
         {
            m_useLongInts = v;
            return *this;
         }

         void solve(OperatorType& stiffness, VectorType& mass, VectorType& solution) const override
         {
            mfem::UMFPackSolver umfpack(m_useLongInts);
            umfpack.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
            umfpack.SetOperator(stiffness);
            umfpack.Mult(mass, solution);
         }

      private:
         bool m_useLongInts;
   };
}

#endif


