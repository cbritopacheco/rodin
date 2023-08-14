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
  UMFPack() -> UMFPack<mfem::SparseMatrix, mfem::Vector>;

  /**
   * @defgroup UMFPackSpecializations UMFPack Template Specializations
   * @brief Template specializations of the UMFPack class.
   * @see UMFPack
   */

  /**
   * @ingroup UMFPackSpecializations
   * @brief UMFPack for use with `mfem::SparseMatrix` and `mfem::Vector`.
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

      void solve(OperatorType& A, VectorType& x, VectorType& b) const override
      {
        mfem::UMFPackSolver umfpack(m_useLongInts);
        umfpack.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
        umfpack.SetOperator(A);
        umfpack.Mult(b, x);
      }

    private:
      bool m_useLongInts;
  };
}

#endif


