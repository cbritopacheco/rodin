/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_UMFPACK_H
#define RODIN_SOLVER_UMFPACK_H

#ifdef RODIN_USE_SUITESPARSE

#include <optional>
#include <functional>

#include "Rodin/Configure.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/SparseMatrix.h"

#include "ForwardDecls.h"
#include "Solver.h"

namespace Rodin::Solver
{
  UMFPack() -> UMFPack<Math::SparseMatrix, Math::Vector>;

  /**
   * @defgroup UMFPackSpecializations UMFPack Template Specializations
   * @brief Template specializations of the UMFPack class.
   * @see UMFPack
   */

  /**
   * @ingroup UMFPackSpecializations
   * @brief UMFPack for use with `Math::SparseMatrix` and `Math::Vector`.
   */
  template <>
  class UMFPack<Math::SparseMatrix, Math::Vector>
    : public SolverBase<Math::SparseMatrix, Math::Vector>
  {
    public:
      using OperatorType = Math::SparseMatrix;
      using VectorType = Math::Vector;

      /**
       * @brief Constructs the UMFPack object with default parameters.
       */
      UMFPack() = default;

      ~UMFPack() = default;

      void solve(OperatorType& A, VectorType& x, VectorType& b) override
      {
        assert(false);
      }
  };
}

#endif // #ifdef RODIN_USE_SUITESPARSE
#endif


