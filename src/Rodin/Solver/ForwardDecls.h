/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_FORWARDDECLS_H
#define RODIN_SOLVER_FORWARDDECLS_H

namespace Rodin::Solver
{
  template <class OperatorType, class VectorType>
  class SolverBase;

  /**
   * @brief Conjugate gradient solver for self-adjoint problems.
   * @see CGSpecializations
   */
  template <class OperatorType, class VectorType>
  class CG;

  /**
   * @brief Direct sparse LLT Cholesky factorizations
   * @see SimplicialLLTSpecializations
   */
  template <class OperatorType, class VectorType>
  class SimplicialLLT;

  /**
   * @brief Direct sparse LDLT Cholesky factorizations without square root
   * @see SimplicialLDLTSpecializations
   */
  template <class OperatorType, class VectorType>
  class SimplicialLDLT;

  /**
   * @brief Sparse supernodal LU factorization for general matrices.
   * @see SparseLUSpecializations
   */
  template <class OperatorType, class VectorType>
  class SparseLU;

  /**
   * @brief Sparse left-looking QR factorization with numerical column
   * pivoting.
   * @see SparseQRSpecializations
   */
  template <class OperatorType, class VectorType>
  class SparseQR;

  template <class OperatorType, class VectorType>
  class LeastSquaresCG;

  template <class OperatorType, class VectorType>
  class BiCGSTAB;

  template <class OperatorType, class VectorType>
  class UMFPack;
}

#endif
