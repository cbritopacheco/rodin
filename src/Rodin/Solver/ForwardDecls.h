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
  /**
   * @brief Abstract base class for linear algebra solvers.
   * @tparam OperatorType Type of operator for the left hand side
   * @tparam VectorType Type of vector for the right hand side
   *
   * Represents an object which can solve sysystems of the type:
   * @f[
   *  Ax = b \: ,
   * @f]
   * where @f$ A @f$ has type @f$ \text{OperatorType} @f$, the solution @f$ x @f$
   * has type @f$ \text{VectorType} @f$, and the right hand side @f$ b @f$ has
   * type @f$ \text{VectorType} @f$.
   */
  template <class OperatorType, class VectorType>
  class SolverBase;

  /**
   * @brief Wrapper class for any Eigen sparse solver.
   * @tparam OperatorType Type of operator for the left hand side
   * @tparam VectorType Type of vector for the right hand side
   * @see EigenSolverSpecializations
   */
  template <class EigenSolverType, class OperatorType, class VectorType>
  class EigenSolver;

  /**
   * @brief Conjugate gradient solver for self-adjoint problems.
   * @tparam OperatorType Type of operator for the left hand side
   * @tparam VectorType Type of vector for the right hand side
   * @see CGSpecializations
   */
  template <class OperatorType, class VectorType>
  class CG;

  /**
   * @brief Direct sparse LLT Cholesky factorizations
   * @tparam OperatorType Type of operator for the left hand side
   * @tparam VectorType Type of vector for the right hand side
   * @see SimplicialLLTSpecializations
   */
  template <class OperatorType, class VectorType>
  class SimplicialLLT;

  /**
   * @brief Direct sparse LDLT Cholesky factorizations without square root
   * @tparam OperatorType Type of operator for the left hand side
   * @tparam VectorType Type of vector for the right hand side
   * @see SimplicialLDLTSpecializations
   */
  template <class OperatorType, class VectorType>
  class SimplicialLDLT;

  /**
   * @brief Robust Cholesky decomposition of a dense matrix with pivoting.
   */
  template <class OperatorType, class VectorType>
  class LDLT;

  /**
   * @brief Robust Cholesky decomposition of a dense matrix with pivoting.
   */
  template <class OperatorType, class VectorType>
  class HouseholderQR;

  /**
   * @brief Sparse supernodal LU factorization for general matrices.
   * @tparam OperatorType Type of operator for the left hand side
   * @tparam VectorType Type of vector for the right hand side
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

#ifdef RODIN_USE_UMFPACK

  template <class OperatorType, class VectorType>
  class UMFPack;

#endif

#ifdef RODIN_USE_CHOLMOD

  template <class OperatorType, class VectorType>
  class Cholmod;

#endif

#ifdef RODIN_USE_PASTIX
#endif

#ifdef RODIN_USE_KLU
#endif

#ifdef RODIN_USE_KLU
#endif

#ifdef RODIN_USE_SUPERLU
#endif

#ifdef RODIN_USE_SPQR
#endif

#ifdef RODIN_USE_PARDISO
#endif

#ifdef RODIN_USE_APPLE_ACCELERATE
#endif
}

#endif
