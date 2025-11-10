/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_FORWARDDECLS_H
#define RODIN_SOLVER_FORWARDDECLS_H

#include "Rodin/Configure.h"

namespace Rodin::Solver
{
  /**
   * @brief Abstract base class for linear algebra solvers.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Represents an object which can solve systems of the type:
   * @f[
   *  Ax = b \: ,
   * @f]
   * where @f$ A @f$ is the operator (matrix), the solution @f$ x @f$
   * is the unknown vector, and the right hand side @f$ b @f$ is a known vector.
   *
   * @see SolverBase for the full implementation details.
   */
  template <class LinearSystem>
  class SolverBase;

  /**
   * @brief Wrapper class for any Eigen sparse solver.
   * @tparam EigenSolverType Type of the underlying Eigen solver
   * @tparam OperatorType Type of operator (matrix) for the left hand side
   * @tparam VectorType Type of vector for the right hand side and solution
   *
   * This class provides a uniform interface to various Eigen solver types,
   * allowing integration with Rodin's solver framework.
   */
  template <class EigenSolverType, class OperatorType, class VectorType>
  class EigenSolver;

  /**
   * @brief Conjugate gradient solver for self-adjoint problems.
   * @tparam LinearSystem Type of linear system to solve
   *
   * The conjugate gradient (CG) method is an iterative algorithm for solving
   * linear systems @f$ Ax = b @f$ where @f$ A @f$ is a symmetric positive
   * definite matrix. It is particularly efficient for large sparse systems.
   *
   * @see CGSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class CG;

  /**
   * @brief Direct sparse LLT Cholesky factorization solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs Cholesky decomposition @f$ A = LL^T @f$ where @f$ L @f$ is a
   * lower triangular matrix. This solver is suitable for symmetric positive
   * definite sparse matrices.
   *
   * @see SimplicialLLTSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SimplicialLLT;

  /**
   * @brief Direct sparse LDLT Cholesky factorization solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs Cholesky decomposition @f$ A = LDL^T @f$ where @f$ L @f$ is a
   * unit lower triangular matrix and @f$ D @f$ is a diagonal matrix. This
   * variant avoids square root operations and is suitable for symmetric
   * positive definite sparse matrices.
   *
   * @see SimplicialLDLTSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SimplicialLDLT;

  /**
   * @brief Robust LDLT Cholesky decomposition for dense matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs robust Cholesky decomposition @f$ A = LDL^T @f$ with pivoting
   * for dense matrices. This solver handles symmetric indefinite matrices
   * and provides numerical stability through pivoting.
   *
   * @see LDLTSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class LDLT;

  /**
   * @brief Householder QR decomposition for dense matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs QR decomposition @f$ A = QR @f$ using Householder reflections,
   * where @f$ Q @f$ is orthogonal and @f$ R @f$ is upper triangular. This
   * solver is suitable for dense matrices and provides good numerical stability.
   *
   * @see HouseholderQRSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class HouseholderQR;

  /**
   * @brief Sparse supernodal LU factorization for general matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs LU decomposition @f$ A = LU @f$ where @f$ L @f$ is lower
   * triangular and @f$ U @f$ is upper triangular. This direct solver is
   * suitable for general (non-symmetric) sparse matrices.
   *
   * @see SparseLUSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SparseLU;

  /**
   * @brief Sparse QR factorization with column pivoting.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs QR decomposition @f$ AP = QR @f$ using a left-looking algorithm
   * with numerical column pivoting, where @f$ P @f$ is a permutation matrix,
   * @f$ Q @f$ is orthogonal, and @f$ R @f$ is upper triangular. This solver
   * is suitable for rank-deficient sparse matrices.
   *
   * @see SparseQRSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SparseQR;

  /**
   * @brief Least-squares conjugate gradient solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Iterative solver for solving least-squares problems @f$ \min_x \|Ax - b\|^2 @f$
   * using the conjugate gradient method applied to the normal equations
   * @f$ A^TAx = A^Tb @f$. Suitable for overdetermined or rectangular systems.
   *
   * @see LeastSquaresCGSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class LeastSquaresCG;

  /**
   * @brief Bi-conjugate gradient stabilized iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * The BiCGSTAB method is an iterative solver for non-symmetric linear systems.
   * It combines ideas from the BiCG algorithm with stabilization to improve
   * convergence and avoid irregular convergence patterns.
   *
   * @see BiCGSTABSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class BiCGSTAB;

  /**
   * @brief Generalized minimal residual iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * GMRES is an iterative method for solving non-symmetric linear systems.
   * It minimizes the residual norm over Krylov subspaces.
   */
  template <class LinearSystem>
  class GMRES;

  /**
   * @brief Deflated generalized minimal residual iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * DGMRES is a variant of GMRES with deflation techniques to improve
   * convergence for difficult problems.
   */
  template <class LinearSystem>
  class DGMRES;

  /**
   * @brief IDR(s)STABL iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * IDR(s)STABL is an induced dimension reduction method with stabilization
   * for solving non-symmetric linear systems.
   */
  template <class LinearSystem>
  class IDRSTABL;

#ifdef RODIN_USE_SPQR
  /**
   * @brief SuiteSparseQR multifrontal sparse QR factorization.
   * @tparam LinearSystem Type of linear system to solve
   *
   * SPQR is a high-performance sparse QR factorization from the SuiteSparse
   * library. It uses multifrontal techniques and is particularly effective
   * for large sparse least-squares problems.
   *
   * @see SPQRSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SPQR;
#endif

#ifdef RODIN_USE_UMFPACK
  /**
   * @brief UMFPACK multifrontal sparse LU factorization.
   * @tparam LinearSystem Type of linear system to solve
   *
   * UMFPack is a high-performance sparse direct solver from the SuiteSparse
   * library that uses unsymmetric multifrontal LU factorization. It is
   * particularly effective for large sparse unsymmetric systems.
   *
   * @see UMFPackSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class UMFPack;
#endif

#ifdef RODIN_USE_CHOLMOD
  /**
   * @brief CHOLMOD solver implementations from SuiteSparse.
   *
   * This namespace contains wrappers for CHOLMOD solvers, which provide
   * high-performance sparse Cholesky factorizations from the SuiteSparse library.
   */
  namespace CHOLMOD
  {
    /**
     * @brief CHOLMOD supernodal LLT Cholesky factorization.
     * @tparam LinearSystem Type of linear system to solve
     *
     * Supernodal LLT factorization for symmetric positive definite sparse
     * matrices using CHOLMOD from SuiteSparse. This solver exploits supernodal
     * structure for improved performance on dense blocks within the sparse matrix.
     */
    template <class LinearSystem>
    class SupernodalLLT;
  }
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
