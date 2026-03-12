/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_FORWARDDECLS_H
#define RODIN_SOLVER_FORWARDDECLS_H

/**
 * @file ForwardDecls.h
 * @brief Forward declarations for all solver classes
 *
 * This file contains forward declarations and brief descriptions for all
 * solver classes in the Rodin::Solver module.
 */

#include "Rodin/Configure.h"

namespace Rodin::Solver
{
  /**
   * @brief Abstract base class for linear system solvers.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Represents an object which can solve systems of the type:
   * @f[
   *  Ax = b \: ,
   * @f]
   * where @f$ A @f$ is the operator (matrix), @f$ x @f$ is the solution vector,
   * and @f$ b @f$ is the right-hand side vector.
   *
   * @see LinearSolverBase for the full implementation.
   */
  template <class LinearSystem>
  class LinearSolverBase;

  /**
   * @brief Wrapper class for any Eigen sparse solver.
   * @tparam EigenSolverType Type of the underlying Eigen solver
   * @tparam OperatorType Type of operator (matrix)
   * @tparam VectorType Type of vector
   *
   * This class provides a uniform interface to various Eigen solver types.
   */
  template <class EigenSolverType, class OperatorType, class VectorType>
  class EigenSolver;

  /**
   * @brief Conjugate gradient solver for symmetric positive definite systems.
   * @tparam LinearSystem Type of linear system to solve
   *
   * The CG method is an iterative solver particularly efficient for large
   * sparse symmetric positive definite systems.
   *
   * @see CGSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class CG;

  /**
   * @brief Direct sparse LLT Cholesky factorization solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs Cholesky decomposition @f$ A = LL^T @f$ for symmetric positive
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
   * Performs Cholesky decomposition @f$ A = LDL^T @f$ without square root
   * for symmetric positive definite sparse matrices.
   *
   * @see SimplicialLDLTSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SimplicialLDLT;

  /**
   * @brief Robust LDLT Cholesky decomposition for dense matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs robust Cholesky decomposition with pivoting for dense matrices.
   *
   * @see LDLTSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class LDLT;

  /**
   * @brief Householder QR decomposition for dense matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs QR decomposition using Householder reflections for dense matrices.
   *
   * @see HouseholderQRSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class HouseholderQR;

  /**
   * @brief Sparse supernodal LU factorization for general matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs LU decomposition for general (non-symmetric) sparse matrices.
   *
   * @see SparseLUSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SparseLU;

  /**
   * @brief Sparse QR factorization with column pivoting.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs QR decomposition with numerical column pivoting for sparse matrices.
   *
   * @see SparseQRSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class SparseQR;

  /**
   * @brief Least-squares conjugate gradient solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Iterative solver for least-squares problems @f$ \min_x \|Ax - b\|^2 @f$.
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
   */
  template <class LinearSystem>
  class GMRES;

  /**
   * @brief Deflated generalized minimal residual iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * DGMRES is a variant of GMRES with deflation techniques.
   */
  template <class LinearSystem>
  class DGMRES;

  template <class LinearSystem>
  class IDRS;

  /**
   * @brief IDR(s)STABL iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * IDR(s)STABL is an induced dimension reduction method with stabilization.
   */
  template <class LinearSystem>
  class IDRSTABL;

#ifdef RODIN_USE_SPQR
  /**
   * @brief SuiteSparseQR multifrontal sparse QR factorization.
   * @tparam LinearSystem Type of linear system to solve
   *
   * High-performance sparse QR factorization from the SuiteSparse library.
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
   * High-performance sparse direct solver from the SuiteSparse library.
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
   * This namespace contains wrappers for CHOLMOD solvers.
   */
  namespace CHOLMOD
  {
    /**
     * @brief CHOLMOD supernodal LLT Cholesky factorization.
     * @tparam LinearSystem Type of linear system to solve
     *
     * Supernodal LLT factorization using CHOLMOD from SuiteSparse.
     */
    template <class LinearSystem>
    class SupernodalLLT;
  }
#endif
}

#endif
