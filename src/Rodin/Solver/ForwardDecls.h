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
   * @brief Abstract base class for Newton-type nonlinear solvers.
   * @tparam LinearSolver Type of the linear solver.
   */
  template <class LinearSolver>
  class NewtonSolverBase;

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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref CG "CG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Iterative conjugate gradient solver for sparse SPD systems. |
   * | @ref CG "CG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Iterative conjugate gradient solver for dense SPD systems. |
   * | @ref CG "CG<PETSc::Math::LinearSystem>" | PETSc KSP-backed conjugate gradient solver. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref SimplicialLLT "SimplicialLLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Eigen sparse LLT factorization for SPD sparse systems. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref SimplicialLDLT "SimplicialLDLT<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Eigen sparse LDLT factorization for SPD sparse systems. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref LDLT "LDLT<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Eigen dense LDLT factorization with pivoting. |
   *
   * @see LDLTSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class LDLT;


  /**
   * @brief Dense LU factorization with partial pivoting for general matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs LU decomposition with partial pivoting for dense matrices.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref PartialPivLU "PartialPivLU<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Eigen dense LU factorization with partial pivoting. |
   *
   * @see PartialPivLUSpecializations for available template specializations.
   */
  template <class LinearSystem>
  class PartialPivLU;

  /**
   * @brief Householder QR decomposition for dense matrices.
   * @tparam LinearSystem Type of linear system to solve
   *
   * Performs QR decomposition using Householder reflections for dense matrices.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref HouseholderQR "HouseholderQR<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Eigen dense Householder QR factorization. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref SparseLU "SparseLU<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Eigen sparse LU factorization for general sparse systems. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref SparseQR "SparseQR<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Eigen sparse QR factorization with column pivoting. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref LeastSquaresCG "LeastSquaresCG<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Least-squares conjugate gradient solver for sparse systems. |
   * | @ref LeastSquaresCG "LeastSquaresCG<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Least-squares conjugate gradient solver for dense systems. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref BiCGSTAB "BiCGSTAB<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Stabilized bi-conjugate gradient solver for sparse non-symmetric systems. |
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
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref GMRES "GMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Restarted GMRES solver for sparse non-symmetric systems. |
   * | @ref GMRES "GMRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Restarted GMRES solver for dense non-symmetric systems. |
   * | @ref GMRES "GMRES<PETSc::Math::LinearSystem>" | PETSc KSP-backed GMRES solver. |
   */
  template <class LinearSystem>
  class GMRES;

  /**
   * @brief Deflated generalized minimal residual iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * DGMRES is a variant of GMRES with deflation techniques.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref DGMRES "DGMRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Deflated GMRES solver for sparse non-symmetric systems. |
   * | @ref DGMRES "DGMRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Deflated GMRES solver for dense non-symmetric systems. |
   */
  template <class LinearSystem>
  class DGMRES;

  /**
   * @brief Induced dimension reduction iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref IDRS "IDRS<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | IDR(s) solver for sparse non-symmetric systems. |
   * | @ref IDRS "IDRS<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | IDR(s) solver for dense non-symmetric systems. |
   */
  template <class LinearSystem>
  class IDRS;

  /**
   * @brief Minimum residual iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref MINRES "MINRES<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | MINRES solver for sparse symmetric indefinite systems. |
   * | @ref MINRES "MINRES<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | MINRES solver for dense symmetric indefinite systems. |
   */
  template <class LinearSystem>
  class MINRES;

  /**
   * @brief IDR(s)STABL iterative solver.
   * @tparam LinearSystem Type of linear system to solve
   *
   * IDR(s)STABL is an induced dimension reduction method with stabilization.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref IDRSTABL "IDRSTABL<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | Stabilized IDR(s) solver for sparse non-symmetric systems. |
   * | @ref IDRSTABL "IDRSTABL<Math::LinearSystem<Math::Matrix<Scalar>, Math::Vector<Scalar>>>" | Stabilized IDR(s) solver for dense non-symmetric systems. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref SPQR "SPQR<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | SuiteSparseQR factorization for sparse systems. |
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
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref UMFPack "UMFPack<Math::LinearSystem<Math::SparseMatrix<Scalar>, Math::Vector<Scalar>>>" | UMFPACK multifrontal sparse LU factorization. |
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
