/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SparseProblem.h
 * @brief Sparse matrix problem specialization.
 *
 * This file defines the SparseProblem class, a specialization of the Problem
 * class that uses sparse matrix storage for efficient handling of large-scale
 * finite element systems. Sparse storage is the standard choice for most FEM
 * problems as system matrices are typically very sparse.
 *
 * ## Sparse Matrix Storage
 * For a typical FEM problem with @f$ N @f$ degrees of freedom:
 * - **Dense storage**: @f$ O(N^2) @f$ memory
 * - **Sparse storage**: @f$ O(N) @f$ memory (for 2D/3D problems)
 *
 * ## Sparsity Pattern
 * The sparsity pattern is determined by the mesh connectivity:
 * @f[
 *   A_{ij} \neq 0 \iff \text{DOFs } i \text{ and } j \text{ share an element}
 * @f]
 *
 * ## Advantages
 * - Memory efficient for large problems
 * - Faster matrix-vector products
 * - Enables iterative solvers (CG, GMRES, etc.)
 * - Essential for 3D problems
 *
 * ## Usage Example
 * ```cpp
 * P1 Vh(mesh);
 * TrialFunction u(Vh);
 * TestFunction v(Vh);
 * 
 * // Explicitly uses sparse storage
 * SparseProblem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * problem.solve(solver);
 * ```
 */
#ifndef RODIN_VARIATIONAL_SPARSEPROBLEM_H
#define RODIN_VARIATIONAL_SPARSEPROBLEM_H

#include "Problem.h"

namespace Rodin::Variational
{
  template <class ... Parameters>
  class SparseProblem;

  /**
   * @defgroup SparseProblemSpecializations SparseProblem Template Specializations
   * @brief Template specializations of the SparseProblem class.
   * @see SparseProblem
   */

  /**
   * @ingroup SparseProblemSpecializations
   * @brief General class to assemble linear systems with `Math::SparseMatrix`
   * and `Math::Vector` types in a serial context.
   */
  template <class LinearSystem, class U, class V>
  class SparseProblem<LinearSystem, U, V> : public Problem<LinearSystem, U, V>
  {
    public:
      using Parent = Problem<LinearSystem, U, V>;
      using Parent::Parent;
      using Parent::operator=;
  };

  /**
   * @ingroup RodinCTAD
   */
  template <class U, class V>
  SparseProblem(U& u, V& v)
    -> SparseProblem<
        Math::LinearSystem<
          Math::SparseMatrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<typename FormLanguage::Traits<U>::FESType>::ScalarType,
              typename FormLanguage::Traits<typename FormLanguage::Traits<V>::FESType>::ScalarType>::Type>,
          Math::Vector<
            typename FormLanguage::Traits<typename FormLanguage::Traits<V>::FESType>::ScalarType>>,
          U, V>;
}

#endif

