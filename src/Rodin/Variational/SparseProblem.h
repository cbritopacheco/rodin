/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SPARSEPROBLEM_H
#define RODIN_VARIATIONAL_SPARSEPROBLEM_H

/**
 * @file
 * @brief Sparse linear system problem for finite element methods.
 *
 * This file defines the SparseProblem class which assembles finite element
 * problems into sparse matrix representations using `Math::SparseMatrix` and
 * `Math::Vector`. Sparse problems are the standard choice for FEM because:
 * - Finite element matrices are naturally sparse
 * - Memory efficiency for large-scale problems
 * - Fast assembly and solve times
 *
 * ## Sparsity Pattern
 * The sparsity of the system matrix @f$ A @f$ reflects the mesh connectivity:
 * - @f$ A_{ij} \neq 0 @f$ only if basis functions @f$ \phi_i @f$ and @f$ \phi_j @f$ overlap
 * - For FEM, typically @f$ O(N) @f$ non-zeros in @f$ N \times N @f$ matrix
 *
 * ## When to Use
 * SparseProblem is the default choice for:
 * - Standard finite element problems
 * - Large-scale simulations
 * - Most practical applications
 *
 * @see DenseProblem, Problem
 */

#include "Problem.h"

namespace Rodin::Variational
{
  template <class ... Parameters>
  class SparseProblem;

  template <class TrialFES, class TestFES>
  SparseProblem(TrialFunction<TrialFES>&, TestFunction<TestFES>&)
    -> SparseProblem<TrialFES, TestFES,
          Math::Matrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
          Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>;

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
  template <class TrialFES, class TestFES>
  class SparseProblem<
    TrialFES, TestFES,
    Math::Matrix<
      typename FormLanguage::Mult<
        typename FormLanguage::Traits<TrialFES>::ScalarType,
        typename FormLanguage::Traits<TestFES>::ScalarType>::Type>,
    Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    : public Problem<
        TrialFES, TestFES,
        Math::SparseMatrix<
          typename FormLanguage::Mult<
            typename FormLanguage::Traits<TrialFES>::ScalarType,
            typename FormLanguage::Traits<TrialFES>::ScalarType>
          ::Type>,
        Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>
    {
      public:
        using Parent = Problem<
          TrialFES, TestFES,
          Math::SparseMatrix<
            typename FormLanguage::Mult<
              typename FormLanguage::Traits<TrialFES>::ScalarType,
              typename FormLanguage::Traits<TrialFES>::ScalarType>
            ::Type>,
          Math::Vector<typename FormLanguage::Traits<TestFES>::ScalarType>>;

        using Parent::Parent;
        using Parent::operator=;
    };
}

#endif

