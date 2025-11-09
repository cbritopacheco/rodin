/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DenseProblem.h
 * @brief Dense matrix problem specialization.
 *
 * This file defines the DenseProblem class, a specialization of the Problem
 * class that uses dense matrix storage. While less memory-efficient than
 * sparse storage, dense matrices can be advantageous for small problems or
 * when matrix structure doesn't favor sparsity.
 *
 * ## Dense Matrix Storage
 * Dense storage stores all @f$ N^2 @f$ matrix entries regardless of zeros:
 * - **Memory**: @f$ O(N^2) @f$
 * - **Access**: @f$ O(1) @f$ for any entry
 * - **Operations**: Cache-friendly for small matrices
 *
 * ## When to Use Dense Storage
 * - Small problems (few hundred DOFs)
 * - Nearly full matrices (rare in standard FEM)
 * - Prototyping and testing
 * - Spectral methods with global basis functions
 * - Problems with dense blocks (e.g., boundary elements)
 *
 * ## Performance Considerations
 * For typical FEM problems:
 * - Sparse: scales to millions of DOFs
 * - Dense: practical only for ~1000 DOFs or less
 *
 * ## Usage Example
 * ```cpp
 * P1 Vh(mesh);  // Small mesh
 * TrialFunction u(Vh);
 * TestFunction v(Vh);
 * 
 * // Use dense storage
 * DenseProblem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * problem.solve(solver);
 * ```
 */
#ifndef RODIN_VARIATIONAL_DENSEPROBLEM_H
#define RODIN_VARIATIONAL_DENSEPROBLEM_H

#include "Problem.h"

namespace Rodin::Variational
{
  template <class ... Parameters>
  class DenseProblem;

  /**
   * @defgroup DenseProblemSpecializations DenseProblem Template Specializations
   * @brief Template specializations of the DenseProblem class.
   * @see DenseProblem
   */

  /**
   * @ingroup DenseProblemSpecializations
   * @brief General class to assemble linear systems with `Math::Matrix`
   * and `Math::Vector` types in a serial context.
   */
  template <class LinearSystem, class U, class V>
  class DenseProblem<LinearSystem, U, V> : public Problem<LinearSystem, U, V>
  {
    public:
      using Parent = Problem<LinearSystem, U, V>;
      using Parent::Parent;
      using Parent::operator=;
  };

  template <class U, class V>
  DenseProblem(U& u, V& v)
    -> DenseProblem<
        Math::LinearSystem<
          Math::Matrix<
            typename FormLanguage::Traits<typename FormLanguage::Traits<U>::FESType>
            ::ScalarType>,
          Math::Vector<
            typename FormLanguage::Traits<typename FormLanguage::Traits<V>::FESType>
            ::ScalarType>>,
          U, V>;
}

#endif

