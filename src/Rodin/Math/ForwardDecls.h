/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_FORWARDDECLS_H
#define RODIN_MATH_FORWARDDECLS_H

#include "Matrix.h"
#include "Vector.h"

#include "SparseMatrix.h"

namespace Rodin::Math
{
  /**
   * @brief Represents a linear system of equations.
   * @tparam Operator Type of the operator
   * @tparam Vector Type of the vector
   *
   * Represents the mathematical expression:
   * @f[
   *   \text{Operator} \text{Vector} = \text{Vector}
   * @f]
   */
  template <class Operator, class Vector>
  class LinearSystem;
}

#endif
