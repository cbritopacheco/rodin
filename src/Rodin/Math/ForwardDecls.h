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
#include "Tensor.h"

#include "SparseMatrix.h"

namespace Rodin::Math
{
  template <class Src, class Dst>
  static void copy(const Src& src, Dst& dst);

  template <class Src, class Dst>
  static void duplicate(const Src& src, Dst& dst);

  template <class Y, class A, class X>
  static void axpy(Y& y, A a, const X& x);
}

#endif
