/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "CG.h"

namespace Rodin::Solver
{
  template class CG<Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>>;
  template class CG<
    Math::LinearSystem<Math::SparseMatrix<Complex>, Math::Vector<Complex>>>;
}
