/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SPARSEPROBLEM_H
#define RODIN_VARIATIONAL_SPARSEPROBLEM_H

#include "Problem.h"

namespace Rodin::Variational
{
  template <class TrialFES, class TestFES>
  using SparseProblem = Problem<TrialFES, TestFES, Math::SparseMatrix<Real>, Math::Vector<Real>>;
}

#endif


