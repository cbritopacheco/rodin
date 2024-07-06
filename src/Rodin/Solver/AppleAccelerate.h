/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_APPLEACCELERATE_H
#define RODIN_SOLVER_APPLEACCELERATE_H

#ifdef RODIN_USE_APPLE_ACCELERATE

#include <Eigen/AppleAccelerate>

namespace Rodin::Solver
{
  template <>
  class AppleAccelerateLDLT<Math::SparseMatrix<Real>, Math::Vector<Real>> final
    : public SolverBase<Math::SparseMatrix<Real>, Math::Vector<Real>>
  {

  };
}

#endif
#endif
