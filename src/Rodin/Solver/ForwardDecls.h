/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_FORWARDDECLS_H
#define RODIN_SOLVER_FORWARDDECLS_H

namespace Rodin::Solver
{
   template <class OperatorType, class VectorType>
   class UMFPack;

   template <class OperatorType, class VectorType>
   class CG;

   template <class OperatorType, class VectorType>
   class SolverBase;
}

#endif
