/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_H
#define RODIN_SOLVER_H

#include "Solver/Solver.h"

#include "Solver/LDLT.h"

// Built-in direct solvers
#include "Solver/SparseLU.h"
#include "Solver/SparseQR.h"
#include "Solver/SimplicialLLT.h"
#include "Solver/SimplicialLDLT.h"

// Built-in iteratives solvers
#include "Solver/CG.h"
#include "Solver/LeastSquaresCG.h"

// SuiteSparse solvers
#include "Solver/UMFPack.h"

#endif
