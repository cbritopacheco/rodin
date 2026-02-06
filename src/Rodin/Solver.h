/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLVER_H
#define RODIN_SOLVER_H

/**
 * @file Solver.h
 * @brief Main header for the Solver module
 *
 * This file includes all solver implementations available in Rodin. The Solver
 * module provides various algorithms for solving linear systems @f$ Ax = b @f$
 * that arise from finite element discretizations.
 *
 * Available solver categories:
 * - **Direct solvers**: Factorization-based methods (LU, Cholesky, QR)
 * - **Iterative solvers**: Krylov subspace methods (CG, BiCGSTAB)
 * - **Specialized solvers**: Least-squares and problem-specific algorithms
 * - **External library solvers**: High-performance solvers from SuiteSparse
 *
 * @see Rodin::Solver namespace for all solver classes
 */

#include "Solver/Solver.h"

#include "Solver/LDLT.h"

// Built-in direct solvers
#include "Solver/SparseLU.h"
#include "Solver/SparseQR.h"
#include "Solver/SimplicialLLT.h"
#include "Solver/SimplicialLDLT.h"

// Built-in iteratives solvers
#include "Solver/CG.h"
#include "Solver/BiCGSTAB.h"
#include "Solver/LeastSquaresCG.h"
#include "Solver/GMRES.h"

// SuiteSparse solvers
#include "Solver/UMFPack.h"
#include "Solver/SPQR.h"
#include "Solver/CHOLMOD.h"

#endif
