/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MATH_H
#define RODIN_MATH_H

/**
 * @file
 * @brief Top level include for the Rodin::Math module.
 *
 * The Math module provides mathematical types and operations used throughout
 * Rodin, including vectors, matrices, sparse matrices, and mathematical constants.
 * These types are built on top of the Eigen library for high-performance
 * linear algebra operations.
 *
 * @see Rodin::Math
 */

#include "Math/ForwardDecls.h"
#include "Math/Traits.h"
#include "Math/Common.h"
#include "Math/Constants.h"
#include "Math/Vector.h"
#include "Math/Matrix.h"
#include "Math/SparseMatrix.h"
#include "Math/SpatialVector.h"
#include "Math/SpatialMatrix.h"
#include "Math/Unit.h"
#include "Math/Rad.h"
#include "Math/Deg.h"
#include "Math/LinearSystem.h"
#include "Math/RootFinding/NewtonRaphson.h"
#include "Math/RungeKutta/RK2.h"
#include "Math/RungeKutta/RK4.h"

#endif
