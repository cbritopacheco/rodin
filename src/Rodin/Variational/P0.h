/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file P0.h
 * @brief Piecewise constant (P0) finite element space.
 *
 * This header aggregates all functionality related to the P0 finite element
 * space, which consists of piecewise constant functions over mesh elements.
 *
 * ## Mathematical Foundation
 * The P0 space is defined as:
 * @f[
 *   P_0(\mathcal{T}_h) = \{u \in L^2(\Omega) : u|_K \in \mathbb{P}_0(K) \, \forall K \in \mathcal{T}_h\}
 * @f]
 * where @f$ \mathbb{P}_0(K) @f$ denotes constant functions on element @f$ K @f$.
 *
 * ## Features
 * - One degree of freedom per element (at element center)
 * - Discontinuous across element boundaries
 * - Suitable for DG methods and element-wise constant approximations
 * - Gradient of P0 functions is zero within each element
 *
 * ## Usage Example
 * ```cpp
 * P0 Vh(mesh);           // Define P0 space
 * GridFunction<P0> u(Vh); // P0 grid function
 * ```
 */
#ifndef RODIN_VARIATIONAL_P0_H
#define RODIN_VARIATIONAL_P0_H

#include "P0/P0.h"
#include "P0/Grad.h"
#include "P0/P0Element.h"
#include "P0/GridFunction.h"
#include "P0/ShapeFunction.h"

#endif
