/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_H
#define RODIN_VARIATIONAL_P0_H

/**
 * @file
 * @brief Piecewise constant (P0) finite element space.
 *
 * This file includes all components related to P0 (piecewise constant)
 * finite element spaces. P0 elements are the simplest finite elements where
 * basis functions are constant on each mesh element:
 * @f[
 *   \phi_K(x) = \begin{cases}
 *     1 & \text{if } x \in K \\
 *     0 & \text{otherwise}
 *   \end{cases}
 * @f]
 *
 * ## Properties
 * - **Continuity**: Discontinuous across elements
 * - **Polynomial degree**: 0
 * - **DOFs per element**: 1 per element
 * - **Applications**: Discontinuous Galerkin, finite volumes, material properties
 */

#include "P0/P0.h"
#include "P0/Grad.h"
#include "P0/P0Element.h"
#include "P0/GridFunction.h"

#endif
