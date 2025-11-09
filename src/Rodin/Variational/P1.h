/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_H
#define RODIN_VARIATIONAL_P1_H

/**
 * @file
 * @brief Piecewise linear (P1) finite element space.
 *
 * This file includes all components related to P1 (piecewise linear/affine)
 * finite element spaces. P1 elements use linear basis functions on each element:
 * @f[
 *   \phi_i(x) = a_i + b_i \cdot x
 * @f]
 * where @f$ a_i @f$ and @f$ b_i @f$ are determined by the nodal values.
 *
 * ## Properties
 * - **Continuity**: @f$ H^1 @f$-conforming (continuous across elements)
 * - **Polynomial degree**: 1
 * - **DOFs per element**: 
 *   - Triangle: 3 (vertices)
 *   - Tetrahedron: 4 (vertices)
 * - **Applications**: Standard FEM for elliptic/parabolic PDEs
 *
 * ## Mathematical Foundation
 * P1 basis functions satisfy the Lagrange property:
 * @f[
 *   \phi_i(x_j) = \delta_{ij}
 * @f]
 * where @f$ x_j @f$ are the mesh vertices.
 */

#include "P1/P1.h"
#include "P1/Div.h"
#include "P1/Grad.h"
#include "P1/Jacobian.h"
#include "P1/Potential.h"
#include "P1/P1Element.h"
#include "P1/Derivative.h"
#include "P1/QuadratureRule.h"
#include "P1/LinearElasticity.h"

#endif

