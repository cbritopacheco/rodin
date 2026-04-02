/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file P1.h
 * @brief Piecewise linear (P1) finite element space.
 *
 * This header aggregates all functionality related to the P1 finite element
 * space, which consists of continuous piecewise linear functions over mesh
 * elements. P1 is the most commonly used conforming finite element space.
 *
 * ## Mathematical Foundation
 * The P1 space is defined as:
 * @f[
 *   P_1(\mathcal{T}_h) = \{u \in C^0(\Omega) : u|_K \in \mathbb{P}_1(K) \, \forall K \in \mathcal{T}_h\}
 * @f]
 * where @f$ \mathbb{P}_1(K) @f$ denotes linear polynomials on element @f$ K @f$.
 *
 * ## Features
 * - One degree of freedom per mesh vertex (nodal values)
 * - Continuous across element boundaries (@f$ C^0 @f$ conforming)
 * - Gradient is constant within each element
 * - Optimal convergence: @f$ O(h^2) @f$ in @f$ L^2 @f$ norm for smooth solutions
 * - Basis functions: Lagrange nodal basis @f$ \phi_i(x_j) = \delta_{ij} @f$
 *
 * ## Usage Example
 *
 * ```cpp
 * P1 Vh(mesh);           // Define P1 space
 * TrialFunction u(Vh);   // P1 trial function
 * TestFunction v(Vh);    // P1 test function
 * 
 * // Poisson problem: -Δu = f
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) - Integral(f, v);
 * ```
 *
 * ## Supported Operations
 * - Gradient (Grad)
 * - Divergence (Div) for vector-valued functions
 * - Jacobian transformations
 * - Potential recovery
 * - Linear elasticity specializations
 */
#ifndef RODIN_VARIATIONAL_P1_H
#define RODIN_VARIATIONAL_P1_H

#include "P1/P1.h"
#include "P1/Div.h"
#include "P1/Grad.h"
#include "P1/Jacobian.h"
#include "P1/Potential.h"
#include "P1/P1Element.h"
#include "P1/Derivative.h"
#include "P1/QuadratureRule.h"
#include "Rodin/Solid/Linear/P1/LinearElasticityIntegral.h"
#include "P1/ShapeFunction.h"

#endif
