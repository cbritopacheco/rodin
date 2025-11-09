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
 * basis functions are constant on each mesh element.
 *
 * # Mathematical Foundation
 *
 * The P0 space consists of functions that are constant on each element:
 * @f[
 *   V_h^0 = \{ v \in L^2(\Omega) : v|_K \in \mathbb{P}_0(K) \ \forall K \in \mathcal{T}_h \}
 * @f]
 *
 * Each basis function @f$ \phi_K @f$ is the characteristic function of element @f$ K @f$:
 * @f[
 *   \phi_K(x) = \begin{cases}
 *     1 & \text{if } x \in K \\
 *     0 & \text{otherwise}
 *   \end{cases}
 * @f]
 *
 * ## Properties
 * - **Continuity**: Discontinuous across elements (no inter-element continuity)
 * - **Polynomial degree**: 0 (constant functions)
 * - **DOFs per element**: 1 degree of freedom per element
 * - **Global DOFs**: Equal to the number of mesh elements
 * - **Local basis**: Single basis function per element (value 1)
 *
 * ## Applications
 * - Discontinuous Galerkin (DG) methods
 * - Finite volume methods
 * - Material property fields (density, conductivity, etc.)
 * - Characteristic functions and indicators
 * - Error indicators for adaptivity
 *
 * # Usage Example
 *
 * ## Scalar P0 Space
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * Mesh mesh;
 * mesh.load("domain.mesh");
 *
 * // Create P0 finite element space
 * P0 Qh(mesh);
 *
 * // Grid function for material properties
 * GridFunction<P0> rho(Qh);  // Density field
 *
 * // Set constant value per element
 * rho.getDOFs() << 1.0, 2.0, 1.5, ...;  // One value per element
 * @endcode
 *
 * ## Discontinuous Galerkin Example
 * @code{.cpp}
 * // P0 spaces for DG formulation
 * P0 Qh(mesh);
 *
 * TrialFunction u(Qh);
 * TestFunction  v(Qh);
 *
 * // DG bilinear form with jumps and averages
 * Problem problem(u, v);
 * problem = Integral(u * v)  // Mass term
 *         + FaceIntegral(Jump(u) * Average(v))  // DG coupling
 *         + BoundaryIntegral(u * v);
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Vector-Valued P0 Space
 * @code{.cpp}
 * // P0 space with 2 components (2D vector field)
 * P0 Qh(mesh, 2);
 *
 * GridFunction<P0> velocity(Qh);
 * // velocity has 2 * numElements DOFs
 * @endcode
 *
 * @see P1, FiniteElementSpace, GridFunction
 * @see TrialFunction, TestFunction, Problem
 */

#include "P0/P0.h"
#include "P0/Grad.h"
#include "P0/P0Element.h"
#include "P0/GridFunction.h"

#endif
