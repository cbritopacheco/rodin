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
 * finite element spaces. P1 elements are the most commonly used conforming
 * finite elements, providing continuous piecewise linear approximations.
 *
 * # Mathematical Foundation
 *
 * The P1 space consists of continuous functions that are affine on each element:
 * @f[
 *   V_h^1 = \{ v \in C^0(\Omega) : v|_K \in \mathbb{P}_1(K) \ \forall K \in \mathcal{T}_h \}
 * @f]
 *
 * On each element @f$ K @f$, functions are linear (affine):
 * @f[
 *   v(x) = a + b \cdot x, \quad x \in K
 * @f]
 *
 * ## Lagrange Property
 * P1 basis functions @f$ \phi_i @f$ are the unique continuous piecewise linear
 * functions satisfying:
 * @f[
 *   \phi_i(x_j) = \delta_{ij}
 * @f]
 * where @f$ x_j @f$ are the mesh vertices.
 *
 * ## Properties
 * - **Continuity**: @f$ H^1 @f$-conforming (continuous across elements)
 * - **Polynomial degree**: 1 (linear/affine functions)
 * - **DOFs per element**: Equal to number of vertices per element
 *   - Triangle: 3 DOFs (at vertices)
 *   - Quadrilateral: 4 DOFs (at vertices)
 *   - Tetrahedron: 4 DOFs (at vertices)
 * - **Global DOFs**: Equal to the number of mesh vertices
 * - **Support**: Each basis function has compact support on adjacent elements
 *
 * ## Applications
 * - Standard finite element method (FEM) for elliptic PDEs
 * - Parabolic PDEs (heat equation, diffusion)
 * - Linear elasticity and structural mechanics
 * - Navier-Stokes equations (velocity/pressure approximation)
 * - General second-order PDEs
 *
 * # Usage Examples
 *
 * ## Poisson Equation with P1 Elements
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * // Create mesh and P1 space
 * Mesh mesh;
 * mesh.load("domain.mesh");
 * mesh.getConnectivity().compute(1, 2);
 *
 * P1 Vh(mesh);
 *
 * // Define trial and test functions
 * TrialFunction u(Vh);
 * TestFunction  v(Vh);
 *
 * // Weak form: Find u such that (grad u, grad v) = (f, v)
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v))
 *         - Integral(f * v)
 *         + DirichletBC(u, 0.0);
 *
 * Solver::SparseLU solver;
 * problem.solve(solver);
 *
 * // Access solution
 * GridFunction<P1>& solution = problem.getSolution();
 * solution.save("solution.case");
 * @endcode
 *
 * ## Vector-Valued P1 Space (Elasticity)
 * @code{.cpp}
 * // 2D elasticity problem
 * P1 Vh(mesh, 2);  // 2 components for displacement vector
 *
 * TrialFunction u(Vh);  // Displacement: u = (u_x, u_y)
 * TestFunction  v(Vh);
 *
 * // Strain tensor
 * auto epsilon = [](auto& w) {
 *   return 0.5 * (Grad(w) + Transpose(Grad(w)));
 * };
 *
 * // Lame parameters
 * double lambda = 1.0, mu = 1.0;
 *
 * Problem problem(u, v);
 * problem = Integral(lambda * Trace(epsilon(u)) * Trace(epsilon(v))
 *                   + 2 * mu * Dot(epsilon(u), epsilon(v)))
 *         + DirichletBC(u, Zero());
 *
 * problem.solve(solver);
 * @endcode
 *
 * ## Accessing Gradient
 * @code{.cpp}
 * P1 Vh(mesh);
 * GridFunction<P1> u(Vh);
 *
 * // Gradient is piecewise constant (in P0)
 * auto grad_u = Grad(u);
 *
 * // Compute L2 norm of gradient
 * double grad_norm = std::sqrt(Integral(Dot(grad_u, grad_u)).compute());
 * @endcode
 *
 * ## Mixed Formulations
 * @code{.cpp}
 * // Stokes problem: velocity in P1, pressure in P0
 * P1 Vh(mesh, 2);  // Velocity space
 * P0 Qh(mesh);     // Pressure space (discontinuous)
 *
 * TrialFunction u(Vh), p(Qh);
 * TestFunction  v(Vh), q(Qh);
 *
 * // Weak form with velocity and pressure
 * // ... (Stokes formulation)
 * @endcode
 *
 * @see P0, FiniteElementSpace, GridFunction
 * @see TrialFunction, TestFunction, Problem
 * @see Grad, Div, Jacobian
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

