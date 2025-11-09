/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H
#define RODIN_VARIATIONAL_H

/**
 * @file
 * @brief Main header for the Variational module - finite element formulations and problem solving.
 *
 * This is the top-level include file for the Rodin::Variational namespace, providing
 * comprehensive facilities for defining and solving variational problems in finite
 * element analysis.
 *
 * # Overview
 *
 * The Variational module implements the mathematical framework for finite element
 * methods, enabling users to:
 * - Define variational formulations of PDEs
 * - Construct bilinear and linear forms
 * - Specify boundary conditions
 * - Assemble and solve discrete systems
 * - Manipulate grid functions (solutions)
 *
 * # Mathematical Foundation
 *
 * ## Weak Formulation
 * The module supports standard weak formulations: Find @f$ u \in V_h @f$ such that
 * @f[
 *   a(u,v) = l(v) \quad \forall v \in V_h
 * @f]
 * where:
 * - @f$ a(u,v) @f$ is a bilinear form
 * - @f$ l(v) @f$ is a linear form
 * - @f$ V_h @f$ is a finite element space
 *
 * ## Discrete System
 * This leads to the discrete system @f$ Au = b @f$ where:
 * @f[
 *   A_{ij} = a(\phi_j, \psi_i), \quad b_i = l(\psi_i)
 * @f]
 * with @f$ \{\phi_j\} @f$ trial basis functions and @f$ \{\psi_i\} @f$ test basis functions.
 *
 * # Usage Example
 *
 * ## Poisson Equation
 * Solving @f$ -\Delta u = f @f$ in @f$ \Omega @f$ with @f$ u = 0 @f$ on @f$ \partial\Omega @f$:
 *
 * @code{.cpp}
 * using namespace Rodin;
 * using namespace Rodin::Variational;
 *
 * // Create mesh
 * Mesh mesh;
 * mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16);
 * mesh.getConnectivity().compute(1, 2);
 *
 * // Define finite element space
 * P1 Vh(mesh);
 *
 * // Define trial and test functions
 * TrialFunction u(Vh);
 * TestFunction v(Vh);
 *
 * // Build variational problem: Find u such that
 * // Integral(Grad(u), Grad(v)) = Integral(f * v) for all v
 * Problem problem(u, v);
 * problem = Integral(Grad(u), Grad(v)) 
 *         - Integral(f * v)
 *         + DirichletBC(u, 0.0);
 *
 * // Solve the system
 * Solver::SparseLU solver;
 * problem.solve(solver);
 *
 * // Access the solution
 * GridFunction& solution = problem.getSolution();
 * @endcode
 *
 * ## Linear Elasticity
 * @code{.cpp}
 * // Vector-valued finite element space
 * P1 Vh(mesh, 2);  // 2D elasticity
 *
 * TrialFunction u(Vh);  // Displacement
 * TestFunction  v(Vh);  // Test function
 *
 * // Strain tensor
 * auto epsilon = [](auto& w) { 
 *   return 0.5 * (Grad(w) + Transpose(Grad(w))); 
 * };
 *
 * // Lam\u00e9 parameters
 * double lambda = 1.0, mu = 1.0;
 *
 * // Weak form: a(u,v) = l(v)
 * Problem problem(u, v);
 * problem = Integral(lambda * Trace(epsilon(u)) * Trace(epsilon(v))
 *                   + 2 * mu * Dot(epsilon(u), epsilon(v)))
 *         - Integral(f, v)
 *         + DirichletBC(u, Zero());
 *
 * problem.solve(solver);
 * @endcode
 *
 * # Key Components
 *
 * ## Finite Element Spaces
 * - P0: Piecewise constant elements (discontinuous)
 * - P1: Piecewise linear elements (continuous)
 * - Higher-order spaces available
 *
 * ## Forms and Operators
 * - **Bilinear forms**: Integral(u, v), BoundaryIntegral(u, v)
 * - **Linear forms**: Integral(f * v), BoundaryIntegral(g * v)
 * - **Differential operators**: Grad, Div, Jacobian, Derivative
 * - **Discontinuous Galerkin**: Jump, Average, face integrals
 *
 * ## Boundary Conditions
 * - DirichletBC: Essential boundary conditions
 * - PeriodicBC: Periodic boundary conditions  
 * - Natural BCs: Through boundary integrals
 *
 * ## Problem Assembly
 * - Problem: Main problem class
 * - SparseProblem: Sparse matrix systems (default)
 * - DenseProblem: Dense matrix systems
 *
 * @see @ref RodinVariational
 * @see Problem, BilinearForm, LinearForm, GridFunction
 * @see TrialFunction, TestFunction, FiniteElementSpace
 */

#include "Variational/ForwardDecls.h"

#include "Variational/GridFunction.h"
#include "Variational/FiniteElementSpace.h"

#include "Variational/ShapeFunction.h"
#include "Variational/TrialFunction.h"
#include "Variational/TestFunction.h"

#include "Variational/Component.h"
#include "Variational/Restriction.h"

#include "Variational/LinearForm.h"
#include "Variational/BilinearForm.h"

#include "Variational/Zero.h"
#include "Variational/Dot.h"
#include "Variational/Pow.h"
#include "Variational/Sqrt.h"
#include "Variational/Sum.h"
#include "Variational/Mult.h"
#include "Variational/Minus.h"
#include "Variational/Division.h"
#include "Variational/UnaryMinus.h"

#include "Variational/Abs.h"
#include "Variational/Exp.h"
#include "Variational/Conjugate.h"

#include "Variational/EQ.h"
#include "Variational/GT.h"
#include "Variational/LT.h"
#include "Variational/GEQ.h"
#include "Variational/LEQ.h"
#include "Variational/AND.h"
#include "Variational/OR.h"

#include "Variational/Div.h"
#include "Variational/Grad.h"
#include "Variational/FaceNormal.h"
#include "Variational/BoundaryNormal.h"
#include "Variational/Jacobian.h"
#include "Variational/Derivative.h"
#include "Variational/Jump.h"
#include "Variational/Average.h"

#include "Variational/Trace.h"
#include "Variational/Transpose.h"
#include "Variational/IdentityMatrix.h"

#include "Variational/Max.h"
#include "Variational/Min.h"

#include "Variational/Re.h"
#include "Variational/Im.h"

#include "Variational/Sine.h"
#include "Variational/Cosine.h"
#include "Variational/Sinh.h"
#include "Variational/Cosh.h"
#include "Variational/Tangent.h"

#include "Variational/Frobenius.h"

#include "Variational/Integral.h"
#include "Variational/FaceIntegral.h"
#include "Variational/BoundaryIntegral.h"
#include "Variational/InterfaceIntegral.h"
#include "Variational/Problem.h"
#include "Variational/DenseProblem.h"

#include "Variational/RealFunction.h"
#include "Variational/VectorFunction.h"
#include "Variational/MatrixFunction.h"
#include "Variational/BooleanFunction.h"

#include "Variational/P0.h"
#include "Variational/P1.h"

#include "Variational/DirichletBC.h"
#include "Variational/PeriodicBC.h"

#include "Variational/Flow.h"
#include "Variational/Potential.h"

#include "Variational/F.h"

#endif
