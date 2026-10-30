/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file P1.h
 * @brief Vertex-based continuous finite element spaces and their form operators.
 *
 * The scalar family has one degree of freedom at each mesh vertex. Shared
 * vertex coefficients identify compatible traces across neighboring cells.
 * Vector and matrix ranges repeat the scalar construction componentwise.
 *
 * ## Reference and physical shape functions
 *
 * On a simplex the reference space is affine.
 * @f$\widehat P=\mathbb P_1(\widehat K)@f$. Product cells use their
 * vertex-based product spaces that are bilinear on quadrilaterals, trilinear on
 * hexahedra and triangle-by-segment on wedges. The pyramid uses its rational
 * reference construction. The physical functions are
 * @f$\phi_i=\widehat\phi_i\circ F_K^{-1}@f$.
 *
 * The physical gradient is constant on affine simplices. Non-affine maps
 * and the other cell families generally produce varying gradients.
 * Approximation estimates require the corresponding mesh regularity,
 * mapping and solution assumptions.
 *
 * ## Usage
 *
 * Given a mesh with boundary connectivity and a source coefficient @c f,
 * the following code constructs the scalar Poisson form.
 * @code{.cpp}
 * P1 space(mesh);
 * TrialFunction u(space);
 * TestFunction v(space);
 * Problem poisson(u, v);
 * poisson = Integral(Grad(u), Grad(v)) - Integral(f, v)
 *         + DirichletBC(u, Zero());
 * @endcode
 *
 * Include @c Rodin/Variational.h for the complete form vocabulary and
 * @c Rodin/Solver.h for solvers. The runnable
 * @ref examples-pdes-poisson "Poisson tutorial" includes the mesh, solve and output.
 *
 * ## Operators
 *
 * This aggregate provides P1 derivatives, gradients, divergence, Jacobians,
 * shape-function operators, potential forms, quadrature dispatch and linear
 * elasticity integrators. Range and geometry restrictions are documented on
 * the corresponding specializations.
 *
 * @see @ref fem-ciarlet, @ref fem-spaces, @ref guides-polytopes
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
