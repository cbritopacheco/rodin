/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file H1.h
 * @brief Positive-degree continuous Lagrange spaces and reference elements.
 *
 * H1 constructs a conforming field from degree-@f$k@f$ reference elements,
 * with @f$k\geq1@f$. Vertex, edge, face and cell degrees of freedom are
 * shared according to the mesh incidence and the family trace maps.
 * Scalar, vector and matrix ranges use the corresponding range specialization.
 *
 * ## Approximation spaces
 *
 * On affine simplices the reference space is @f$\mathbb P_k@f$; on tensor
 * cells it is @f$\mathbb Q_k@f$. Wedges combine triangular and segment
 * factors, while pyramids use their compatible rational construction.
 * The physical basis is transported by the geometric map,
 * @f$\phi_i=\widehat\phi_i\circ F_K^{-1}@f$.
 *
 * For scalar elements the local counts are:
 *
 * | Geometry | Scalar degrees of freedom |
 * | :--- | :--- |
 * | Segment | @f$k+1@f$ |
 * | Triangle | @f$(k+1)(k+2)/2@f$ |
 * | Quadrilateral | @f$(k+1)^2@f$ |
 * | Tetrahedron | @f$(k+1)(k+2)(k+3)/6@f$ |
 * | Wedge | @f$(k+1)^2(k+2)/2@f$ |
 * | Hexahedron | @f$(k+1)^3@f$ |
 * | Pyramid | @f$(k+1)(k+2)(2k+3)/6@f$ |
 *
 * Reference nodal duality, geometric mapping and global coefficient numbering
 * are separate objects. Use the element and space accessors when evaluating
 * fields or transferring their coefficients. Polynomial approximation rates
 * depend on solution regularity, mesh quality, geometry and quadrature.
 *
 * ## Usage
 *
 * A quadratic scalar space on a prepared mesh is constructed by selecting
 * its degree at compile time.
 * @code{.cpp}
 * H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
 * TrialFunction u(quadratic);
 * TestFunction v(quadratic);
 * @endcode
 *
 * @ref Variational/FiniteElementSpaces.cpp "The finite element space example"
 * constructs and samples constant, linear and quadratic fields on one mesh.
 * @ref Rodin::Variational::RealH1Element "RealH1Element<K>" exposes a scalar
 * reference element; range-specific aliases and space specializations provide
 * the other supported field types.
 *
 * @see @ref fem-ciarlet, @ref fem-spaces, @ref guides-ordering
 */
#ifndef RODIN_VARIATIONAL_H1_H
#define RODIN_VARIATIONAL_H1_H

#include "H1/H1.h"
#include "H1/H1Element.h"
#include "H1/Derivative.h"
#include "H1/Grad.h"
#include "H1/Jacobian.h"
#include "H1/Div.h"
#include "H1/ShapeFunction.h"
#include "H1/QuadratureRule.h"

#endif
