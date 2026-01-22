/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file H1.h
 * @brief H1-conforming Lagrange finite element space.
 *
 * This header provides access to the H1-conforming Lagrange finite element
 * family. These are high-order-stable elements with continuous (H¹) basis
 * functions across element boundaries.
 *
 * ## Mathematical Foundation
 * The H1-conforming Lagrange space of degree k is defined as:
 * @f[
 *   V_h^k = \{u \in H^1(\Omega) : u|_K \in \mathbb{P}_k(K) \, \forall K \in \mathcal{T}_h\}
 * @f]
 * where @f$ \mathbb{P}_k(K) @f$ denotes polynomials of degree at most k on element @f$ K @f$.
 *
 * ## Features
 * - Supports arbitrary polynomial degree k (K=0, 1, 2, 3, 4, 5, 6, ...)
 * - H¹-conforming (continuous across element boundaries)
 * - Lagrange nodal basis functions: @f$ \phi_i(x_j) = \delta_{ij} @f$
 * - Optimal convergence: @f$ O(h^{k+1}) @f$ in @f$ L^2 @f$ norm, @f$ O(h^k) @f$ in @f$ H^1 @f$ norm
 * - Supports all standard geometries: Segment, Triangle, Quadrilateral, Tetrahedron, Wedge
 * - High-order-stable implementation
 *
 * ## Degree of Freedom (DOF) Formulas
 * - Segment: @f$ k+1 @f$
 * - Triangle: @f$ \frac{(k+1)(k+2)}{2} @f$
 * - Quadrilateral: @f$ (k+1)^2 @f$
 * - Tetrahedron: @f$ \frac{(k+1)(k+2)(k+3)}{6} @f$
 * - Wedge: @f$ (k+1) \cdot \frac{(k+1)(k+2)}{2} @f$
 *
 * ## Usage Example
 * ```cpp
 * // H1 element of degree 2 (quadratic)
 * RealH1Element<2> h1_2(Polytope::Type::Triangle);
 * std::cout << "DOFs: " << h1_2.getCount() << std::endl;  // 6
 * 
 * // H1 element of degree 3 (cubic)
 * RealH1Element<3> h1_3(Polytope::Type::Segment);
 * auto basis = h1_3.getBasis(0);
 * auto value = basis(Math::Vector<Real>{{0.5}});
 * 
 * // Complex-valued H1 element of degree 4
 * ComplexH1Element<4> h1_c(Polytope::Type::Quadrilateral);
 * 
 * // Vector-valued H1 element for elasticity
 * VectorH1Element<2> h1_v(Polytope::Type::Tetrahedron, 3);
 *
 * // H1 finite element space of degree 2
 * H1<2> Vh(mesh);
 * TrialFunction u(Vh);
 * TestFunction v(Vh);
 * ```
 */
#ifndef RODIN_VARIATIONAL_H1_H
#define RODIN_VARIATIONAL_H1_H

#include "H1/H1.h"
#include "H1/H1Element.h"
#include "H1/Grad.h"
#include "H1/Jacobian.h"
#include "H1/Div.h"

#endif
