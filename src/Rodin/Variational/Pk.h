/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Pk.h
 * @brief Arbitrary degree Lagrange (Pk) finite element space.
 *
 * This header aggregates all functionality related to the Pk finite element
 * space, which consists of continuous piecewise polynomial functions of
 * degree k over mesh elements.
 *
 * ## Mathematical Foundation
 * The Pk space of degree k is defined as:
 * @f[
 *   P_k(\mathcal{T}_h) = \{u \in C^0(\Omega) : u|_K \in \mathbb{P}_k(K) \, \forall K \in \mathcal{T}_h\}
 * @f]
 * where @f$ \mathbb{P}_k(K) @f$ denotes polynomials of degree at most k on element @f$ K @f$.
 *
 * ## Features
 * - Supports arbitrary polynomial degree k (P0, P1, P2, P3, P4, P5, P6, ...)
 * - Continuous across element boundaries (@f$ C^0 @f$ conforming)
 * - Lagrange nodal basis functions: @f$ \phi_i(x_j) = \delta_{ij} @f$
 * - Optimal convergence: @f$ O(h^{k+1}) @f$ in @f$ L^2 @f$ norm, @f$ O(h^k) @f$ in @f$ H^1 @f$ norm
 * - Supports all standard geometries: Segment, Triangle, Quadrilateral, Tetrahedron, Wedge
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
 * // P2 element (quadratic)
 * RealPkElement<2> p2(Polytope::Type::Triangle);
 * std::cout << "DOFs: " << p2.getCount() << std::endl;  // 6
 * 
 * // P3 element (cubic)
 * RealPkElement<3> p3(Polytope::Type::Segment);
 * auto basis = p3.getBasis(0);
 * auto value = basis(Math::Vector<Real>{{0.5}});
 * 
 * // Complex-valued P4 element
 * ComplexPkElement<4> p4c(Polytope::Type::Quadrilateral);
 * 
 * // Vector-valued P2 element for elasticity
 * VectorPkElement<2> p2v(Polytope::Type::Tetrahedron);
 * ```
 *
 * ## Consistency
 * - `PkElement<0>` is equivalent to `P0Element`
 * - `PkElement<1>` is equivalent to `P1Element`
 */
#ifndef RODIN_VARIATIONAL_PK_H
#define RODIN_VARIATIONAL_PK_H

#include "Pk/PkElement.h"

#endif
