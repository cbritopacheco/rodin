/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearElasticity.h
 * @brief Linear elasticity formulations and integrators.
 *
 * This header aggregates functionality for linear elasticity problems, which
 * model the deformation of elastic bodies under small strains.
 *
 * ## Mathematical Foundation
 * Linear elasticity is governed by the equilibrium equation:
 * @f[
 *   -\nabla \cdot \boldsymbol{\sigma}(\mathbf{u}) = \mathbf{f}
 * @f]
 * where:
 * - @f$ \mathbf{u} @f$ is the displacement field
 * - @f$ \boldsymbol{\sigma} @f$ is the stress tensor
 * - @f$ \mathbf{f} @f$ is the body force
 *
 * ## Constitutive Relations
 * For isotropic linear elastic materials:
 * @f[
 *   \boldsymbol{\sigma} = \lambda \text{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2\mu \boldsymbol{\varepsilon}
 * @f]
 * where:
 * - @f$ \boldsymbol{\varepsilon} = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T) @f$ is the strain tensor
 * - @f$ \lambda, \mu @f$ are the Lamé parameters
 *
 * ## Weak Formulation
 * Find @f$ \mathbf{u} \in V @f$ such that:
 * @f[
 *   \int_\Omega \boldsymbol{\sigma}(\mathbf{u}) : \boldsymbol{\varepsilon}(\mathbf{v}) \, dx = \int_\Omega \mathbf{f} \cdot \mathbf{v} \, dx
 * @f]
 * for all @f$ \mathbf{v} \in V @f$.
 *
 * ## Usage Example
 * ```cpp
 * // Define vector-valued P1 space (dimension = spatial dimension)
 * P1 Vh(mesh, mesh.getSpaceDimension());
 * 
 * TrialFunction u(Vh);
 * TestFunction v(Vh);
 * 
 * // Lamé parameters
 * Real lambda = 1.0, mu = 1.0;
 * 
 * // Linear elasticity bilinear form
 * auto a = LinearElasticityIntegral(lambda, mu, u, v);
 * ```
 */
#ifndef RODIN_VARIATIONAL_LINEARELASTICITY_H
#define RODIN_VARIATIONAL_LINEARELASTICITY_H

#include "LinearElasticity/LinearElasticityIntegral.h"

#endif
