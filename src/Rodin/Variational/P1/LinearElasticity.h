#ifndef RODIN_VARIATIONAL_P1_LINEARELASTICITY_H
#define RODIN_VARIATIONAL_P1_LINEARELASTICITY_H

/**
 * @file LinearElasticity.h
 * @brief P1 linear elasticity integrators.
 *
 * This file aggregates linear elasticity functionality for P1 vector spaces,
 * including specialized bilinear form integrators for elasticity problems.
 *
 * ## Linear Elasticity Formulation
 * The linear elasticity bilinear form is:
 * @f[
 *   a(\mathbf{u}, \mathbf{v}) = \int_\Omega 2\mu \, \boldsymbol{\varepsilon}(\mathbf{u}) : \boldsymbol{\varepsilon}(\mathbf{v}) + \lambda (\nabla \cdot \mathbf{u})(\nabla \cdot \mathbf{v}) \, dx
 * @f]
 * where:
 * - @f$ \boldsymbol{\varepsilon}(\mathbf{u}) = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T) @f$ is the strain tensor
 * - @f$ \mu @f$ is the shear modulus (Lamé's second parameter)
 * - @f$ \lambda @f$ is Lamé's first parameter
 *
 * @see LinearElasticityIntegral, P1
 */

#include "LinearElasticity/LinearElasticityIntegral.h"

#endif

