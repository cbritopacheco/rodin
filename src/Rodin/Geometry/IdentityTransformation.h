/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_IDENTITYTRANSFORMATION_H
#define RODIN_GEOMETRY_IDENTITYTRANSFORMATION_H

/**
 * @file
 * @brief Identity transformation for polytopes.
 */

#include "PolytopeTransformation.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Polytope identity transformation.
   *
   * This class represents the identity transformation @f$ x(r) = r @f$ where
   * reference coordinates are equal to physical coordinates. The Jacobian
   * is the identity matrix and is constant everywhere.
   *
   * This is typically used when the reference element and physical element
   * coincide, or when no geometric transformation is needed.
   *
   * @see PolytopeTransformation, IsoparametricTransformation
   */
  class IdentityTransformation final : public PolytopeTransformation
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PolytopeTransformation;

      /**
       * @brief Constructs an identity transformation.
       * @param[in] sdim Spatial dimension
       *
       * Creates an identity transformation in @f$ \mathbb{R}^{\text{sdim}} @f$.
       */
      IdentityTransformation(size_t sdim)
        : Parent(sdim, sdim)
      {}

      /**
       * @brief Copy constructor.
       */
      IdentityTransformation(const IdentityTransformation& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       */
      IdentityTransformation(IdentityTransformation&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Gets the polynomial order of the transformation.
       * @returns 1 (linear transformation)
       */
      size_t getOrder() const override
      {
        return 1;
      }

      /**
       * @brief Gets the polynomial order of the Jacobian.
       * @returns 0 (constant Jacobian)
       */
      size_t getJacobianOrder() const override
      {
        return 0;
      }

      /**
       * @brief Applies the identity transformation.
       * @param[in] rc Reference coordinates
       * @param[out] pc Physical coordinates (equal to reference coordinates)
       *
       * Sets @f$ pc = rc @f$.
       */
      void transform(const Math::SpatialVector<Real>& rc, Math::SpatialVector<Real>& pc) const override
      {
        pc = rc;
      }

      /**
       * @brief Computes the Jacobian of the identity transformation.
       * @param[in] rc Reference coordinates (unused)
       * @param[out] res Jacobian matrix (identity matrix)
       *
       * Sets @p res to the identity matrix.
       */
      void jacobian(const Math::SpatialVector<Real>& rc, Math::SpatialMatrix<Real>& res) const override
      {
        res.setIdentity();
      }

      /**
       * @brief Creates a copy of this transformation.
       * @returns Pointer to a new IdentityTransformation object
       */
      IdentityTransformation* copy() const noexcept override
      {
        return new IdentityTransformation(*this);
      }
  };
}