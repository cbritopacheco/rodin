/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_IDENTITYTRANSFORMATION_H
#define RODIN_GEOMETRY_IDENTITYTRANSFORMATION_H


#include "PolytopeTransformation.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Polytope Identity transformation.
   */
  class IdentityTransformation final : public PolytopeTransformation
  {
    public:
      using Parent = PolytopeTransformation;

      IdentityTransformation(size_t sdim)
        : Parent(sdim, sdim)
      {}

      IdentityTransformation(const IdentityTransformation& other)
        : Parent(other)
      {}

      IdentityTransformation(IdentityTransformation&& other)
        : Parent(std::move(other))
      {}

      size_t getOrder() const override
      {
        return 1;
      }

      size_t getJacobianOrder() const override
      {
        return 0;
      }

      void transform(const Math::SpatialVector<Real>& rc, Math::SpatialVector<Real>& pc) const override
      {
        pc = rc;
      }

      void jacobian(const Math::SpatialVector<Real>& rc, Math::SpatialMatrix<Real>& res) const override
      {
        res.setIdentity();
      }

      IdentityTransformation* copy() const noexcept override
      {
        return new IdentityTransformation(*this);
      }
  };
}

#endif

