/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_IDENTITYTRANSFORMATION_H
#define RODIN_GEOMETRY_IDENTITYTRANSFORMATION_H

#include "Rodin/Geometry/Simplex.h"

#include "SimplexTransformation.h"

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
        : m_sdim(sdim)
      {}

      IdentityTransformation(const IdentityTransformation& other)
        : Parent(other),
          m_sdim(other.m_sdim)
      {}

      IdentityTransformation(IdentityTransformation&& other)
        : Parent(std::move(other)),
          m_sdim(std::move(other.m_sdim))
      {}

      inline
      Math::Vector<Real> transform(const Math::Vector<Real>& rc) const override
      {
        return rc;
      }

      inline
      Math::Matrix<Real> jacobian(const Math::Vector<Real>& rc) const override
      {
        return Math::Matrix<Real>::Identity(m_sdim == 0 ? 1 : m_sdim, rc.size());
      }

    private:
      const size_t m_sdim;
  };
}

#endif

