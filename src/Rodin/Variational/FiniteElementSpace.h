/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTSPACE_H
#define RODIN_VARIATIONAL_FINITEELEMENTSPACE_H

#include <variant>

#include "Rodin/Types.h"
#include "Rodin/Utility.h"
#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"
#include "FiniteElement.h"

namespace Rodin::Variational
{
  class FiniteElementSpaceBase
  {
    public:
      FiniteElementSpaceBase() = default;

      FiniteElementSpaceBase(const FiniteElementSpaceBase&) = default;

      FiniteElementSpaceBase(FiniteElementSpaceBase&&) = default;

      FiniteElementSpaceBase& operator=(FiniteElementSpaceBase&&) = default;

      virtual size_t getSize() const = 0;

      virtual size_t getVectorDimension() const = 0;

      virtual const Geometry::MeshBase& getMesh() const = 0;

      virtual Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const = 0;
  };

  inline
  bool operator==(const FiniteElementSpaceBase& lhs, const FiniteElementSpaceBase& rhs)
  {
    return &lhs == &rhs;
  }
}

#endif
