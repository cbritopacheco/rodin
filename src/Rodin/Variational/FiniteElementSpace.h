/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTSPACE_H
#define RODIN_VARIATIONAL_FINITEELEMENTSPACE_H

#include <variant>

#include <mfem.hpp>

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

      size_t getOrder(const Geometry::Polytope& simplex) const
      {
        if (simplex.getDimension() == simplex.getMesh().getDimension())
        {
          return getHandle().GetElementOrder(simplex.getIndex());
        }
        else if (simplex.getDimension() == simplex.getMesh().getDimension() - 1)
        {
          return getHandle().GetFaceOrder(simplex.getIndex());
        }
        else
        {
          assert(false);
          return 0;
        }
      }

      /**
       * @brief Gets the vector dimensions
       */
      size_t getVectorDimension() const;

      mfem::Array<int> getEssentialTrueDOFs(const std::set<Geometry::Attribute>& bdrAttr) const;

      mfem::Array<int> getEssentialTrueDOFs(const std::set<Geometry::Attribute>& bdrAttr, size_t component) const;

      virtual size_t getSize() const = 0;

      virtual mfem::Array<int> getDOFs(const Geometry::Polytope& element) const = 0;

      virtual const Geometry::MeshBase& getMesh() const = 0;

      virtual const FiniteElementCollectionBase& getFiniteElementCollection() const = 0;

      virtual mfem::FiniteElementSpace& getHandle() const = 0;
  };

  inline
  bool operator==(const FiniteElementSpaceBase& lhs, const FiniteElementSpaceBase& rhs)
  {
    return &lhs == &rhs;
  }
}

#endif
