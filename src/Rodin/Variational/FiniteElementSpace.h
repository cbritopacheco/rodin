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
      constexpr
      FiniteElementSpaceBase() = default;

      constexpr
      FiniteElementSpaceBase(const FiniteElementSpaceBase&) = default;

      constexpr
      FiniteElementSpaceBase(FiniteElementSpaceBase&&) = default;

      constexpr
      FiniteElementSpaceBase& operator=(FiniteElementSpaceBase&&) = default;

      void update();

      /**
       * @returns Order of the highest dimensional finite element.
       */
      size_t getOrder() const;

      size_t getNumberOfDofs() const;

      /**
       * @brief Gets the vector dimensions
       */
      size_t getVectorDimension() const;

      mfem::Array<int> getEssentialTrueDOFs(const std::set<int>& bdrAttr) const;

      mfem::Array<int> getEssentialTrueDOFs(const std::set<int>& bdrAttr, int component) const;

      virtual size_t getSize() const = 0;

      virtual mfem::Array<int> getDOFs(
          const Geometry::Simplex& element) const = 0;

      virtual FiniteElement getFiniteElement(
          const Geometry::Simplex& element) const = 0;

      virtual Geometry::MeshBase& getMesh() = 0;

      virtual const Geometry::MeshBase& getMesh() const = 0;

      virtual const FiniteElementCollectionBase& getFiniteElementCollection() const = 0;

      virtual mfem::FiniteElementSpace& getHandle() = 0;

      virtual const mfem::FiniteElementSpace& getHandle() const = 0;
  };
}

#endif
