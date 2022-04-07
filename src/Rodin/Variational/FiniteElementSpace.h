/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FINITEELEMENTSPACE_H
#define RODIN_VARIATIONAL_FINITEELEMENTSPACE_H

#include <mfem.hpp>

#include "Rodin/Mesh.h"
#include "Rodin/Utility.h"

#include "ForwardDecls.h"

namespace Rodin::Variational
{
   class FiniteElementSpaceBase
   {
      public:
         void update();

         /**
          * @returns Order of the highest dimensional finite element.
          */
         int getOrder() const;

         int getNumberOfDofs() const;

         /**
          * @brief Gets the vector dimensions
          */
         int getVectorDimension() const;

         mfem::Array<int> getEssentialTrueDOFs(const std::set<int>& bdrAttr) const;

         mfem::Array<int> getEssentialTrueDOFs(const std::set<int>& bdrAttr, int component) const;

         virtual MeshBase& getMesh() = 0;

         virtual const MeshBase& getMesh() const = 0;

         virtual mfem::FiniteElementSpace& getFES() = 0;

         virtual const mfem::FiniteElementSpace& getFES() const = 0;

         virtual mfem::FiniteElementCollection& getFEC() = 0;

         virtual const mfem::FiniteElementCollection& getFEC() const = 0;
   };
}

#endif
