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

         virtual const MeshBase& getMesh() const = 0;

         virtual const FiniteElementCollectionBase& getFiniteElementCollection() const = 0;

         virtual mfem::FiniteElementSpace& getHandle() = 0;

         virtual const mfem::FiniteElementSpace& getHandle() const = 0;
   };

   template <class FEC>
   class FiniteElementSpace<FEC, Parallel::Trait::Serial>
      : public FiniteElementSpaceBase
   {
      public:
         FiniteElementSpace(
               Mesh<Parallel::Trait::Serial>& mesh,
               int vdim = 1, int order = 1, typename FEC::Basis basis = FEC::DefaultBasis)
            :  m_mesh(mesh),
               m_fec(order, mesh.getDimension(), basis),
               m_fes(&mesh.getHandle(), &m_fec.getHandle(), vdim)
         {}

         const MeshBase& getMesh() const override
         {
            return m_mesh;
         }

         const FEC& getFiniteElementCollection() const override
         {
            return m_fec;
         }

         mfem::FiniteElementSpace& getHandle() override
         {
            return m_fes;
         }

         const mfem::FiniteElementSpace& getHandle() const override
         {
            return m_fes;
         }

      private:
         FEC m_fec;
         Mesh<Parallel::Trait::Serial>& m_mesh;
         mfem::FiniteElementSpace m_fes;
   };
}

#endif
