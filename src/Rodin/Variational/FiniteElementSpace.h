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
   /**
    * @brief Abstract base class for finite element spaces.
    */
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

         virtual bool isParallel() const = 0;

         virtual MeshBase& getMesh() = 0;

         virtual const MeshBase& getMesh() const = 0;

         virtual const FiniteElementCollectionBase& getFiniteElementCollection() const = 0;

         virtual mfem::FiniteElementSpace& getHandle() = 0;

         virtual const mfem::FiniteElementSpace& getHandle() const = 0;
   };

   /**
    * @brief Represents a finite element space defined on some element
    * collection.
    * @tparam FEC Finite element collection
    */
   template <class FEC>
   class FiniteElementSpace<FEC, Traits::Serial>
      : public FiniteElementSpaceBase
   {
      public:
         /**
          * @brief Constructs a finite element space supported on the give
          * mesh.
          * @param[in] mesh Support of finite element space
          * @param[in] vdim Vector dimensions of the finite element functions
          * @param[in] order Element order
          * @param[in] basis Basis of the finite element space
          */
         FiniteElementSpace(
               Mesh<Traits::Serial>& mesh,
               int vdim = 1, int order = 1, typename FEC::Basis basis = FEC::DefaultBasis)
            :  m_mesh(mesh),
               m_fec(order, mesh.getDimension(), basis),
               m_fes(&mesh.getHandle(), &m_fec.getHandle(), vdim)
         {}

         bool isParallel() const override
         {
            return false;
         }

         Mesh<Traits::Serial>& getMesh() override
         {
            return m_mesh;
         }

         const Mesh<Traits::Serial>& getMesh() const override
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
         Mesh<Traits::Serial>& m_mesh;
         mfem::FiniteElementSpace m_fes;
   };
}

#endif
