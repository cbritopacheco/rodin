/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHPRINTER_H
#define RODIN_MESH_MESHPRINTER_H

#include <string>
#include <optional>
#include <mfem.hpp>

#include "Rodin/Mesh.h"
#include "Rodin/IO/Printer.h"

#include "ForwardDecls.h"

namespace Rodin::IO
{
   template <class Trait>
   class MeshPrinterBase : public IO::Printer<Rodin::Mesh<Trait>>
   {
      public:
         MeshPrinterBase(const Rodin::Mesh<Traits::Serial>& mesh)
            : m_mesh(mesh)
         {}

      protected:
         const Rodin::Mesh<Traits::Serial>& getObject() const override
         {
            return m_mesh;
         }

      private:
         const Rodin::Mesh<Traits::Serial>& m_mesh;
   };

   template <>
   class MeshPrinter<FileFormat::MFEM, Traits::Serial>
      : public MeshPrinterBase<Traits::Serial>
   {
      public:
         MeshPrinter(const Rodin::Mesh<Traits::Serial>& mesh)
            : MeshPrinterBase(mesh)
         {}

         void print(std::ostream& os) override
         {
            getObject().getHandle().Print(os);
         }
   };

   template <>
   class MeshPrinter<FileFormat::GMSH, Traits::Serial>
      : public MeshPrinterBase<Traits::Serial>
   {
      public:
         MeshPrinter(const Rodin::Mesh<Traits::Serial>& mesh)
            : MeshPrinterBase(mesh)
         {}

         void print(std::ostream& os) override;
   };

   template <>
   class MeshPrinter<FileFormat::MEDIT, Traits::Serial>
      : public MeshPrinterBase<Traits::Serial>
   {
      public:
         MeshPrinter(const Rodin::Mesh<Traits::Serial>& mesh)
            : MeshPrinterBase(mesh)
         {}

         void print(std::ostream& os) override;
   };
}

#endif

