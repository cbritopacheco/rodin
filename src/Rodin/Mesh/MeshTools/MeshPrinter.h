/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHTOOLS_PRINTER_H
#define RODIN_MESH_MESHTOOLS_PRINTER_H

#include <string>
#include <optional>
#include <mfem.hpp>

#include "Rodin/IO/Printer.h"

#include "ForwardDecls.h"

namespace Rodin::MeshTools
{
   class MeshPrinterBase : public IO::Printer<mfem::Mesh>
   {
      public:
         MeshPrinterBase(mfem::Mesh& mesh)
            : m_mesh(mesh)
         {}

         const mfem::Mesh& getMesh() const
         {
            return m_mesh;
         }

         virtual IO::Status print(std::ostream& os) = 0;

      private:
         mfem::Mesh& m_mesh;
   };

   template <>
   class MeshPrinter<MeshFormat::MFEM> : public MeshPrinterBase
   {
      public:
         MeshPrinter(mfem::Mesh& mesh)
            : MeshPrinterBase(mesh)
         {}

         IO::Status print(std::ostream& os) override
         {
            getMesh().Print(os);
            return {true, {}};
         }
   };

   template <>
   class MeshPrinter<MeshFormat::GMSH> : public MeshPrinterBase
   {
      public:
         MeshPrinter(mfem::Mesh& mesh)
            : MeshPrinterBase(mesh)
         {}

         IO::Status print(std::ostream& os) override;
   };
}

#endif

