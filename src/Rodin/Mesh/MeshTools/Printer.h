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

#include "ForwardDecls.h"

namespace Rodin::MeshTools
{
   class PrinterBase
   {
      public:
         PrinterBase(mfem::Mesh& mesh)
            : m_mesh(mesh)
         {}

         struct Error
         {
            std::string message;
         };

         struct Status
         {
            bool success;
            std::optional<Error> error;
         };

         const mfem::Mesh& getMesh() const
         {
            return m_mesh;
         }

         virtual Status print(std::ostream& os) = 0;

      private:
         mfem::Mesh& m_mesh;
   };

   template <>
   class Printer<MeshFormat::MFEM> : public PrinterBase
   {
      public:
         Printer(mfem::Mesh& mesh)
            : PrinterBase(mesh)
         {}

         Status print(std::ostream& os) override
         {
            getMesh().Print(os);
            return {true, {}};
         }
   };

   template <>
   class Printer<MeshFormat::GMSH> : public PrinterBase
   {
      public:
         Printer(mfem::Mesh& mesh)
            : PrinterBase(mesh)
         {}

         Status print(std::ostream& os) override;
   };
}

#endif

