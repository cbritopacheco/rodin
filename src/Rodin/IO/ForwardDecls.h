/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_FORWARDDECLS_H
#define RODIN_IO_FORWARDDECLS_H

#include <string>
#include <optional>
#include <boost/filesystem.hpp>

namespace Rodin::IO
{
   template <class T>
   class Loader;

   template <class T>
   class Printer;

   struct Error
   {
      std::string message;
   };

   struct Status
   {
      bool success;
      std::optional<Error> error;
   };

   enum class MeshFormat
   {
      MFEM,
      GMSH,
      MEDIT
   };

   enum class GridFunctionFormat
   {
      MFEM,
      MEDIT
   };

   template <MeshFormat fmt, class Trait>
   class MeshLoader;

   template <MeshFormat fmt, class Trait>
   class MeshPrinter;

   template <GridFunctionFormat fmt, class FEC, class Trait>
   class GridFunctionLoader;

   template <GridFunctionFormat fmt, class FEC, class Trait>
   class GridFunctionPrinter;
}

#endif

