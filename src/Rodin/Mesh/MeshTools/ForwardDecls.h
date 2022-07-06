/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHTOOLS_FORWARDDECLS_H
#define RODIN_MESH_MESHTOOLS_FORWARDDECLS_H

#include "../ForwardDecls.h"

namespace Rodin::MeshTools
{
   class LoaderBase;

   template <MeshFormat fmt>
   class Loader;

   class PrinterBase;

   template <MeshFormat fmt>
   class Printer;
}

#endif

