/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_FORWARDDECLS_H
#define RODIN_MESH_FORWARDDECLS_H

#include "Rodin/Parallel/ForwardDecls.h"

namespace Rodin
{
   class MeshBase;

   template <Parallel::Trait ParallelOrSerial = Parallel::Trait::Serial>
   class Mesh;

   class SubMesh;
}

#endif
