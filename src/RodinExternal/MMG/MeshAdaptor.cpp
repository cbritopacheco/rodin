/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MeshAdaptor.h"

#include "MMG5.h"

namespace Rodin::External::MMG
{
  int MeshAdaptor::adaptMMG3D(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    assert(mesh->np > 0 && mesh->ne > 0);
    return MMG3D_mmg3dlib(mesh, sol);
  }
}
