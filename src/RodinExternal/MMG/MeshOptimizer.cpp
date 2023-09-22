/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MeshOptimizer.h"

namespace Rodin::External::MMG
{
  int Optimizer::optimizeMMG2D(MMG5_pMesh mesh)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    MMG5_pSol sol = nullptr;
    MMG5_SAFE_CALLOC(sol, 1, MMG5_Sol,
       Alert::Exception() << "Could not allocate MMG5_Sol." << Alert::Raise);
    if (!MMG2D_Set_solSize(mesh, sol, MMG5_Vertex, 0, MMG5_Scalar))
      Alert::Exception() << "Could not set solution size." << Alert::Raise;
    MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_optim, 1);
    const ReturnCode rc = MMG2D_mmg2dlib(mesh, sol);
    if (sol->m)
      MMG5_DEL_MEM(mesh, sol->m);
    MMG5_SAFE_FREE(sol);
    return rc;
  }

  int Optimizer::optimizeMMG3D(MMG5_pMesh mesh)
  {
    assert(mesh->np > 0 && mesh->ne > 0);
    MMG5_pSol sol = nullptr;
    MMG5_SAFE_CALLOC(sol, 1, MMG5_Sol,
       Alert::Exception() << "Could not allocate MMG5_Sol." << Alert::Raise);
    if (!MMG3D_Set_solSize(mesh, sol, MMG5_Vertex, 0, MMG5_Scalar))
      Alert::Exception() << "Could not set solution size." << Alert::Raise;
    MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_optim, 1);
    const ReturnCode rc = MMG3D_mmg3dlib(mesh, sol);
    if (sol->m)
      MMG5_DEL_MEM(mesh, sol->m);
    MMG5_SAFE_FREE(sol);
    return rc;
  }

  int Optimizer::optimizeMMGS(MMG5_pMesh mesh)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    MMG5_pSol sol = nullptr;
    MMG5_SAFE_CALLOC(sol, 1, MMG5_Sol,
       Alert::Exception() << "Could not allocate MMG5_Sol." << Alert::Raise);
    if (!MMGS_Set_solSize(mesh, sol, MMG5_Vertex, 0, MMG5_Scalar))
      Alert::Exception() << "Could not set solution size." << Alert::Raise;
    MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_optim, 1);
    const ReturnCode rc = MMGS_mmgslib(mesh, sol);
    if (sol->m)
      MMG5_DEL_MEM(mesh, sol->m);
    MMG5_SAFE_FREE(sol);
    return rc;
  }
}
