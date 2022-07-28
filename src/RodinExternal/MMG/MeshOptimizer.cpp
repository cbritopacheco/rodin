/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MeshOptimizer.h"

namespace Rodin::External::MMG
{
  int MeshOptimizer::optimizeMMG2D(MMG5_pMesh mesh)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    MMG5_pSol sol = nullptr;
    MMG5_SAFE_CALLOC(sol, 1, MMG5_Sol,
        Alert::Exception("Could not allocate MMG5_Sol.").raise());
    if (!MMG2D_Set_solSize(mesh, sol, MMG5_Vertex, 0, MMG5_Scalar))
    Alert::Exception("Could not set solution size.").raise();
    MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_optim, 1);
    return MMG2D_mmg2dlib(mesh, sol);
  }

  int MeshOptimizer::optimizeMMG3D(MMG5_pMesh mesh)
  {
    assert(mesh->np > 0 && mesh->ne > 0);
    MMG5_pSol sol = nullptr;
    MMG5_SAFE_CALLOC(sol, 1, MMG5_Sol,
        Alert::Exception("Could not allocate MMG5_Sol.").raise());
    if (!MMG3D_Set_solSize(mesh, sol, MMG5_Vertex, 0, MMG5_Scalar))
    Alert::Exception("Could not set solution size.").raise();
    MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_optim, 1);
    return MMG3D_mmg3dlib(mesh, sol);
  }

  int MeshOptimizer::optimizeMMGS(MMG5_pMesh mesh)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    MMG5_pSol sol = nullptr;
    MMG5_SAFE_CALLOC(sol, 1, MMG5_Sol,
        Alert::Exception("Could not allocate MMG5_Sol.").raise());
    if (!MMGS_Set_solSize(mesh, sol, MMG5_Vertex, 0, MMG5_Scalar))
    Alert::Exception("Could not set solution size.").raise();
    MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_optim, 1);
    return MMGS_mmgslib(mesh, sol);
  }
}
