/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert/MemberFunctionException.h"
#include "Adapt.h"

namespace Rodin::MMG
{
  void Adapt::adapt(MMG::Mesh& mesh, const RealGridFunction& sizeMap)
  {
    MMG5_pMesh mmgMesh = MMG::V5::Mesh::toV5(mesh);
    MMG5_pSol mmgSol = MMG::V5::Solution::create(mmgMesh, 1);
    MMG::V5::Solution::copy(sizeMap, mmgSol);

    Base::setParameters(mmgMesh);

    bool isSurface = mesh.isSurface();
    ReturnCode retcode = MMG5_STRONGFAILURE;
    switch (mmgMesh->dim)
    {
      case 2:
      {
        assert(!isSurface);
        retcode = adaptMMG2D(mmgMesh, mmgSol);
        break;
      }
      case 3:
      {
        if (isSurface)
          retcode = adaptMMGS(mmgMesh, mmgSol);
        else
          retcode = adaptMMG3D(mmgMesh, mmgSol);
        break;
      }
    }

    if (retcode != MMG5_SUCCESS)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to adapt the mesh." << Alert::Raise;
    }

    mesh = MMG::V5::Mesh::fromV5(mmgMesh);

    MMG::V5::Solution::destroy(mmgSol);
    MMG::V5::Mesh::destroy(mmgMesh);
  }

  ReturnCode Adapt::adaptMMG2D(MMG5_pMesh mesh, MMG5_pSol size)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    assert(mesh->dim == 2);
    return MMG2D_mmg2dlib(mesh, size);
  }

  ReturnCode Adapt::adaptMMG3D(MMG5_pMesh mesh, MMG5_pSol size)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    assert(mesh->dim == 3);
    return MMG3D_mmg3dlib(mesh, size);
  }

  ReturnCode Adapt::adaptMMGS(MMG5_pMesh mesh, MMG5_pSol size)
  {
    assert(mesh->np > 0 && mesh->nt > 0);
    return MMGS_mmgslib(mesh, size);
  }
}
