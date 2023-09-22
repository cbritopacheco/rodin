#include "Adapt.h"

namespace Rodin::External::MMG
{
  void Adapt::adapt(MMG::Mesh& mesh, const ScalarGridFunction& sizeMap)
  {
    MMG5_pMesh mmgMesh = rodinToMesh(mesh);
    MMG5_pSol mmgSol = createSolution(mmgMesh, 1);
    copySolution(sizeMap, mmgSol);

    MMG5::setParameters(mmgMesh);

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
      Alert::Exception()
        << "Failed to adapt the mesh."
        << Alert::Raise;
    }

    mesh = meshToRodin(mmgMesh);

    destroySolution(mmgSol);
    destroyMesh(mmgMesh);
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
