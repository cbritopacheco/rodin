#include "ImplicitDomainMesher.h"

namespace Rodin::External::MMG
{
  ImplicitDomainMesher& ImplicitDomainMesher::surface(bool meshTheSurface)
  {
    m_meshTheSurface = meshTheSurface;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setLevelSet(double ls)
  {
    m_ls = ls;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setRMC(double rmc)
  {
    m_rmc = rmc;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setBoundaryReference(
      const MaterialReference& ref)
  {
    m_isoref = ref;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::split(
      const MaterialReference& ref, const Split& s)
  {
    m_split[ref] = s;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::noSplit(
      const MaterialReference& ref)
  {
    m_split[ref] = NoSplit;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setSplit(const SplitMap& split)
  {
    assert(split.size() > 0);
    m_split = split;
    return *this;
  }

  int ImplicitDomainMesher::discretizeMMG2D(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      MMG2D_Set_dparameter(mesh, sol, MMG2D_DPARAM_rmc, *m_rmc);
    if (m_split.size() > 0)
    {
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_numberOfMat, m_split.size());
      for (const auto& v : m_split)
      {
        const auto& ref = v.first;
        const auto& split = v.second;

        std::visit(Utility::Overloaded{
          [&](const NoSplitT&)
          {
            if (!MMG2D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_NoSplit, ref, ref))
            {
              Alert::Exception() << "Could not set the multi-material reference lookup."
                                 << Alert::Raise;
            }
          },
          [&](const Split& s)
          {
            if (!MMG2D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_Split, s.interior, s.exterior))
            {
              Alert::Exception() << "Could not set the multi-material reference lookup."
                                 << Alert::Raise;
            }
          }
        }, split);
      }
    }

    if (m_isoref)
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_isoref, *m_isoref);
    if (m_meshTheSurface)
    {
      Alert::Exception()
        << "Meshing the surface for a 2D mesh is not supported."
        << Alert::Raise;
    }
    else
    {
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_iso, 1);
    }
    MMG2D_Set_dparameter(mesh, sol, MMG2D_DPARAM_ls, m_ls);
    return MMG2D_mmg2dls(mesh, sol, nullptr);
  }


  int ImplicitDomainMesher::discretizeMMG3D(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      MMG3D_Set_dparameter(mesh, sol, MMG3D_DPARAM_rmc, *m_rmc);
    if (m_isoref)
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_isoref, *m_isoref);
    if (m_split.size() > 0)
    {
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_numberOfMat, m_split.size());
      for (const auto& v : m_split)
      {
        const auto& ref = v.first;
        const auto& split = v.second;

        std::visit(Utility::Overloaded{
          [&](const NoSplitT&)
          {
            if (!MMG3D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_NoSplit, ref, ref))
            {
              Alert::Exception() << "Could not set the multi-material reference lookup."
                                 << Alert::Raise;
            }
          },
          [&](const Split& s)
          {
            if (!MMG3D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_Split, s.interior, s.exterior))
            {
              Alert::Exception() << "Could not set the multi-material reference lookup."
                                 << Alert::Raise;
            }
          }
        }, split);
      }
    }

    if (m_meshTheSurface)
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_isosurf, 1);
    else
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_iso, 1);

    MMG3D_Set_dparameter(mesh, sol, MMG3D_DPARAM_ls, m_ls);
    return MMG3D_mmg3dls(mesh, sol, nullptr);
  }

  int ImplicitDomainMesher::discretizeMMGS(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      Alert::Warning("Warning RMC option is not supported for surfaces").raise();
    if (m_isoref)
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_isoref, *m_isoref);
    if (m_split.size() > 0)
      Alert::Warning("Warning material splitting is not supported for surfaces").raise();


    if (m_meshTheSurface)
    {
      Alert::Exception()
        << "Meshing the surface for a surface mesh is not supported."
        << Alert::Raise;
    }
    else
    {
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_iso, 1);
    }
    MMGS_Set_dparameter(mesh, sol, MMGS_DPARAM_ls, m_ls);
    return MMGS_mmgsls(mesh, sol, nullptr);
  }
}
