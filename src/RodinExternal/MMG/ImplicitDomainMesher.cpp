#include "ImplicitDomainMesher.h"

namespace Rodin::External::MMG
{
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

  void ImplicitDomainMesher::discretizeMMG2D(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      MMG2D_Set_dparameter(mesh, sol, MMG2D_DPARAM_rmc, *m_rmc);
    if (m_isoref)
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_isoref, *m_isoref);
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
    MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_iso, 1);
    MMG2D_Set_dparameter(mesh, sol, MMG2D_DPARAM_ls, m_ls);
    MMG2D_mmg2dls(mesh, sol, nullptr);
  }


  void ImplicitDomainMesher::discretizeMMG3D(MMG5_pMesh mesh, MMG5_pSol sol)
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
    MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_iso, 1);
    MMG3D_Set_dparameter(mesh, sol, MMG3D_DPARAM_ls, m_ls);
    MMG3D_mmg3dls(mesh, sol, nullptr);
  }

  void ImplicitDomainMesher::discretizeMMGS(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      Alert::Warning("Warning RMC option is not supported for surfaces").raise();
    if (m_isoref)
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_isoref, *m_isoref);
    if (m_split.size() > 0)
      Alert::Warning("Warning splitting option is not supported for surfaces").raise();
    MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_iso, 1);
    MMGS_Set_dparameter(mesh, sol, MMGS_DPARAM_ls, m_ls);
    MMGS_mmgsls(mesh, sol, nullptr);
  }
}
