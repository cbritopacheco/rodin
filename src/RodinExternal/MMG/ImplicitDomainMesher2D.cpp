/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"
#include "Rodin/Utility/Overloaded.h"

#include "ImplicitDomainMesher2D.h"

namespace Rodin::External::MMG
{
  ImplicitDomainMesher2D::ImplicitDomainMesher2D()
    : m_ls(0.0)
  {}

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setLevelSet(double ls)
  {
    m_ls = ls;
    return *this;
  }

  ImplicitDomainMesher2D::Discretization
  ImplicitDomainMesher2D::discretize(ScalarSolution2D& ls)
  {
    Mesh2D mesh(ls.getMesh());
    ScalarSolution2D sol(ls);

    if (mesh.count<Mesh2D::Vertex>() == 0)
    {
       Alert::Exception("Mesh vertex count is zero. Nothing to optimize.").raise();
    }
    else
    {
      if (m_hmin)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMG2D_DPARAM_hmin, *m_hmin);
      }
      if (m_hmax)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMG2D_DPARAM_hmax, *m_hmax);
      }
      if (m_hgrad)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMG2D_DPARAM_hgrad, *m_hgrad);
      }
      if (m_hausd)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMG2D_DPARAM_hausd, *m_hausd);
      }
      if (m_rmc)
      {
        MMG2D_Set_dparameter(mesh.getHandle(), sol.getHandle(),
            MMG2D_DPARAM_rmc, *m_rmc);
      }
      if (m_split.size() > 0)
      {
        MMG2D_Set_iparameter(mesh.getHandle(), sol.getHandle(),
            MMG2D_IPARAM_numberOfMat, m_split.size());
        for (const auto& v : m_split)
        {
          const auto& ref = v.first;
          const auto& split = v.second;

          std::visit(Utility::Overloaded{
            [&](const NoSplitT&)
            {
              if (!MMG2D_Set_multiMat(mesh.getHandle(), sol.getHandle(),
                    ref, MMG5_MMAT_NoSplit, ref, ref))
              {
                Alert::Exception()
                  << "Could not set the multi-material reference lookup"
                  << Alert::Raise;
              }
            },
            [&](const Split& s)
            {
              if (!MMG2D_Set_multiMat(mesh.getHandle(), sol.getHandle(),
                    ref, MMG5_MMAT_Split, s.interior, s.exterior))
              {
                Alert::Exception()
                  << "Could not set the multi-material reference lookup"
                  << Alert::Raise;
              }
            }
          }, split);
        }
      }
      MMG2D_Set_iparameter(mesh.getHandle(), sol.getHandle(),
          MMG2D_IPARAM_iso, 1);
      MMG2D_Set_dparameter(mesh.getHandle(), sol.getHandle(),
            MMG2D_DPARAM_ls, m_ls);
      MMG2D_mmg2dls(mesh.getHandle(), ls.getHandle(), NULL);
    }
    return {std::move(mesh), std::move(sol)};
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::split(
      const MaterialReference& ref, const Split& s)
  {
    m_split[ref] = s;
    return *this;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::noSplit(const MaterialReference& ref)
  {
    m_split[ref] = NoSplit;
    return *this;
  }

  const SplitMap& ImplicitDomainMesher2D::getSplit() const
  {
    return m_split;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setSplit(const SplitMap& split)
  {
    assert(split.size() > 0);
    m_split = split;
    return *this;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setRMC(double rmc)
  {
    assert(rmc > 0);
    m_rmc = rmc;
    return *this;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setHMin(double hmin)
  {
    m_hmin = hmin;
    return *this;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setHMax(double hmax)
  {
    m_hmax = hmax;
    return *this;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setHausdorff(double hausd)
  {
    m_hausd = hausd;
    return *this;
  }

  ImplicitDomainMesher2D& ImplicitDomainMesher2D::setGradation(double hgrad)
  {
    m_hgrad = hgrad;
    return *this;
  }
}
