/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
      MMG2D_Set_iparameter(mesh.getHandle(), sol.getHandle(),
          MMG2D_IPARAM_iso, 1);
      MMG2D_Set_dparameter(mesh.getHandle(), sol.getHandle(),
            MMG2D_DPARAM_ls, m_ls);
      MMG2D_mmg2dls(mesh.getHandle(), ls.getHandle(), NULL);
    }
    return {std::move(mesh), std::move(sol)};
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
