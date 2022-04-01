/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"
#include "Rodin/Utility/Overloaded.h"

#include "ImplicitDomainMesherS.h"

namespace Rodin::External::MMG
{
  ImplicitDomainMesherS::ImplicitDomainMesherS()
    : m_ls(0.0)
  {}

  ImplicitDomainMesherS& ImplicitDomainMesherS::setLevelSet(double ls)
  {
    m_ls = ls;
    return *this;
  }

  ImplicitDomainMesherS& ImplicitDomainMesherS::setBoundaryReference(
      const MaterialReference& ref)
  {
    m_isoref = ref;
    return *this;
  }

  ImplicitDomainMesherS::Discretization
  ImplicitDomainMesherS::discretize(ScalarSolutionS& ls)
  {
    MeshS mesh(ls.getMesh());
    auto sol = ScalarSolutionS(ls).setMesh(mesh);

    if (mesh.count(MeshS::Entity::Vertex) == 0)
    {
       Alert::Exception("Mesh vertex count is zero. Nothing to optimize.").raise();
    }
    else
    {
      if (m_hmin)
      {
        MMGS_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMGS_DPARAM_hmin, *m_hmin);
      }
      if (m_hmax)
      {
        MMGS_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMGS_DPARAM_hmax, *m_hmax);
      }
      if (m_hgrad)
      {
        MMGS_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMGS_DPARAM_hgrad, *m_hgrad);
      }
      if (m_hausd)
      {
        MMGS_Set_dparameter(
            mesh.getHandle(), sol.getHandle(), MMGS_DPARAM_hausd, *m_hausd);
      }
      if (m_isoref)
      {
        MMGS_Set_iparameter(mesh.getHandle(), sol.getHandle(),
            MMGS_IPARAM_isoref, *m_isoref);
      }
      MMGS_Set_iparameter(mesh.getHandle(), sol.getHandle(), MMGS_IPARAM_iso, 1);
      MMGS_Set_dparameter(mesh.getHandle(), sol.getHandle(), MMGS_DPARAM_ls, m_ls);
      MMGS_mmgsls(mesh.getHandle(), sol.getHandle(), NULL);
    }
    return {std::move(mesh), std::move(sol)};
  }

  ImplicitDomainMesherS& ImplicitDomainMesherS::setHMin(double hmin)
  {
    m_hmin = hmin;
    return *this;
  }

  ImplicitDomainMesherS& ImplicitDomainMesherS::setHMax(double hmax)
  {
    m_hmax = hmax;
    return *this;
  }

  ImplicitDomainMesherS& ImplicitDomainMesherS::setHausdorff(double hausd)
  {
    m_hausd = hausd;
    return *this;
  }

  ImplicitDomainMesherS& ImplicitDomainMesherS::setGradation(double hgrad)
  {
    m_hgrad = hgrad;
    return *this;
  }
}

