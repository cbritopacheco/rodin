/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MeshOptimizer2D.h"

namespace Rodin::External::MMG
{
  MeshOptimizer2D& MeshOptimizer2D::setHMin(double hmin)
  {
    m_hmin = hmin;
    return *this;
  }

  MeshOptimizer2D& MeshOptimizer2D::setHMax(double hmax)
  {
    m_hmax = hmax;
    return *this;
  }

  MeshOptimizer2D& MeshOptimizer2D::setGradation(double gradation)
  {
    m_hgrad = gradation;
    return *this;
  }

  MeshOptimizer2D& MeshOptimizer2D::setHausdorff(double hausd)
  {
    m_hausd = hausd;
    return *this;
  }

  ScalarSolution2D MeshOptimizer2D::optimize(Mesh2D& mesh)
  {
    ScalarSolution2D metric(mesh);
    if (mesh.count<Mesh2D::Vertex>() == 0)
    {
       Alert::Exception("Mesh vertex count is zero. Nothing to optimize.").raise();
    }
    else
    {
      if (m_hmin)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), metric.getHandle(), MMG2D_DPARAM_hmin, *m_hmin);
      }
      if (m_hmax)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), metric.getHandle(), MMG2D_DPARAM_hmax, *m_hmax);
      }
      if (m_hgrad)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), metric.getHandle(), MMG2D_DPARAM_hgrad, *m_hgrad);
      }
      if (m_hausd)
      {
        MMG2D_Set_dparameter(
            mesh.getHandle(), metric.getHandle(), MMG2D_DPARAM_hausd, *m_hausd);
      }
      MMG2D_Set_iparameter(mesh.getHandle(), metric.getHandle(),
          MMG2D_IPARAM_optim, 1);
      MMG2D_mmg2dlib(mesh.getHandle(), metric.getHandle());
    }
    return metric;
  }
}

