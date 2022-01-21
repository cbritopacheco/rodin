/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert/Exception.h"
#include "Rodin/Utility/Overloaded.h"

#include "MeshAdaptor2D.h"

namespace Rodin::External::MMG
{
  MeshAdaptor2D& MeshAdaptor2D::setHMin(double hmin)
  {
    m_hmin = hmin;
    return *this;
  }

  MeshAdaptor2D& MeshAdaptor2D::setHMax(double hmax)
  {
    m_hmax = hmax;
    return *this;
  }

  MeshAdaptor2D& MeshAdaptor2D::setGradation(double gradation)
  {
    m_hgrad = gradation;
    return *this;
  }

  MeshAdaptor2D& MeshAdaptor2D::setHausdorff(double hausd)
  {
    m_hausd = hausd;
    return *this;
  }

  MeshAdaptor2D& MeshAdaptor2D::setMetric(const ScalarSolution2D& metric)
  {
    m_metric = metric;
    return *this;
  }

  void MeshAdaptor2D::adapt(Mesh2D& mesh) const
  {
    if (mesh.count<Mesh2D::Vertex>() == 0)
    {
       Alert::Exception("Mesh vertex count is zero. Nothing to optimize.").raise();
    }
    else
    {
      std::visit(Utility::Overloaded
      {
        [&mesh, this](const auto& metric)
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
        MMG2D_mmg2dlib(mesh.getHandle(), metric.getHandle());
        }
      }, m_metric ? *m_metric : ScalarSolution2D<>(mesh));
    }
  }
}
