/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "P1.h"

namespace Rodin::Variational
{
  const std::array<Math::Matrix, 3> P1Element::s_dofs =
  {
    Math::Matrix{{0}},
    Math::Matrix{{0}, {1}},
    Math::Matrix{{0, 1, 0},
                 {0, 0, 1}}
  };

  P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>>
  ::P1(const Geometry::Mesh<Context::Serial>& mesh)
    : m_mesh(mesh)
  {
    m_elements.resize(mesh.getDimension() + 1);
    for (size_t d = 1; d < mesh.getDimension() + 1; d++)
    {
      const size_t count = mesh.getConnectivity().getCount(d);
      for (size_t i = 0; i < count; i++)
        m_elements[d].emplace_back(mesh.getConnectivity().getGeometry(d, i));
    }
  }

  P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>>
  ::P1(const Geometry::Mesh<Context::Serial>& mesh, size_t vdim)
    : m_mesh(mesh), m_vdim(vdim)
  {
    m_elements.resize(mesh.getDimension() + 1);
    for (size_t d = 1; d < mesh.getDimension() + 1; d++)
    {
      const size_t count = mesh.getConnectivity().getCount(d);
      for (size_t i = 0; i < count; i++)
        m_elements[d].emplace_back(mesh.getConnectivity().getGeometry(d, i));
    }
  }
}
