/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "SubMesh.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Serial>::Builder::Builder()
  {}

  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::setReference(
      Mesh<Context::Serial>::Builder&& build, SubMesh<Context::Serial>& mesh)
  {
    const size_t dim = mesh.getDimension();
    m_mbuild.emplace(std::move(build));
    m_ref.emplace(std::ref(mesh));
    m_s2ps.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::include(size_t dim, std::set<Index> indices)
  {
    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    assert(m_mbuild.has_value());
    auto& build = m_mbuild.value();

    assert(m_s2ps.size() == ref.getDimension() + 1);
    const auto& parent = m_ref->get().getParent();

    if (dim == ref.getDimension())
    {
      for (const auto& idx : indices)
      {
        auto simplex = parent.getSimplex(dim, idx);

        // Add simplex vertices to the resulting mesh
        const Array<Index>& pvs =
          parent.getConnectivity(dim, 0).getIncidence(idx);

        Array<Index> vs(pvs.size());
        for (Index i = 0; i < static_cast<Index>(vs.size()); i++)
        {
          Index pvidx = pvs[i];
          if (m_s2ps[0].right.count(pvidx) == 0) // Only add vertex if it is not in the map
          {
            build.vertex(parent.getVertex(pvidx)->coordinates());
            vs[i] = m_sidx[0]++;
            m_s2ps[0].insert({ vs[i], pvidx });
          }
          else // Else get the id of the vertex in the submesh
          {
            vs[i] = m_s2ps[0].right.at(pvidx);
          }
        }

        // Add element with the new vertex ordering
        build.element(simplex->getGeometry(), vs, simplex->getAttribute());
        m_s2ps[dim].insert({ m_sidx[dim]++, idx });
      }

      if (dim == parent.getDimension()) // We are not in the surface case
      {
        for (auto it = parent.getBoundary(); !it.end(); ++it)
        {
          int el1 = -1, el2 = -1;
          parent.getHandle().GetFaceElements(it->getIndex(), &el1, &el2);
          if (indices.count(el1) || indices.count(el2))
          {
            assert(el1 >= 0 || el2 >= 0);
            mfem::Array<int> pvs;
            parent.getHandle().GetFaceVertices(it->getIndex(), pvs);
            Array<Index> vs(pvs.Size());
            for (int i = 0; i < pvs.Size(); i++)
              vs[i] = m_s2ps[0].right.at(pvs[i]);
            build.face(it->getGeometry(), vs, it->getAttribute());
            m_sidx[ref.getDimension() - 1]++;
          }
        }
      }
    }
    else if (dim == ref.getDimension() - 1)
    {
      assert(false);
    }
    else
    {
      assert(false);
    }
    return *this;
  }

  void SubMesh<Context::Serial>::Builder::finalize()
  {
    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    assert(m_mbuild.has_value());
    auto& build = m_mbuild.value();

    build.finalize();

    ref.m_s2ps = std::move(m_s2ps);
  }
}
