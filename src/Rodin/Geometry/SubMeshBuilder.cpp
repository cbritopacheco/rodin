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
  SubMesh<Context::Serial>::Builder::include(size_t dim, std::set<Index> indices)
  {
    assert(m_mesh.has_value());
    assert(dim < m_indices.size());
    auto& mesh = m_mesh->get();
    assert(mesh.m_s2ps.size() == mesh.getDimension() + 1);
    auto& parent = m_mesh->get().getParent();
    if (dim == mesh.getDimension())
    {
      for (const auto& idx : indices)
      {
        auto simplex = parent.getSimplex(dim, idx);

        // Add simplex vertices to the resulting mesh
        mfem::Array<int> pvs;
        if (dim == parent.getDimension())
        {
          parent.getHandle().GetElementVertices(idx, pvs);
        }
        else if (dim == parent.getDimension() - 1)
        {
          parent.getHandle().GetFaceVertices(idx, pvs);
        }
        else
        {
          assert(false);
        }

        mfem::Array<int> sv(pvs.Size());
        for (int i = 0; i < sv.Size(); i++)
        {
          int pvid = pvs[i];
          if (mesh.m_s2ps[0].right.count(pvid) == 0) // Only add vertex if it is not in the map
          {
            sv[i] = mesh.getHandle().AddVertex(parent.getHandle().GetVertex(pvid));
            mesh.m_s2ps[0].insert({static_cast<size_t>(sv[i]), static_cast<size_t>(pvid)});
          }
          else // Else get the id of the vertex in the submesh
          {
            sv[i] = mesh.m_s2ps[0].right.at(pvid);
          }
        }

        // Add element with the new vertex ordering
        mfem::Element* newEl =
          mesh.getHandle().NewElement(static_cast<int>(simplex->getGeometry()));
        newEl->SetVertices(sv);
        newEl->SetAttribute(simplex->getAttribute());
        size_t seid = mesh.getHandle().AddElement(newEl);
        mesh.m_s2ps[dim].insert({seid, idx});
      }

      if (dim == parent.getDimension()) // we are not in the surface case
      {
        for (int i = 0; i < parent.getHandle().GetNBE(); i++)
        {
          Index faceIdx = parent.getHandle().GetBdrFace(i);
          Attribute attr = parent.getFaceAttribute(faceIdx);
          int el1 = -1, el2 = -1;
          parent.getHandle().GetFaceElements(faceIdx, &el1, &el2);
          if (indices.count(el1) || indices.count(el2))
          {
            assert(el1 >= 0 || el2 >= 0);
            mfem::Array<int> vs;
            parent.getHandle().GetFaceVertices(faceIdx, vs);
            for (int j = 0; j < vs.Size(); j++)
              vs[j] = mesh.m_s2ps[0].right.at(vs[j]);
            mfem::Element* newEl =
              mesh.getHandle().NewElement(parent.getHandle().GetFaceGeometry(faceIdx));
            newEl->SetVertices(vs);
            newEl->SetAttribute(attr);
            mesh.getHandle().AddBdrElement(newEl);
          }
        }
      }
    }
    else if (dim == mesh.getDimension() - 1)
    {
      assert(false);
    }
    else
    {
      assert(false);
    }
    m_indices[dim] = std::move(indices);
    return *this;
  }

  void SubMesh<Context::Serial>::Builder::finalize()
  {
    assert(m_mesh.has_value());
    auto& mesh = m_mesh->get();


    mesh.getHandle().FinalizeTopology();
    mesh.getHandle().Finalize(false, true);
    for (int i = 0; i < mesh.getHandle().GetNBE(); i++)
      mesh.m_f2b[mesh.getHandle().GetBdrElementEdgeIndex(i)] = i;
  }
}
