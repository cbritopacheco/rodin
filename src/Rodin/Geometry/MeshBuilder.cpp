/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  Mesh<Context::Serial>::Builder::Builder()
  {}

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::setMesh(Mesh<Context::Serial>& mesh)
  {
    m_mesh = mesh;
    m_connectivity.emplace(
        std::pair(mesh.getDimension(), 0), Connectivity(mesh.getDimension(), 0));
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::vertex(const std::vector<double>& x)
  {
    assert(m_mesh.has_value());
    auto& mesh = m_mesh->get();
    auto sdim = mesh.getSpaceDimension();
    if (x.size() != sdim)
    {
      Alert::Exception()
        << "Vertex dimension is different from space dimension"
        << " (" << x.size() << " != " << sdim << ")"
        << Alert::Raise;
    }
    mesh.getHandle().AddVertex(x.data());
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::element(
      Type geom,
      const std::vector<Index>& vs, Attribute attr)
  {
    assert(m_mesh.has_value());
    auto& mesh = m_mesh->get();
    mfem::Element* el = mesh.getHandle().NewElement(static_cast<int>(geom));
    for (int i = 0; i < el->GetNVertices(); i++)
      el->GetVertices()[i] = vs[i];
    el->SetAttribute(attr);
    mesh.getHandle().AddElement(el);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::face(
      Type geom,
      const std::vector<Index>& vs, Attribute attr)
  {
    assert(m_mesh.has_value());
    auto& mesh = m_mesh->get();
    mfem::Element* el = mesh.getHandle().NewElement(static_cast<int>(geom));
    for (int i = 0; i < el->GetNVertices(); i++)
      el->GetVertices()[i] = vs[i];
    el->SetAttribute(attr);
    mesh.getHandle().AddBdrElement(el);
    return *this;
  }

  void Mesh<Context::Serial>::Builder::finalize()
  {
    assert(m_mesh.has_value());
    auto& mesh = m_mesh->get();

    mesh.m_connectivity = std::move(m_connectivity);

    mesh.getHandle().FinalizeTopology();
    mesh.getHandle().Finalize(false, true);

    std::vector<Index> indices;
    std::vector<size_t> offsets(mesh.getHandle().GetNE() + 1);
    offsets[0] = 0;

    for (int i = 0; i < mesh.getHandle().GetNE(); i++)
    {
      auto offset = mesh.getHandle().GetElement(i)->GetNVertices();
      for (int j = 0; j < offset; j++)
        indices.push_back(mesh.getHandle().GetElement(i)->GetVertices()[j]);
      offsets[i + 1] = offsets[i] + offset;
    }

    mesh.m_connectivity.clear();
    auto [it, inserted] = mesh.m_connectivity.emplace(
        std::pair(mesh.getDimension(), 0),
        Connectivity(mesh.getDimension(), 0));
    assert(inserted);

    it->second.setOffsets(offsets);
    it->second.setOffsets(indices);

    for (int i = 0; i < mesh.getHandle().GetNBE(); i++)
      mesh.m_f2b[mesh.getHandle().GetBdrElementEdgeIndex(i)] = i;
  }
}
