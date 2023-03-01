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
  Mesh<Context::Serial>::Builder::setReference(Mesh<Context::Serial>& mesh)
  {
    m_dim = mesh.getDimension();
    m_sdim = mesh.getSpaceDimension();

    // Track the object
    m_ref.emplace(std::ref(mesh));

    // Set counts to zero
    m_count.resize(m_dim + 1, 0);

    // Emplace empty connectivity objects
    m_connectivity.resize(m_dim + 1);
    for (size_t i = 0; i < m_connectivity.size(); i++)
    {
      m_connectivity[i].reserve(m_connectivity.size());
      for (size_t j = 0; j < m_connectivity.size(); j++)
      {
        m_connectivity[i].push_back(Connectivity(i, j));
      }
    }

    // Emplace tranformation vectors
    m_transformations.resize(m_dim + 1);

    // Emplace the implementation
    m_impl.reset(new mfem::Mesh(m_dim, 0, 0, 0, m_sdim));

    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::vertex(const Math::Vector& x)
  {
    auto sdim = m_impl->SpaceDimension();
    assert(x.size() >= 0);
    if (static_cast<decltype(sdim)>(x.size()) != sdim)
    {
      Alert::Exception()
        << "Vertex dimension is different from space dimension"
        << " (" << x.size() << " != " << sdim << ")"
        << Alert::Raise;
    }
    m_impl->AddVertex(x.data());
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::element(Type geom, const Array<Index>& vs, Attribute attr)
  {
    assert(m_ref.has_value());
    const auto& ref = m_ref->get();
    mfem::Element* el = m_impl->NewElement(static_cast<int>(geom));
    std::copy_n(vs.begin(), el->GetNVertices(), el->GetVertices());
    el->SetAttribute(attr);
    const Index idx = m_impl->AddElement(el);
    m_connectivity[ref.getDimension()][0].connect(idx, vs);
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::face(Type geom, const Array<Index>& vs, Attribute attr)
  {
    mfem::Element* el = m_impl->NewElement(static_cast<int>(geom));
    std::copy_n(vs.begin(), el->GetNVertices(), el->GetVertices());
    el->SetAttribute(attr);
    m_impl->AddBdrElement(el);
    return *this;
  }

  void Mesh<Context::Serial>::Builder::finalize()
  {
    m_impl->FinalizeTopology();
    m_impl->Finalize(false, true);

    // TODO: Compute counts of all simplices
    m_count[m_dim] = m_impl->GetNE();
    m_count[m_dim - 1] = m_impl->GetNumFaces();
    m_count[0] = m_impl->GetNV();

    for (size_t d = 0; d < m_count.size(); d++)
      m_transformations[d].resize(m_count[d]);

    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    ref.m_impl = std::move(m_impl);
    ref.m_count = std::move(m_count);
    ref.m_connectivity = std::move(m_connectivity);
    ref.m_transformations = std::move(m_transformations);

    for (int i = 0; i < ref.getHandle().GetNBE(); i++)
      ref.m_f2b[ref.getHandle().GetBdrElementEdgeIndex(i)] = i;
  }
}
