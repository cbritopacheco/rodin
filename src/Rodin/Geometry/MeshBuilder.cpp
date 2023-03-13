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
    m_count.initialize(m_dim);

    // Set indexes
    m_attrs.initialize(m_dim);
    m_transformations.initialize(m_dim);

    // Emplace empty connectivity objects
    m_connectivity.initialize(m_dim);

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
    Index idx = m_impl->AddVertex(x.data());
    m_connectivity.connect({0, 0}, {idx});
    return *this;
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::element(Type geom, const Array<Index>& vs, Attribute attr)
  {
    assert(getGeometryDimension(geom) == m_ref->get().getDimension());
    return simplex(geom, vs, attr);
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::face(Type geom, const Array<Index>& vs, Attribute attr)
  {
    assert(getGeometryDimension(geom) == m_ref->get().getDimension() - 1);
    return simplex(geom, vs, attr);
  }

  Mesh<Context::Serial>::Builder&
  Mesh<Context::Serial>::Builder::simplex(Type geom, const Array<Index>& vs, Attribute attr)
  {
    assert(m_ref.has_value());
    const auto& ref = m_ref->get();

    assert(getGeometryDimension(geom) <= ref.getDimension());
    const size_t dim = getGeometryDimension(geom);

    if (dim == ref.getDimension())
    {
      assert(m_ref.has_value());
      const auto& ref = m_ref->get();
      mfem::Element* el = m_impl->NewElement(static_cast<int>(geom));
      std::copy_n(vs.begin(), el->GetNVertices(), el->GetVertices());
      el->SetAttribute(attr);
      m_impl->AddElement(el);
      m_connectivity.connect({ref.getDimension(), 0}, vs);
    }
    else if (dim == ref.getDimension() - 1)
    {
      mfem::Element* el = m_impl->NewElement(static_cast<int>(geom));
      std::copy_n(vs.begin(), el->GetNVertices(), el->GetVertices());
      el->SetAttribute(attr);
      m_impl->AddBdrElement(el);
    }
    else
    {
      assert(false);
    }
    return *this;
  }

  void Mesh<Context::Serial>::Builder::finalize()
  {
    m_impl->FinalizeTopology();
    m_impl->Finalize(false, true);

    // TODO: Compute counts of all simplices
    m_count.at(m_dim) = m_impl->GetNE();
    m_count.at(m_dim - 1) = m_impl->GetNumFaces();
    m_count.at(0) = m_impl->GetNV();

    // TODO: Verify that tracking these attributes is correct
    m_attrs.reserve(m_dim, m_count.at(m_dim));
    for (size_t i = 0; i < m_count.at(m_dim); i++)
      m_attrs.track(m_dim, i, m_impl->GetAttribute(i));

    m_attrs.reserve(m_dim - 1, m_impl->GetNBE());
    for (int i = 0; i < m_impl->GetNBE(); i++)
      m_attrs.track(m_dim - 1, m_impl->GetBdrFace(i), m_impl->GetBdrAttribute(i));

    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    ref.m_impl = std::move(m_impl);

    ref.m_count = std::move(m_count);
    ref.m_connectivity = std::move(m_connectivity);

    ref.m_attrs = std::move(m_attrs);
    ref.m_transformations = std::move(m_transformations);
  }
}
