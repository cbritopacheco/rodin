#include "SubMesh.h"

#include "Element.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Serial>::SubMesh(const MeshBase& parent)
    : m_parent(parent)
  {}

  SubMesh<Context::Serial>::SubMesh(const SubMesh& other)
    : Mesh(other),
      m_parent(other.m_parent),
      m_s2ps(other.m_s2ps)
  {}

  SubMesh<Context::Serial>::SubMesh(SubMesh&& other)
    : Mesh(std::move(other)),
      m_parent(std::move(other.m_parent)),
      m_s2ps(std::move(other.m_s2ps))
  {}

  const MeshBase& SubMesh<Context::Serial>::getParent() const
  {
    return m_parent.get();
  }

  SubMesh<Context::Serial>::Builder
  SubMesh<Context::Serial>::initialize(size_t dim)
  {
    m_s2ps.resize(dim + 1);
    getHandle() = mfem::Mesh(dim, 0, 0, 0, getParent().getSpaceDimension());
    SubMesh<Context::Serial>::Builder build;
    build.setMesh(*this);
    return build;
  }
}
