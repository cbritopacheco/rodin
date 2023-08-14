#include "SubMesh.h"

#include "Simplex.h"

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
  SubMesh<Context::Serial>::initialize(size_t dim, size_t sdim)
  {
    assert(sdim == getParent().getSpaceDimension());
    SubMesh<Context::Serial>::Builder build;
    auto mbuild = Mesh<Context::Serial>::initialize(dim, getParent().getSpaceDimension());
    build.setReference(std::move(mbuild), *this);
    return build;
  }
}
