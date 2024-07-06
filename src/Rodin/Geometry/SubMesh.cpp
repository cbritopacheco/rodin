#include "Rodin/Alert/MemberFunctionException.h"

#include "SubMesh.h"

#include "Polytope.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Sequential>::SubMesh(std::reference_wrapper<const Mesh<Context>> parent)
    : m_parent(parent)
  {
    if (parent.get().isSubMesh())
    {
      const auto& submesh = parent.get().asSubMesh();
      m_ancestors = submesh.getAncestors();
    }
    m_ancestors.push_front(parent.get());
  }

  SubMesh<Context::Sequential>::SubMesh(const SubMesh& other)
    : Parent(other),
      m_parent(other.m_parent),
      m_s2ps(other.m_s2ps),
      m_ancestors(other.m_ancestors)
  {}

  SubMesh<Context::Sequential>::SubMesh(SubMesh&& other)
    : Parent(std::move(other)),
      m_parent(std::move(other.m_parent)),
      m_s2ps(std::move(other.m_s2ps)),
      m_ancestors(std::move(other.m_ancestors))
  {}

  const Mesh<Context::Sequential>& SubMesh<Context::Sequential>::getParent() const
  {
    return m_parent.get();
  }

  std::optional<Point> SubMesh<Context::Sequential>::restriction(const Point& p) const
  {
    const auto& polytope = p.getPolytope();
    const size_t d = polytope.getDimension();
    Index i = polytope.getIndex();
    const auto& ancestors = getAncestors();
    Deque<std::reference_wrapper<const SubMeshBase>> descendants;
    descendants.push_front(*this);
    for (auto it = ancestors.begin(); it != ancestors.end(); ++it)
    {
      if (it->get() == polytope.getMesh())
      {
        break;
      }
      else if (!it->get().isSubMesh())
      {
        // Invalid restriction.
        // The SubMesh is not a descendant of the Mesh to which the Point belongs to.
        return {};
      }
      else
      {
        descendants.push_front(it->get().asSubMesh());
      }
    }
    for (auto it = descendants.begin(); it != descendants.end(); ++it)
    {
      const auto& polytopeMap = it->get().getPolytopeMap(d);
      auto find = polytopeMap.right.find(i);
      if (find == polytopeMap.right.end())
      {
        // Invalid restriction.
        // Could not find Polytope(d, i) in the SubMesh to parent Mesh map.
        return {};
      }
      i = find->get_left();
    }
    std::unique_ptr<Polytope> childPolytope(getPolytope(d, i).release());
    return Point(
        std::move(*childPolytope),
        getPolytopeTransformation(d, i),
        std::cref(p.getReferenceCoordinates()), p.getPhysicalCoordinates());
  }
}
