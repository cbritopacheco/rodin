#include "Rodin/Alert/MemberFunctionException.h"

#include "SubMesh.h"

#include "Polytope.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Serial>::SubMesh(std::reference_wrapper<const Mesh<Context::Serial>> parent)
    : m_parent(parent)
  {}

  SubMesh<Context::Serial>::SubMesh(const SubMesh& other)
    : Parent(other),
      m_parent(other.m_parent),
      m_s2ps(other.m_s2ps)
  {}

  SubMesh<Context::Serial>::SubMesh(SubMesh&& other)
    : Parent(std::move(other)),
      m_parent(std::move(other.m_parent)),
      m_s2ps(std::move(other.m_s2ps))
  {}

  const Mesh<Context::Serial>& SubMesh<Context::Serial>::getParent() const
  {
    return m_parent.get();
  }

  Point SubMesh<Context::Serial>::inclusion(const Point& p) const
  {
    const auto& parent = getParent();
    const auto& polytope = p.getPolytope();
    assert(*this == polytope.getMesh());
    const size_t d = polytope.getDimension();
    const Index i = polytope.getIndex();
    const Index pi = getPolytopeMap(d).left.at(i);
    return Point(std::move(*parent.getPolytope(d, pi)), parent.getPolytopeTransformation(d, pi),
        std::cref(p.getReferenceCoordinates()), p.getPhysicalCoordinates());
  }

  Point SubMesh<Context::Serial>::restriction(const Point& p) const
  {
    const auto& child = *this;
    const auto& polytope = p.getPolytope();
    assert(getParent() == polytope.getMesh());
    const size_t d = polytope.getDimension();
    const Index i = polytope.getIndex();
    const auto& pmap = getPolytopeMap(d);
    const auto find = pmap.right.find(i);
    if (find == pmap.right.end())
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Invalid restriction. Could not find " << Alert::Notation::Polytope(d, i)
        << " in the SubMesh to parent Mesh map."
        << Alert::Raise;
    }
    const Index ci = find->get_left();
    return Point(std::move(*child.getPolytope(d, ci)), child.getPolytopeTransformation(d, ci),
        std::cref(p.getReferenceCoordinates()), p.getPhysicalCoordinates());
  }
}
