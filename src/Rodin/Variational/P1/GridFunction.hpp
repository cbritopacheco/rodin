#ifndef RODIN_VARIATIONAL_P1_GRIDFUNCTION_HPP
#define RODIN_VARIATIONAL_P1_GRIDFUNCTION_HPP

#include "GridFunction.h"

namespace Rodin::Variational
{
  template <class Range, class Context, class Mesh>
  auto GridFunction<P1<Range, Context, Mesh>>
  ::getValue(const Geometry::Point& p) const
  {
    const auto& fes = this->getFiniteElementSpace();
    const auto& fesMesh = fes.getMesh();
    const auto& polytope = p.getPolytope();
    assert(fesMesh == polytope.getMesh());
    const size_t d = polytope.getDimension();
    const Index i = polytope.getIndex();
    const auto& fe = fes.getFiniteElement(d, i);
    const auto& r = p.getCoordinates(Geometry::Point::Coordinates::Reference);
    if constexpr (std::is_same_v<RangeType, Scalar>)
    {
      Scalar res = 0;
      for (Index local = 0; local < fe.getCount(); local++)
      {
        const auto& basis = fe.getBasis(local);
        res += getValue({d, i}, local) * basis(r);
      }
      return res;
    }
    else if constexpr (std::is_same_v<RangeType, Math::Vector>)
    {
      Math::Vector res;
      getValueByReference(res, p);
      return res;
    }
    else
    {
      assert(false);
      return void();
    }
  }

  template <class Range, class Context, class Mesh>
  void GridFunction<P1<Range, Context, Mesh>>
  ::getValueByReference(Math::Vector& res, const Geometry::Point& p) const
  {
    static_assert(std::is_same_v<RangeType, Math::Vector>);
    const auto& fes = this->getFiniteElementSpace();
    const auto& polytope = p.getPolytope();
    const size_t d = polytope.getDimension();
    const Index i = polytope.getIndex();
    const auto& fe = fes.getFiniteElement(d, i);
    const auto& r = p.getCoordinates(Geometry::Point::Coordinates::Reference);
    const size_t vdim = fes.getVectorDimension();
    const size_t dofs = fe.getCount();
    res.resize(vdim);
    res.setZero();
    for (Index local = 0; local < dofs; local++)
    {
      const auto& basis = fe.getBasis(local);
      res += getValue({d, i}, local).coeff(local % vdim) * basis(r);
    }
  }
}

#endif
