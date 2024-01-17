/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_AVERAGE_H
#define RODIN_VARIATIONAL_AVERAGE_H

#include "ForwardDecls.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
  template <class Derived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Average<ShapeFunctionBase<Derived, FESType, SpaceType>> final
    : public ShapeFunctionBase<Average<ShapeFunctionBase<Derived, FESType, SpaceType>>>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using Operand = ShapeFunctionBase<Derived, FESType, SpaceType>;
      using Parent = ShapeFunctionBase<Average<ShapeFunctionBase<Derived, FESType, SpaceType>>>;
      using Range = typename FormLanguage::Traits<Operand>::RangeType;

      constexpr
      Average(const Operand& op)
        : Parent(op.getFiniteElementSpace()),
          m_op(op.copy())
      {}

      constexpr
      Average(const Average& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      constexpr
      Average(Average&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      inline
      constexpr
      const Operand& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape();
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      inline
      auto getTensorBasis(const Geometry::Point& p) const
      {
        assert(p.getPolytope().isFace());
        const auto& face = p.getPolytope();
        const size_t d = face.getDimension();
        const auto& mesh = face.getMesh();
        const auto& inc = mesh.getConnectivity().getIncidence({ d, d + 1 }, face.getIndex() );
        assert(inc.size() == 2);
        const Index idx1 = *inc.begin();
        const Index idx2 = *std::next(inc.begin());
        const auto it1 = mesh.getPolytope(d + 1, idx1);
        const auto it2 = mesh.getPolytope(d + 1, idx2);
        const auto& pc = p.getPhysicalCoordinates();
        const Math::SpatialVector rc1 = it1->getTransformation().inverse(pc);
        const Math::SpatialVector rc2 = it2->getTransformation().inverse(pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
        const auto& lhs = this->object(getOperand().getTensorBasis(p1));
        const auto& rhs = this->object(getOperand().getTensorBasis(p2));
        return 0.5 * (lhs + rhs);
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      inline Average* copy() const noexcept override
      {
        return new Average(*this);
      }

    private:
      std::unique_ptr<Operand> m_op;
  };

  template <class Derived, class FESType, ShapeFunctionSpaceType SpaceType>
  Average(const ShapeFunctionBase<Derived, FESType, SpaceType>&)
    -> Average<ShapeFunctionBase<Derived, FESType, SpaceType>>;
}

#endif
