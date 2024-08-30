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
  template <class FunctionDerived>
  class Average<FunctionBase<FunctionDerived>> final
    : public FunctionBase<Average<FunctionBase<FunctionDerived>>>
  {
    public:
      using OperandType = FunctionBase<FunctionDerived>;

      using Parent = FunctionBase<Average<FunctionBase<FunctionDerived>>>;

      constexpr
      Average(const OperandType& op)
        : m_operand(op.copy())
      {}

      constexpr
      Average(const Average& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Average(Average&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      RangeShape getRangeShape() const
      {
        return getOperand().getRangeShape();
      }

      constexpr
      const auto& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      constexpr
      Average& traceOf(Geometry::Attribute attr)
      {
        return *this;
      }

      constexpr
      Average& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        return *this;
      }

      template <class T>
      void getValue(T& res, const Geometry::Point& p) const
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
        const Math::SpatialVector<Real> rc1 = it1->getTransformation().inverse(pc);
        const Math::SpatialVector<Real> rc2 = it2->getTransformation().inverse(pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
        getOperand().getValue(res, p1);
        res += getOperand().getValue(p2);
        res *= 0.5;
      }

      auto getValue(const Geometry::Point& p) const
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
        const Math::SpatialVector<Real> rc1 = it1->getTransformation().inverse(pc);
        const Math::SpatialVector<Real> rc2 = it2->getTransformation().inverse(pc);
        const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
        const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
        return 0.5 * (getOperand().getValue(p1) + getOperand().getValue(p2));
      }

      Average* copy() const noexcept override
      {
        return new Average(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class Derived>
  Average(const FunctionBase<Derived>&) -> Average<FunctionBase<Derived>>;

  // template <class Derived, class FES, ShapeFunctionSpaceType Space>
  // class Average<ShapeFunctionBase<Derived, FES, Space>> final
  //   : public ShapeFunctionBase<Average<ShapeFunctionBase<Derived, FES, Space>>>
  // {
  //   public:
  //     using FESType = FES;
  //     static constexpr ShapeFunctionSpaceType SpaceType = Space;

  //     using OperandType = ShapeFunctionBase<Derived, FESType, SpaceType>;

  //     using RangeType = typename FormLanguage::Traits<OperandType>::RangeType;

  //     using Parent = ShapeFunctionBase<Average<ShapeFunctionBase<Derived, FESType, SpaceType>>>;

  //     constexpr
  //     Average(const OperandType& op)
  //       : Parent(op.getFiniteElementSpace()),
  //         m_operand(op.copy())
  //     {}

  //     constexpr
  //     Average(const Average& other)
  //       : Parent(other),
  //         m_operand(other.m_operand->copy())
  //     {}

  //     constexpr
  //     Average(Average&& other)
  //       : Parent(std::move(other)),
  //         m_operand(std::move(other.m_operand))
  //     {}

  //     inline
  //     constexpr
  //     const OperandType& getOperand() const
  //     {
  //       assert(m_operand);
  //       return *m_operand;
  //     }

  //     inline
  //     constexpr
  //     const auto& getLeaf() const
  //     {
  //       return getOperand().getLeaf();
  //     }

  //     inline
  //     constexpr
  //     RangeShape getRangeShape() const
  //     {
  //       return getOperand().getRangeShape();
  //     }

  //     inline
  //     constexpr
  //     size_t getDOFs(const Geometry::Polytope& element) const
  //     {
  //       return getOperand().getDOFs(element);
  //     }

  //     inline
  //     const Geometry::Point& getPoint() const
  //     {
  //       return m_operand->getPoint();
  //     }

  //     inline
  //     Average& setPoint(const Geometry::Point& p)
  //     {
  //       m_operand->setPoint(p);
  //       return *this;
  //     }

  //     inline
  //     auto getBasis(size_t local, const Geometry::Point& p) const
  //     {
  //       assert(p.getPolytope().isFace());
  //       const auto& face = p.getPolytope();
  //       const size_t d = face.getDimension();
  //       const auto& mesh = face.getMesh();
  //       const auto& inc = mesh.getConnectivity().getIncidence({ d, d + 1 }, face.getIndex() );
  //       assert(inc.size() == 2);
  //       const Index idx1 = *inc.begin();
  //       const Index idx2 = *std::next(inc.begin());
  //       const auto it1 = mesh.getPolytope(d + 1, idx1);
  //       const auto it2 = mesh.getPolytope(d + 1, idx2);
  //       const auto& pc = p.getPhysicalCoordinates();
  //       const Math::SpatialVector<Real> rc1 = it1->getTransformation().inverse(pc);
  //       const Math::SpatialVector<Real> rc2 = it2->getTransformation().inverse(pc);
  //       const Geometry::Point p1(std::cref(*it1), std::cref(rc1), pc);
  //       const Geometry::Point p2(std::cref(*it2), std::cref(rc2), pc);
  //       const auto& lhs = this->object(getOperand().getBasis(local, p1));
  //       const auto& rhs = this->object(getOperand().getBasis(local, p2));
  //       return 0.5 * (lhs + rhs);
  //     }

  //     inline
  //     constexpr
  //     const auto& getFiniteElementSpace() const
  //     {
  //       return getOperand().getFiniteElementSpace();
  //     }

  //     inline Average* copy() const noexcept override
  //     {
  //       return new Average(*this);
  //     }

  //   private:
  //     std::unique_ptr<OperandType> m_operand;
  // };

  // template <class Derived, class FESType, ShapeFunctionSpaceType SpaceType>
  // Average(const ShapeFunctionBase<Derived, FESType, SpaceType>&)
  //   -> Average<ShapeFunctionBase<Derived, FESType, SpaceType>>;
}

#endif
