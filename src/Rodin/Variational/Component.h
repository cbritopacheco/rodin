#ifndef RODIN_VARIATIONAL_COMPONENT_H
#define RODIN_VARIATIONAL_COMPONENT_H

#include "Rodin/Utility.h"

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "TrialFunction.h"

namespace Rodin::FormLanguage
{
  template <class OperandDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Component<Variational::ShapeFunctionBase<OperandDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr const Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunctionBase<OperandDerived, FESType, SpaceType>;

    using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

    using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

    using RangeType = ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Represents the component (or entry) of a vectorial FunctionBase
   * instance.
   */
  template <class OperandDerived>
  class Component<FunctionBase<OperandDerived>, size_t> final
    : public RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t>>
  {
    public:
      using OperandType = FunctionBase<OperandDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t>>;


      constexpr
      Component(const OperandType& fn, size_t component)
        : m_fn(fn.copy()), m_idx(component)
      {}

      constexpr
      Component(const Component& other)
        : Parent(other),
          m_fn(other.m_fn->copy()),
          m_idx(other.m_idx)
      {}

      constexpr
      Component(Component&& other)
        : Parent(std::move(other)),
          m_fn(std::move(other.m_fn)),
          m_idx(std::move(other.m_idx))
      {}

      constexpr
      size_t getIndex() const
      {
        return m_idx;
      }

      constexpr
      const OperandType& getOperand() const
      {
        assert(m_fn);
        return *m_fn;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).coeff(m_idx);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<OperandType> m_fn;
      const size_t m_idx;
  };

  template <class OperandDerived>
  Component(const FunctionBase<OperandDerived>&, size_t) -> Component<FunctionBase<OperandDerived>, size_t>;

  template <class OperandDerived>
  class Component<FunctionBase<OperandDerived>, size_t, size_t> final
    : public RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t, size_t>>
  {
    public:
      using OperandType = FunctionBase<OperandDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Component<FunctionBase<OperandDerived>, size_t, size_t>>;

      constexpr
      Component(const OperandType& fn, size_t i, size_t j)
        : m_fn(fn.copy()), m_i(i), m_j(j)
      {}

      constexpr
      Component(const Component& other)
        : Parent(other),
          m_fn(other.m_fn->copy()),
          m_i(other.m_i),
          m_j(other.m_j)
      {}

      constexpr
      Component(Component&& other)
        : Parent(std::move(other)),
          m_fn(std::move(other.m_fn)),
          m_i(std::move(other.m_i)),
          m_j(std::move(other.m_j))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        assert(m_fn);
        return *m_fn;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).coeff(m_i, m_j);
      }


      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<OperandType> m_fn;
      const size_t m_i;
      const size_t m_j;
  };

  template <class OperandDerived>
  Component(const FunctionBase<OperandDerived>&, size_t, size_t)
    -> Component<FunctionBase<OperandDerived>, size_t, size_t>;

  /**
   * @brief Represents the component (or entry) of a vectorial GridFunction.
   */
  template <class FES>
  class Component<GridFunction<FES>> final
    : public RealFunctionBase<Component<GridFunction<FES>>>
  {
    public:
      using OperandType = GridFunction<FES>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Component<OperandType>>;

      /**
       * @brief Constructs the component object from a GridFunction and its
       * component index.
       */
      constexpr
      Component(OperandType& u, size_t component)
        : m_u(u), m_idx(component)
      {}

      constexpr
      Component(const Component& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      constexpr
      Component(Component&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      OperandType& getGridFunction()
      {
        return m_u.get();
      }

      constexpr
      const GridFunction<FES>& getGridFunction() const
      {
        return m_u.get();
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return m_u.get().getValue(p).coeff(m_idx);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::reference_wrapper<OperandType> m_u;
      const size_t m_idx;
  };

  template <class FES>
  Component(GridFunction<FES>&, size_t) -> Component<GridFunction<FES>>;

  /**
   * @brief Represents the component (or entry) of a vectorial ShapeFunction.
   */
  template <class OperandDerived, class FES, ShapeFunctionSpaceType Space>
  class Component<ShapeFunctionBase<OperandDerived, FES, Space>> final
    : public ShapeFunctionBase<Component<ShapeFunctionBase<OperandDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;
      static constexpr const ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunctionBase<OperandDerived, FESType, SpaceType>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandRangeType>::ScalarType;

      using RangeType = ScalarType;

      using Parent = ShapeFunctionBase<Component<OperandType>, FES, Space>;

      static_assert(std::is_same_v<OperandRangeType, Math::Vector<ScalarType>>);

      /**
       * @brief Constructs the component object from a TrialFunction and its
       * component index.
       */
      Component(const OperandType& u, size_t component)
        : Parent(u.getFiniteElementSpace()),
          m_u(u.copy()),
          m_idx(component)
      {}

      Component(const Component& other)
        : Parent(other),
          m_u(other.m_u->copy()),
          m_idx(other.m_idx)
      {}

      Component(Component&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_idx(std::move(other.m_idx))
      {}

      constexpr
      size_t getIndex() const
      {
        return m_idx;
      }

      constexpr
      const OperandType& getOperand() const
      {
        assert(m_u);
        return *m_u;
      }

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& polytope) const
      {
        const size_t d = polytope.getDimension();
        const size_t i = polytope.getIndex();
        return this->getFiniteElementSpace().getFiniteElement(d, i).getCount();
      }

      const Geometry::Point& getPoint() const
      {
        return m_u->getPoint();
      }

      Component& setPoint(const Geometry::Point& p)
      {
        m_u->setPoint(p);
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return getOperand().getBasis(local).coeff(m_idx);
      }

      Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<OperandType> m_u;
      const size_t m_idx;
  };

  template <class OperandDerived, class FES, ShapeFunctionSpaceType Space>
  Component(const ShapeFunctionBase<OperandDerived, FES, Space>&, size_t)
    -> Component<ShapeFunctionBase<OperandDerived, FES, Space>>;
}

#endif

