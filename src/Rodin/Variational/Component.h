#ifndef RODIN_VARIATIONAL_COMPONENT_H
#define RODIN_VARIATIONAL_COMPONENT_H

#include "Rodin/Utility.h"

#include "ForwardDecls.h"
#include "GridFunction.h"
#include "TrialFunction.h"

namespace Rodin::Variational
{
  /**
   * @brief Represents the component (or entry) of a vectorial ShapeFunction.
   */
  template <class OperandDerived, class FES, ShapeFunctionSpaceType Space>
  class Component<ShapeFunction<OperandDerived, FES, Space>> final
  {
    public:
      using Operand = ShapeFunction<OperandDerived, FES, Space>;

      /**
       * @brief Constructs the component object from a TrialFunction and its
       * component index.
       */
      Component(Operand& u, size_t component)
        : m_u(u),
          m_idx(component)
      {}

      Component(const Component& other)
        : m_u(other.m_u),
          m_idx(other.m_idx)
      {}

      Component(Component&& other)
        : m_u(other.m_u),
          m_idx(other.m_idx)
      {}

      inline
      constexpr
      size_t getIndex() const
      {
        return m_idx;
      }

      inline
      constexpr
      Operand& getShapeFunction() const
      {
        return m_u.get();
      }

    private:
      std::reference_wrapper<Operand> m_u;
      const size_t m_idx;
  };

  template <class OperandDerived, class FES, ShapeFunctionSpaceType Space>
  Component(ShapeFunction<OperandDerived, FES, Space>&, size_t)
    -> Component<ShapeFunction<OperandDerived, FES, Space>>;

  /**
   * @brief Represents the component (or entry) of a vectorial FunctionBase
   * instance.
   */
  template <class OperandDerived>
  class Component<FunctionBase<OperandDerived>, size_t> final
    : public ScalarFunctionBase<Component<FunctionBase<OperandDerived>, size_t>>
  {
    public:
      using Operand = FunctionBase<OperandDerived>;
      using Parent = ScalarFunctionBase<Component<FunctionBase<OperandDerived>, size_t>>;
      using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;

      constexpr
      Component(const Operand& fn, size_t component)
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

      inline
      constexpr
      size_t getIndex() const
      {
        return m_idx;
      }

      inline
      constexpr
      const Operand& getOperand() const
      {
        assert(m_fn);
        return *m_fn;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).coeff(m_idx);
      }

      inline Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<Operand> m_fn;
      const size_t m_idx;
  };

  template <class OperandDerived>
  Component(const FunctionBase<OperandDerived>&, size_t) -> Component<FunctionBase<OperandDerived>, size_t>;

  template <class OperandDerived>
  class Component<FunctionBase<OperandDerived>, size_t, size_t> final
    : public ScalarFunctionBase<Component<FunctionBase<OperandDerived>, size_t, size_t>>
  {
    public:
      using Operand = FunctionBase<OperandDerived>;
      using Parent = ScalarFunctionBase<Component<FunctionBase<OperandDerived>, size_t, size_t>>;
      using OperandRange = typename FormLanguage::Traits<Operand>::RangeType;

      constexpr
      Component(const Operand& fn, size_t i, size_t j)
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

      inline
      constexpr
      const Operand& getOperand() const
      {
        assert(m_fn);
        return *m_fn;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return getOperand().getValue(p).coeff(m_i, m_j);
      }


      inline Component* copy() const noexcept override
      {
        return new Component(*this);
      }

    private:
      std::unique_ptr<Operand> m_fn;
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
    : public ScalarFunctionBase<Component<GridFunction<FES>>>
  {
    public:
      using Operand = GridFunction<FES>;
      using Parent = ScalarFunctionBase<Component<Operand>>;

      /**
       * @brief Constructs the component object from a GridFunction and its
       * component index.
       */
      constexpr
      Component(GridFunction<FES>& u, size_t component)
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

      inline
      constexpr
      GridFunction<FES>& getGridFunction()
      {
        return m_u.get();
      }

      inline
      constexpr
      const GridFunction<FES>& getGridFunction() const
      {
        return m_u.get();
      }

      inline
      Scalar getValue(const Geometry::Point& p) const
      {
        Math::Vector v;
        m_u.get().getValue(v, p);
        return v.coeff(m_idx);
      }

      inline Component* copy() const noexcept override
      {
        return new Component(*this);
      }

      // template <class NestedDerived,
      //          typename = std::enable_if_t<Utility::IsSpecialization<FES, H1>::Value>>
      // inline
      // constexpr
      // auto& projectOnBoundary(const FunctionBase<NestedDerived>& fn,
      //                         Geometry::Attribute attr)
      // {
      //   return projectOnBoundary(fn, FlatSet<Geometry::Attribute>{attr});
      // }

      // template <class NestedDerived,
      //          typename = std::enable_if_t<Utility::IsSpecialization<FES, H1>::Value>>
      // auto& projectOnBoundary(const FunctionBase<NestedDerived>& fn,
      //                         const FlatSet<Geometry::Attribute>& attrs = {})
      // {
      //   assert(false);
      //   return *this;
      // }

    private:
      std::reference_wrapper<GridFunction<FES>> m_u;
      const size_t m_idx;
  };

  template <class FES>
  Component(GridFunction<FES>&, size_t) -> Component<GridFunction<FES>>;
}

#endif

