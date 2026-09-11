/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file UnaryMinus.h
 * @brief Negation operator for functions, shape functions, and integrators.
 *
 * Provides the unary minus operator (-) for negating functions, shape functions,
 * and form integrators. For a function @f$ f @f$, the negation is:
 * @f[
 *   (-f)(x) = -f(x)
 * @f]
 */

#ifndef RODIN_VARIATIONAL_UNARYMINUS_H
#define RODIN_VARIATIONAL_UNARYMINUS_H

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/FormLanguage/List.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "ShapeFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c UnaryMinus over a shape function: exposes the finite
  /// element space and the shape function space.
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::UnaryMinus<Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup UnaryMinusSpecializations UnaryMinus Template Specializations
   * @brief Template specializations of the UnaryMinus class.
   * @see @ref UnaryMinus
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref UnaryMinus "UnaryMinus<FunctionBase<NestedDerived>>" | Pointwise negation of a function expression. |
   * | @ref UnaryMinus "UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>" | Negation of a test or trial shape-function expression. |
   * | @ref UnaryMinus "UnaryMinus<LinearFormIntegratorBase<Number>>" | Negation of a linear-form integrator. |
   * | @ref UnaryMinus "UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>" | Negation of a list of linear-form integrators. |
   * | @ref UnaryMinus "UnaryMinus<LocalBilinearFormIntegratorBase<Number>>" | Negation of a local bilinear-form integrator. |
   * | @ref UnaryMinus "UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>>" | Negation of a list of local bilinear-form integrators. |
   */

  /**
   * @ingroup UnaryMinusSpecializations
   * @brief Negation of a function.
   *
   * Represents the negation:
   * @f[
   *   (-f)(x) = -f(x)
   * @f]
   */
  template <class NestedDerived>
  class UnaryMinus<FunctionBase<NestedDerived>> final
    : public FunctionBase<UnaryMinus<FunctionBase<NestedDerived>>>
  {
    public:
      /// @brief Operand type.
      using OperandType = FunctionBase<NestedDerived>;

      /// @brief Range type of the operand.
      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      /// @brief Parent class type.
      using Parent = FunctionBase<UnaryMinus<OperandType>>;

      constexpr
      /// @brief Constructs the expression from its operand.
      UnaryMinus(const OperandType& op)
        : m_op(op.copy())
      {}

      constexpr
      /// @brief Copy constructor.
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      constexpr
      /// @brief Move constructor.
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      /**
       * @brief Evaluates the negated function at a point.
       * @param p Point at which to evaluate
       * @returns Negated function value @f$ -f(p) @f$
       */
      template <class Point>
      constexpr
      auto getValue(const Point& p) const
      {
        const auto v = getOperand().getValue(p);
        if constexpr (requires { (-v).eval(); })
          return (-v).eval();
        else
          return -v;
      }

      constexpr
      /// @brief Restricts the trace of the expression to a mesh attribute.
      UnaryMinus& traceOf(Geometry::Attribute attr)
      {
        Parent::traceOf(attr);
        m_op->traceOf(attr);
        return *this;
      }

      constexpr
      /// @brief Restricts the trace of the expression to a mesh attribute.
      UnaryMinus& traceOf(const FlatSet<Geometry::Attribute>& attrs)
      {
        Parent::traceOf(attrs);
        m_op->traceOf(attrs);
        return *this;
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return getOperand().getOrder(polytope);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<OperandType> m_op;
  };

  /**
   * @brief Deduction guide for UnaryMinus of FunctionBase.
   */
  template <class NestedDerived>
  UnaryMinus(const FunctionBase<NestedDerived>&) -> UnaryMinus<FunctionBase<NestedDerived>>;

  /**
   * @brief Applies unary minus to a function.
   * @param op Function to negate
   * @returns Negated function
   */
  template <class NestedDerived>
  constexpr
  auto operator-(const FunctionBase<NestedDerived>& op)
  {
    return UnaryMinus(op);
  }

  /**
   * @ingroup UnaryMinusSpecializations
   * @brief Negation of a shape function.
   *
   * Represents the negation:
   * @f[
   *   (-u)_i = -u_i
   * @f]
   * where @f$ u_i @f$ are the basis functions.
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Operand type.
      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>;

      constexpr
      /// @brief Constructs the expression from its operand.
      UnaryMinus(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy())
      {}

      constexpr
      /// @brief Copy constructor.
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      /// @brief Move constructor.
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        return *m_operand;
      }

      constexpr
      /// @brief Gets the operand in the shape function expression.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Gets the number of degrees of freedom on an element.
       * @param element Polytope element
       * @returns Number of DOFs
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /// @brief Sets the integration point the expression is evaluated at.
      UnaryMinus& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_operand->setIntegrationPoint(ip);
        return *this;
      }

      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_operand->getIntegrationPoint();
      }

      /**
       * @brief Evaluates the negated basis function.
       * @param local Local DOF index
       * @returns Negated basis function value @f$ -u_i @f$
       */
      constexpr
      auto getBasis(size_t local) const
      {
        decltype(auto) v = getOperand().getBasis(local);
        if constexpr (requires { (-v).eval(); })
          return (-v).eval();
        else
          return -v;
      }

      /// @brief Gets the finite element space.
      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return getOperand().getOrder(polytope);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for UnaryMinus of ShapeFunctionBase.
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  UnaryMinus(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> UnaryMinus<ShapeFunctionBase<NestedDerived, FES, Space>>;


  /**
   * @brief Applies unary minus to a shape function.
   * @param op Shape function to negate
   * @returns Negated shape function
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  constexpr
  auto operator-(const ShapeFunctionBase<NestedDerived, FES, Space>& op)
  {
    return UnaryMinus(op);
  }

  /**
   * @ingroup UnaryMinusSpecializations
   * @brief Negation of a linear form integrator.
   */
  template <class Number>
  class UnaryMinus<LinearFormIntegratorBase<Number>> final
    : public LinearFormIntegratorBase<Number>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Number;

      /// @brief Operand type.
      using OperandType = LinearFormIntegratorBase<ScalarType>;

      /// @brief Parent class type.
      using Parent = LinearFormIntegratorBase<ScalarType>;

      /// @brief Constructs the expression from its operand.
      UnaryMinus(const OperandType& op)
        : Parent(op),
          m_op(op.copy())
      {}

      /// @brief Copy constructor.
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      /// @brief Move constructor.
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      /// @brief Returns the integration region.
      Geometry::Region getRegion() const override
      {
        return getOperand().getRegion();
      }

      /// @brief Returns the polytope the expression is bound to.
      const Geometry::Polytope& getPolytope() const override
      {
        return getOperand().getPolytope();
      }

      /// @brief Binds the expression to a polytope.
      UnaryMinus& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_op->setPolytope(polytope);
        return *this;
      }

      /**
       * @brief Integrates the negated linear form.
       * @param local Local DOF index
       * @returns Negated integral value
       */
      ScalarType integrate(size_t local) override
      {
        return -m_op->integrate(local);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<OperandType> m_op;
  };

  /**
   * @brief Deduction guide for UnaryMinus of LinearFormIntegratorBase.
   */
  template <class Number>
  UnaryMinus(const LinearFormIntegratorBase<Number>&)
    -> UnaryMinus<LinearFormIntegratorBase<Number>>;

  /**
   * @brief Applies unary minus to a linear form integrator.
   * @param lfi Linear form integrator to negate
   * @returns Negated linear form integrator
   */
  template <class Number>
  UnaryMinus<LinearFormIntegratorBase<Number>> operator-(const LinearFormIntegratorBase<Number>& lfi)
  {
    return UnaryMinus(lfi);
  }

  /**
   * @ingroup UnaryMinusSpecializations
   * @brief Negation of a list of linear form integrators.
   */
  template <class Number>
  class UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>
    : public FormLanguage::List<LinearFormIntegratorBase<Number>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Number;

      /// @brief Linear form integrator base type.
      using LinearFormIntegratorBaseType = LinearFormIntegratorBase<Number>;

      /// @brief Operand type.
      using OperandType = FormLanguage::List<LinearFormIntegratorBaseType>;

      /// @brief Parent class type.
      using Parent = FormLanguage::List<LinearFormIntegratorBaseType>;

      /// @brief Constructs the expression from its operand.
      UnaryMinus(const OperandType& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<LinearFormIntegratorBaseType>(p));
      }

      /// @brief Copy constructor.
      UnaryMinus(const UnaryMinus& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };

  /**
   * @brief Deduction guide for UnaryMinus of List<LinearFormIntegratorBase>.
   */
  template <class Number>
  UnaryMinus(const FormLanguage::List<LinearFormIntegratorBase<Number>>&)
    -> UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>;

  /**
   * @brief Applies unary minus to a list of linear form integrators.
   * @param op List to negate
   * @returns Negated integrator list
   */
  template <class Number>
  UnaryMinus<FormLanguage::List<LinearFormIntegratorBase<Number>>>
  operator-(const FormLanguage::List<LinearFormIntegratorBase<Number>>& op)
  {
    return UnaryMinus(op);
  }

  /**
   * @ingroup UnaryMinusSpecializations
   * @brief Negation of a bilinear form integrator.
   */
  template <class Number>
  class UnaryMinus<LocalBilinearFormIntegratorBase<Number>>
    : public LocalBilinearFormIntegratorBase<Number>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Number;

      /// @brief Operand type.
      using OperandType = LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Parent class type.
      using Parent = LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Constructs the expression from its operand.
      UnaryMinus(const OperandType& op)
        : Parent(op),
          m_op(op.copy())
      {}

      /// @brief Copy constructor.
      UnaryMinus(const UnaryMinus& other)
        : Parent(other),
          m_op(other.m_op->copy())
      {}

      /// @brief Move constructor.
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other)),
          m_op(std::move(other.m_op))
      {}

      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        assert(m_op);
        return *m_op;
      }

      /// @brief Returns the integration region.
      Geometry::Region getRegion() const override
      {
        return getOperand().getRegion();
      }

      /// @brief Returns the polytope the expression is bound to.
      const Geometry::Polytope& getPolytope() const override
      {
        return getOperand().getPolytope();
      }

      /// @brief Binds the expression to a polytope.
      UnaryMinus& setPolytope(const Geometry::Polytope& polytope) override
      {
        m_op->setPolytope(polytope);
        return *this;
      }

      /**
       * @brief Integrates the negated bilinear form.
       * @param tr Trial DOF index
       * @param te Test DOF index
       * @returns Negated integral value
       */
      ScalarType integrate(size_t tr, size_t te) override
      {
        return -m_op->integrate(tr, te);
      }

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }

    private:
      std::unique_ptr<OperandType> m_op;
  };

  /**
   * @brief Deduction guide for UnaryMinus of LocalBilinearFormIntegratorBase.
   */
  template <class Number>
  UnaryMinus(const LocalBilinearFormIntegratorBase<Number>&)
    -> UnaryMinus<LocalBilinearFormIntegratorBase<Number>>;

  /**
   * @brief Applies unary minus to a bilinear form integrator.
   * @param op Bilinear form integrator to negate
   * @returns Negated bilinear form integrator
   */
  template <class Number>
  constexpr
  auto
  operator-(const LocalBilinearFormIntegratorBase<Number>& op)
  {
    return UnaryMinus(op);
  }

  /**
   * @ingroup UnaryMinusSpecializations
   * @brief Negation of a list of bilinear form integrators.
   */
  template <class Number>
  class UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>>
    : public FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Number;

      /// @brief Local bilinear form integrator base type.
      using LocalBilinearFormIntegratorBaseType =
        LocalBilinearFormIntegratorBase<ScalarType>;

      /// @brief Operand type.
      using OperandType =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      /// @brief Parent class type.
      using Parent =
        FormLanguage::List<LocalBilinearFormIntegratorBaseType>;

      /// @brief Constructs the expression from its operand.
      UnaryMinus(const OperandType& op)
      {
        for (const auto& p : op)
          add(UnaryMinus<LocalBilinearFormIntegratorBaseType>(p));
      }

      /// @brief Copy constructor.
      UnaryMinus(const UnaryMinus& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      UnaryMinus(UnaryMinus&& other)
        : Parent(std::move(other))
      {}

      UnaryMinus* copy() const noexcept override
      {
        return new UnaryMinus(*this);
      }
  };

  /**
   * @brief Deduction guide for UnaryMinus of List<LocalBilinearFormIntegratorBase>.
   */
  template <class Number>
  UnaryMinus(const FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>&)
    -> UnaryMinus<FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>>;

  /**
   * @brief Applies unary minus to a list of bilinear form integrators.
   * @param op List to negate
   * @returns Negated integrator list
   */
  template <class Number>
  constexpr
  auto
  operator-(const FormLanguage::List<LocalBilinearFormIntegratorBase<Number>>& op)
  {
    return UnaryMinus(op);
  }
}

#endif
