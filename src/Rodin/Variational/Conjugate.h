/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_CONJUGATE_H
#define RODIN_VARIATIONAL_CONJUGATE_H

/**
 * @file Conjugate.h
 * @brief Complex conjugate operation for functions.
 */

#include "ForwardDecls.h"

#include "Rodin/Math/Common.h"

#include "Function.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Conjugate<Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup ConjugateSpecializations Conjugate Template Specializations
   * @brief Template specializations of the Conjugate class.
   * @see Conjugate
   */

  /**
   * @ingroup ConjugateSpecializations
   * @brief Complex conjugate of a function.
   *
   * For a complex-valued function @f$ f: \Omega \to \mathbb{C} @f$ with
   * @f$ f(x) = u(x) + iv(x) @f$, the conjugate is:
   * @f[
   *    \overline{f}(x) = u(x) - iv(x)
   * @f]
   *
   * For real-valued functions, the conjugate equals the original function.
   *
   * @tparam NestedDerived Type of the operand function
   *
   * @see Re, Im, ComplexFunction
   */
  template <class NestedDerived>
  class Conjugate<FunctionBase<NestedDerived>>
    : public FunctionBase<Conjugate<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = FunctionBase<Conjugate<OperandType>>;

      /**
       * @brief Constructs conjugate from a function.
       * @param[in] v Function to conjugate
       */
      Conjugate(const OperandType& v)
        : m_v(v.copy())
      {}

      Conjugate(const Conjugate& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      Conjugate(Conjugate&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      /**
       * @brief Evaluates the conjugate at a point.
       * @param[in] p Point at which to evaluate
       * @returns Conjugate value @f$ \overline{f(p)} @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::conj(this->object(getOperand().getValue(p)));
      }

      /**
       * @brief Gets the operand function.
       * @returns Reference to the operand
       */
      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Conjugate* copy() const noexcept override
      {
        return new Conjugate(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  /**
   * @brief CTAD for Conjugate.
   */
  template <class NestedDerived>
  Conjugate(const FunctionBase<NestedDerived>&) -> Conjugate<FunctionBase<NestedDerived>>;

  /**
   * @ingroup ConjugateSpecializations
   * @brief Specialization for shape functions.
   * @tparam NestedDerived Nested derived type
   * @tparam FES Finite element space type
   * @tparam Space Shape function space type (Trial or Test)
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType =
        ShapeFunctionBase<NestedDerived, FES, Space>;

      using RangeType =
        typename FormLanguage::Traits<OperandType>::RangeType;

      using Parent =
        ShapeFunctionBase<Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>>>;

      /**
       * @brief Constructs conjugate of a shape function.
       * @param[in] op Shape function to conjugate
       */
      constexpr
      Conjugate(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy())
      {}

      constexpr
      Conjugate(const Conjugate& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Conjugate(Conjugate&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Gets the operand shape function.
       * @returns Reference to the operand
       */
      constexpr
      const OperandType& getOperand() const
      {
        return *m_operand;
      }

      /**
       * @brief Gets the leaf shape function.
       * @returns Reference to the underlying shape function
       */
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Gets degrees of freedom count.
       * @param[in] element Element polytope
       * @returns Number of DOFs on element
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Sets the evaluation point.
       * @param[in] p Point for evaluation
       * @returns Reference to this
       */
      constexpr
      Conjugate& setPoint(const Geometry::Point& p)
      {
        m_operand->setPoint(p);
        return *this;
      }

      /**
       * @brief Gets the evaluation point.
       * @returns Current evaluation point
       */
      const Geometry::Point& getPoint() const
      {
        return m_operand->getPoint();
      }

      /**
       * @brief Evaluates conjugate of basis function.
       * @param[in] local Local DOF index
       * @returns Conjugate of basis function value
       */
      constexpr
      decltype(auto) getBasis(size_t local) const
      {
        return Math::conj(this->object(this->getOperand().getBasis(local)));
      }

      /**
       * @brief Gets the finite element space.
       * @returns Reference to the FE space
       */
      const FES& getFiniteElementSpace() const
      {
        return getOperand().getFiniteElementSpace();
      }

      Conjugate* copy() const noexcept override
      {
        return new Conjugate(*this);
      }
    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief CTAD for Conjugate on shape functions.
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Conjugate(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Conjugate<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif



