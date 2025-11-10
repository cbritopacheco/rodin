/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Trace.h
 * @brief Matrix trace operator for matrix-valued functions.
 *
 * This file defines the Trace class, which computes the trace (sum of diagonal
 * elements) of matrix-valued functions in variational formulations.
 *
 * ## Mathematical Foundation
 * For a square matrix @f$ A \in \mathbb{R}^{n \times n} @f$, the trace is:
 * @f[
 *   \text{tr}(A) = \sum_{i=1}^n A_{ii}
 * @f]
 *
 * ## Properties
 * The trace operator satisfies:
 * - **Linearity**: @f$ \text{tr}(A + B) = \text{tr}(A) + \text{tr}(B) @f$
 * - **Scalar multiplication**: @f$ \text{tr}(\alpha A) = \alpha \text{tr}(A) @f$
 * - **Cyclic property**: @f$ \text{tr}(ABC) = \text{tr}(CAB) = \text{tr}(BCA) @f$
 * - **Transpose invariance**: @f$ \text{tr}(A^T) = \text{tr}(A) @f$
 *
 * ## Applications
 * - Linear elasticity: @f$ \text{tr}(\boldsymbol{\varepsilon}) @f$ (volumetric strain)
 * - Continuum mechanics: invariants of stress/strain tensors
 * - Fluid dynamics: divergence from velocity gradient
 * - General tensor operations
 *
 * ## Usage Example
 * ```cpp
 * // Volumetric strain in linear elasticity
 * auto strain = 0.5 * (Jacobian(u) + Transpose(Jacobian(u)));
 * auto volumetric_strain = Trace(strain);  // tr(ε)
 * ```
 */
#ifndef RODIN_VARIATIONAL_TRACE_H
#define RODIN_VARIATIONAL_TRACE_H

#include "ForwardDecls.h"

#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Trace<Variational::ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup TraceSpecializations Trace Template Specializations
   * @brief Template specializations of the Trace class.
   * @see Trace
   */

  /**
   * @ingroup TraceSpecializations
   * @brief Trace of a FunctionBase instance.
   */
  template <class NestedDerived>
  class Trace<FunctionBase<NestedDerived>> final
    : public FunctionBase<Trace<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = FunctionBase<Trace<OperandType>>;

      /**
       * @brief Constructs the trace of a matrix function.
       * @param[in] m Square matrix function
       *
       * Creates the trace operator @f$ \text{tr}(A) = \sum_{i=1}^n A_{ii} @f$
       * for a square matrix @f$ A @f$.
       */
      constexpr
      Trace(const OperandType& m)
        : m_operand(m.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Trace to copy
       */
      constexpr
      Trace(const Trace& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Trace to move from
       */
      constexpr
      Trace(Trace&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Evaluates the trace at a point.
       * @param[in] p Point at which to evaluate
       * @return Trace value @f$ \text{tr}(A(p)) = \sum_i A_{ii}(p) @f$
       *
       * Computes the sum of diagonal entries of the matrix at the given point.
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->getOperand().getValue(p).trace();
      }

      /**
       * @brief Gets the operand matrix function.
       * @return Reference to the matrix being traced
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Sets the trace domain (implementation-specific).
       * @param[in] args Arguments forwarded to operand's traceOf method
       * @return Reference to this trace operator
       */
      template <class ... Args>
      constexpr
      Trace& traceOf(const Args& ... args)
      {
        m_operand->traceOf(args...);
        return *this;
      }

      /**
       * @brief Creates a polymorphic copy of this trace operator.
       * @return Pointer to a new copy
       */
      Trace* copy() const noexcept override
      {
        return new Trace(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived>
  Trace(const FunctionBase<NestedDerived>&) -> Trace<FunctionBase<NestedDerived>>;

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Trace<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Trace<ShapeFunctionBase<NestedDerived, FES, Space>>>
  {
    public:
      using FESType = FES;
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Trace<OperandType>>;

      constexpr
      Trace(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_operand(operand.copy())
      {}

      constexpr
      Trace(const Trace& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Trace(Trace&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      constexpr
      const OperandType& getOperand() const
      {
        return *m_operand;
      }

      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      const Geometry::Point& getPoint() const
      {
        return m_operand->getPoint();
      }

      Trace& setPoint(const Geometry::Point& p)
      {
        m_operand->setPoint(p);
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        return this->getOperand().getBasis(local).transpose();
      }

      const FES& getFiniteElementSpace() const
      {
        return this->getOperand().getFiniteElementSpace();
      }

      Trace* copy() const noexcept override
      {
        return new Trace(*this);
      }
    private:
      std::unique_ptr<OperandType> m_operand;
  };

  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Trace(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Trace<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
