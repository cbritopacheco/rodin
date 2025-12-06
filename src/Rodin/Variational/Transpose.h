/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Matrix and ShapeFunction transpose operations.
 */

#ifndef RODIN_VARIATIONAL_TRANSPOSE_H
#define RODIN_VARIATIONAL_TRANSPOSE_H

#include "ShapeFunction.h"
#include "MatrixFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup TransposeSpecializations Transpose Template Specializations
   * @brief Template specializations of the Transpose class.
   * @see Transpose
   */

  /**
   * @brief Transpose of a matrix-valued function.
   *
   * Given a matrix function @f$ A: \Omega \to \mathbb{R}^{m \times n} @f$,
   * computes its transpose:
   * @f[
   *    A^T(x) = (A(x))^T \in \mathbb{R}^{n \times m}
   * @f]
   *
   * @tparam NestedDerived Type of the matrix function
   * @ingroup TransposeSpecializations
   */
  template <class NestedDerived>
  class Transpose<FunctionBase<NestedDerived>> final
    : public FunctionBase<Transpose<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using Parent = FunctionBase<Transpose<OperandType>>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      /**
       * @brief Constructs the transpose of a matrix function.
       * @param m Matrix function to transpose
       */
      constexpr
      Transpose(const OperandType& m)
        : m_operand(m.copy())
      {}

      constexpr
      Transpose(const Transpose& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Transpose(Transpose&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Gets the underlying matrix function.
       * @returns Reference to the operand
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Evaluates the transpose at a point.
       * @param p Point at which to evaluate
       * @returns Transposed matrix @f$ A^T(p) @f$
       */
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getOperand().getValue(p)).transpose();
      }

      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for function transpose.
   */
  template <class NestedDerived>
  Transpose(const FunctionBase<NestedDerived>&) -> Transpose<FunctionBase<NestedDerived>>;

  /**
   * @brief Transpose of a matrix-valued ShapeFunction.
   *
   * Computes the transpose of matrix trial or test functions in
   * finite element formulations.
   *
   * @tparam NestedDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Trial or test function space
   * @ingroup TransposeSpecializations
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  class Transpose<ShapeFunctionBase<NestedDerived, FES, Space>> final
    : public ShapeFunctionBase<Transpose<ShapeFunctionBase<NestedDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;

      using OperandType = ShapeFunctionBase<NestedDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Transpose<OperandType>, FES, Space>;

      /**
       * @brief Constructs the transpose of a ShapeFunction.
       * @param op Shape function to transpose
       */
      constexpr
      Transpose(const OperandType& op)
        : Parent(op.getFiniteElementSpace()),
          m_operand(op.copy())
      {}

      constexpr
      Transpose(const Transpose& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}

      constexpr
      Transpose(Transpose&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}

      /**
       * @brief Gets the leaf (underlying trial/test function).
       * @returns Reference to the leaf function
       */
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      /**
       * @brief Gets number of degrees of freedom on a polytope.
       * @param simplex Mesh polytope
       * @returns Number of DOFs
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& simplex) const
      {
        return getOperand().getDOFs(simplex);
      }

      /**
       * @brief Gets the underlying shape function.
       * @returns Reference to the operand
       */
      constexpr
      const OperandType& getOperand() const
      {
        assert(m_operand);
        return *m_operand;
      }

      /**
       * @brief Gets the current evaluation point.
       * @returns Reference to the point
       */
      const Geometry::Point& getPoint() const
      {
        return m_operand->getPoint();
      }

      /**
       * @brief Sets the evaluation point.
       * @param p Point to evaluate at
       * @returns Reference to this object
       */
      Transpose& setPoint(const Geometry::Point& p)
      {
        m_operand->setPoint(p);
        return *this;
      }

      /**
       * @brief Gets the transposed basis function for local DOF.
       * @param local Local DOF index
       * @returns Transposed basis function matrix
       */
      constexpr
      auto getBasis(size_t local) const
      {
        return this->object(getOperand().getBasis(local)).transpose();
      }

      /**
       * @brief Gets the finite element space.
       * @returns Reference to the FE space
       */
      constexpr
      const FES& getFiniteElementSpace() const
      {
        return m_operand.getFiniteElementSpace();
      }

      Transpose* copy() const noexcept override
      {
        return new Transpose(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
  };

  /**
   * @brief Deduction guide for ShapeFunction transpose.
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Transpose(const ShapeFunctionBase<NestedDerived, FES, Space>&)
    -> Transpose<ShapeFunctionBase<NestedDerived, FES, Space>>;
}

#endif
