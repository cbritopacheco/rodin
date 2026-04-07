/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FROBENIUS_H
#define RODIN_VARIATIONAL_FROBENIUS_H

/**
 * @file
 * @brief Frobenius norm operation for functions.
 *
 * Computes the Frobenius norm (Euclidean norm) of vectors and matrices:
 * - For scalars: @f$ \|x\| = |x| @f$
 * - For vectors: @f$ \|\mathbf{v}\| = \sqrt{\sum_i v_i^2} @f$
 * - For matrices: @f$ \|A\|_F = \sqrt{\sum_{ij} A_{ij}^2} @f$
 */

#include "Rodin/Math/Common.h"

#include "ForwardDecls.h"
#include "RealFunction.h"
#include "Function.h"

namespace Rodin::Variational
{
  /**
   * @defgroup FrobeniusSpecializations Frobenius Template Specializations
   * @brief Template specializations of the Frobenius class.
   * @see Frobenius
   */

  /**
   * @ingroup FrobeniusSpecializations
   * @brief Frobenius norm of a function.
   *
   * Computes the Frobenius (Euclidean) norm:
   * - For scalars: @f$ \|\cdot\| : \mathbb{R} \to \mathbb{R}, \quad \|x\| = |x| @f$
   * - For vectors @f$ \mathbf{v} \in \mathbb{R}^n @f$:
   *   @f[
   *      \|\mathbf{v}\| = \sqrt{\sum_{i=1}^{n} v_i^2}
   *   @f]
   * - For matrices @f$ A \in \mathbb{R}^{m \times n} @f$:
   *   @f[
   *      \|A\|_F = \sqrt{\sum_{i=1}^{m}\sum_{j=1}^{n} A_{ij}^2}
   *   @f]
   *
   * @tparam NestedDerived Type of the function to take the norm of
   *
   * @note The result is always a real-valued scalar function.
   */
  template <class NestedDerived>
  class Frobenius<FunctionBase<NestedDerived>>
    : public RealFunctionBase<Frobenius<FunctionBase<NestedDerived>>>
  {
    public:
      using OperandType = FunctionBase<NestedDerived>;

      using OperandRangeType = typename FormLanguage::Traits<OperandType>::RangeType;

      using ScalarType = typename FormLanguage::Traits<OperandType>::ScalarType;

      using Parent = RealFunctionBase<Frobenius<OperandType>>;

      /**
       * @brief Constructs a Frobenius norm from a function.
       * @param v Function to compute the norm of
       */
      Frobenius(const OperandType& v)
        : m_v(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param other Frobenius norm object to copy
       */
      Frobenius(const Frobenius& other)
        : Parent(other),
          m_v(other.m_v->copy())
      {}

      /**
       * @brief Move constructor.
       * @param other Frobenius norm object to move from
       */
      Frobenius(Frobenius&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v))
      {}

      /**
       * @brief Evaluates the Frobenius norm at a point.
       * @param p Point at which to evaluate
       * @returns Norm value @f$ \|f(x)\| @f$ at point @f$ x @f$
       */
      auto getValue(const Geometry::Point& p) const
      {
        if constexpr (std::is_same_v<OperandRangeType, Real>)
        {
          return Math::abs(this->getOperand().getValue(p));
        }
        else
        {
          static thread_local OperandRangeType s_v;
          s_v = this->getOperand().getValue(p);
          return s_v.norm();
        }
      }

      /**
       * @brief Gets the operand function.
       * @returns Reference to the function whose norm is computed
       */
      const OperandType& getOperand() const
      {
        assert(m_v);
        return *m_v;
      }

      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return GetOrderIfConstant(getOperand(), polytope);
      }

      /**
       * @brief Creates a copy of the Frobenius norm operation.
       * @returns Pointer to copied object
       */
      Frobenius* copy() const noexcept override
      {
        return new Frobenius(*this);
      }

    private:
      std::unique_ptr<OperandType> m_v;
  };

  /**
   * @brief Deduction guide for Frobenius norm.
   */
  template <class NestedDerived>
  Frobenius(const FunctionBase<NestedDerived>&) -> Frobenius<FunctionBase<NestedDerived>>;
}

#endif
