/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MULT_H
#define RODIN_VARIATIONAL_MULT_H

#include <memory>
#include <type_traits>

#include "Rodin/Alert.h"
#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"
#include "ShapeFunction.h"
#include "Grad.h"

#include "Exceptions.h"

namespace Rodin::Variational
{
  /**
   * @defgroup MultSpecializations Mult Template Specializations
   * @brief Template specializations of the Mult class.
   * @see Mult
   */

  /**
   * @ingroup MultSpecializations
   * @brief Multiplication of two FunctionBase instances.
   */
  template <>
  class Mult<FunctionBase, FunctionBase> : public FunctionBase
  {
    public:
      using Parent = FunctionBase;
      using LHS = FunctionBase;
      using RHS = FunctionBase;

      Mult(const FunctionBase& lhs, const FunctionBase& rhs);

      Mult(const Mult& other);

      Mult(Mult&& other);

      RangeShape getRangeShape() const override;

      Mult& traceOf(Geometry::Attribute attr) override;

      virtual FunctionBase& getLHS()
      {
        return *m_lhs;
      }

      virtual FunctionBase& getRHS()
      {
        return *m_rhs;
      }

      virtual const FunctionBase& getLHS() const
      {
        return *m_lhs;
      }

      virtual const FunctionBase& getRHS() const
      {
        return *m_rhs;
      }

      virtual FunctionValue getValue(const Geometry::Point& p) const override
      {
        if (m_lhs->getRangeType() == RangeType::Scalar)
        {
          Scalar s = m_lhs->getValue(p);
          FunctionValue value = m_rhs->getValue(p);
          value *= s;
          return value;
        }
        else if (m_rhs->getRangeType() == RangeType::Scalar)
        {
          Scalar s = m_rhs->getValue(p);
          FunctionValue value = m_lhs->getValue(p);
          value *= s;
          return value;
        }
        else
        {
          assert(false);
          return 0.0;
        }
      }

      Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }
    private:
      std::unique_ptr<FunctionBase> m_lhs;
      std::unique_ptr<FunctionBase> m_rhs;
  };
  Mult(const FunctionBase&, const FunctionBase&) -> Mult<FunctionBase, FunctionBase>;

  Mult<FunctionBase, FunctionBase>
  operator*(const FunctionBase& lhs, const FunctionBase& rhs);

  template <class T>
  std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, FunctionBase>>
  operator*(T lhs, const FunctionBase& rhs)
  {
    return Mult(ScalarFunction(lhs), rhs);
  }

  template <class T>
  std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, FunctionBase>>
  operator*(const FunctionBase& lhs, T rhs)
  {
    return Mult(lhs, ScalarFunction(rhs));
  }

  /**
   * @ingroup MultSpecializations
   * @brief Left Multiplication by a FunctionBase of a ShapeFunctionBase
   *
   * Represents the following expression:
   * @f[
   *   f A(u)
   * @f]
   * where @f$ f @f$ is the function (scalar, vector or matrix valued), and
   * $A(u)$ is the shape operator (scalar, vector, or matrix valued).
   */
  template <ShapeFunctionSpaceType Space>
  class Mult<FunctionBase, ShapeFunctionBase<Space>>
    : public ShapeFunctionBase<Space>
  {
    public:
      using Parent = ShapeFunctionBase<Space>;
      using LHS = FunctionBase;
      using RHS = ShapeFunctionBase<Space>;

      constexpr
      Mult(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Mult(const Mult& other)
        :  ShapeFunctionBase<Space>(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Mult(Mult&& other)
        :  ShapeFunctionBase<Space>(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      const ShapeFunctionBase<Space>& getLeaf() const override
      {
        return m_rhs->getLeaf();
      }

      int getRows() const override
      {
        return m_rhs->getRangeShape().height();
      }

      int getColumns() const override
      {
        return m_rhs->getRangeShape().width();
      }

      int getDOFs(const Geometry::Simplex& element) const override
      {
        return m_rhs->getDOFs(element);
      }

      FiniteElementSpaceBase& getFiniteElementSpace() override
      {
        return m_rhs->getFiniteElementSpace();
      }

      const FiniteElementSpaceBase& getFiniteElementSpace() const override
      {
        return m_rhs->getFiniteElementSpace();
      }

      virtual FunctionBase& getLHS()
      {
        return *m_lhs;
      }

      virtual ShapeFunctionBase<Space>& getRHS()
      {
        return *m_rhs;
      }

      virtual const FunctionBase& getLHS() const
      {
        return *m_lhs;
      }

      virtual const ShapeFunctionBase<Space>& getRHS() const
      {
        return *m_rhs;
      }

      TensorBasis getOperator(
          ShapeComputator& compute, const Geometry::Point& p) const override
      {
        const auto& fe = getFiniteElementSpace().getFiniteElement(p.getSimplex());
        switch (m_lhs->getRangeType())
        {
          case RangeType::Scalar:
          {
            return m_lhs->getValue(p).scalar() * m_rhs->getOperator(compute, p);
          }
          default:
          {
            assert(false); // Unimplemented
            break;
          }
        }
      }

      virtual Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }
    private:
      std::unique_ptr<FunctionBase> m_lhs;
      std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
  };
  template <ShapeFunctionSpaceType Space>
  Mult(const FunctionBase&, const ShapeFunctionBase<Space>&)
    -> Mult<FunctionBase, ShapeFunctionBase<Space>>;

  template <ShapeFunctionSpaceType Space>
  Mult<FunctionBase, ShapeFunctionBase<Space>>
  operator*(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class T, ShapeFunctionSpaceType Space>
  std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, ShapeFunctionBase<Space>>>
  operator*(T v, const ShapeFunctionBase<Space>& rhs)
  {
    return Mult(ScalarFunction(v), rhs);
  }

  /* ||--- OPTIMIZATIONS ----------------------------------------------------
   * Mult<FunctionBase, ShapeFunctionBase<Space>>
   * ---------------------------------------------------------------------->>
   */

  /**
   * @ingroup MultSpecializations
   * @brief Left Multiplication of a ShapeFunction by a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *   f u
   * @f]
   * where @f$ f @f$ is a function (scalar or or matrix valued).
   */
  template <class FES, ShapeFunctionSpaceType Space>
  class Mult<FunctionBase, ShapeFunction<FES, Space>>
    : public Mult<FunctionBase, ShapeFunctionBase<Space>>
  {
    public:
      using Parent = Mult<FunctionBase, ShapeFunctionBase<Space>>;
      using LHS = FunctionBase;
      using RHS = ShapeFunction<FES, Space>;

      constexpr
      Mult(const FunctionBase& lhs, const ShapeFunction<FES, Space>& rhs)
        : Parent(lhs, rhs)
      {}

      constexpr
      Mult(const Mult& other)
        : Parent(other)
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other))
      {}

      virtual FunctionBase& getLHS() override
      {
        return static_cast<FunctionBase&>(Parent::getLHS());
      }

      virtual const FunctionBase& getLHS() const override
      {
        return static_cast<const FunctionBase&>(Parent::getLHS());
      }

      virtual ShapeFunction<FES, Space>& getRHS() override
      {
        return static_cast<ShapeFunction<FES, Space>&>(Parent::getRHS());
      }

      virtual const ShapeFunction<FES, Space>& getRHS() const override
      {
        return static_cast<const ShapeFunction<FES, Space>&>(Parent::getRHS());
      }

      virtual Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }
  };
  template <class FES, ShapeFunctionSpaceType Space>
  Mult(const FunctionBase&, const ShapeFunction<FES, Space>& rhs)
    -> Mult<FunctionBase, ShapeFunction<FES, Space>>;

  template <class FES, ShapeFunctionSpaceType Space>
  Mult<FunctionBase, ShapeFunction<FES, Space>>
  operator*(const FunctionBase& lhs, const ShapeFunction<FES, Space>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class T, class FES, ShapeFunctionSpaceType Space>
  std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, ShapeFunction<FES, Space>>>
  operator*(T v, const ShapeFunction<FES, Space>& rhs)
  {
    return Mult(ScalarFunction(v), rhs);
  }

  /**
   * @ingroup MultSpecializations
   * @brief Left Multiplication of the gradient of a ShapeFunction by a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *   f \nabla u
   * @f]
   * where @f$ f @f$ is a function (scalar or matrix valued).
   */
  template <class FES, ShapeFunctionSpaceType Space>
  class Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>
    : public Mult<FunctionBase, ShapeFunctionBase<Space>>
  {
    public:
      using Parent = Mult<FunctionBase, ShapeFunctionBase<Space>>;
      using LHS = FunctionBase;
      using RHS = Grad<ShapeFunction<FES, Space>>;

      constexpr
      Mult(const FunctionBase& lhs, const Grad<ShapeFunction<FES, Space>>& rhs)
        : Parent(lhs, rhs)
      {}

      constexpr
      Mult(const Mult& other)
        : Parent(other)
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other))
      {}

      virtual FunctionBase& getLHS() override
      {
        return static_cast<FunctionBase&>(Parent::getLHS());
      }

      virtual const FunctionBase& getLHS() const override
      {
        return static_cast<const FunctionBase&>(Parent::getLHS());
      }

      virtual Grad<ShapeFunction<FES, Space>>& getRHS() override
      {
        return static_cast<Grad<ShapeFunction<FES, Space>>&>(Parent::getRHS());
      }

      virtual const Grad<ShapeFunction<FES, Space>>& getRHS() const override
      {
        return static_cast<const Grad<ShapeFunction<FES, Space>>&>(Parent::getRHS());
      }

      virtual Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }
  };
  template <class FES, ShapeFunctionSpaceType Space>
  Mult(const FunctionBase&, const Grad<ShapeFunction<FES, Space>>&)
    -> Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>;

  template <class FES, ShapeFunctionSpaceType Space>
  Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>
  operator*(const FunctionBase& lhs, const Grad<ShapeFunction<FES, Space>>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class T, class FES, ShapeFunctionSpaceType Space>
  std::enable_if_t<
    std::is_arithmetic_v<T>, Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>>
  operator*(T v, const Grad<ShapeFunction<FES, Space>>& rhs)
  {
    return Mult(ScalarFunction(v), rhs);
  }

  /**
   * @ingroup MultSpecializations
   * @brief Left Multiplication of the gradient of a ShapeFunction by the
   * Jacobian of a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *   f \mathbf{J} u
   * @f]
   * where @f$ f @f$ is a function (scalar or matrix valued).
   */
  template <class FES, ShapeFunctionSpaceType Space>
  class Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>
    : public Mult<FunctionBase, ShapeFunctionBase<Space>>
  {
    public:
      using Parent = Mult<FunctionBase, ShapeFunctionBase<Space>>;
      using LHS = FunctionBase;
      using RHS = Jacobian<ShapeFunction<FES, Space>>;

      constexpr
      Mult(const FunctionBase& lhs, const Jacobian<ShapeFunction<FES, Space>>& rhs)
        : Parent(lhs, rhs)
      {}

      constexpr
      Mult(const Mult& other)
        : Parent(other)
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other))
      {}

      virtual FunctionBase& getLHS() override
      {
        return static_cast<FunctionBase&>(Parent::getLHS());
      }

      virtual const FunctionBase& getLHS() const override
      {
        return static_cast<const FunctionBase&>(Parent::getLHS());
      }

      virtual Jacobian<ShapeFunction<FES, Space>>& getRHS() override
      {
        return static_cast<Jacobian<ShapeFunction<FES, Space>>&>(Parent::getRHS());
      }

      virtual const Jacobian<ShapeFunction<FES, Space>>& getRHS() const override
      {
        return static_cast<const Jacobian<ShapeFunction<FES, Space>>&>(Parent::getRHS());
      }

      virtual Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }
  };
  template <class FES, ShapeFunctionSpaceType Space>
  Mult(const FunctionBase&, const Jacobian<ShapeFunction<FES, Space>>&)
    -> Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>;

  template <class FES, ShapeFunctionSpaceType Space>
  Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>
  operator*(const FunctionBase& lhs, const Jacobian<ShapeFunction<FES, Space>>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class T, class FES, ShapeFunctionSpaceType Space>
  std::enable_if_t<
    std::is_arithmetic_v<T>, Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>>
  operator*(T v, const Jacobian<ShapeFunction<FES, Space>>& rhs)
  {
    return Mult(ScalarFunction(v), rhs);
  }

  /* <<-- OPTIMIZATIONS -----------------------------------------------------
   * Mult<FunctionBase, ShapeFunctionBase<Space>>
   * ----------------------------------------------------------------------||
   */
}

#endif
