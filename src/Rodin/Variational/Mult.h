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
  template <class LHSDerived, class RHSDerived>
  class Mult<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public FunctionBase<Mult<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = FunctionBase<Mult<LHS, RHS>>;

      constexpr
      Mult(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeType() == RangeType::Scalar
            || rhs.getRangeType() == RangeType::Scalar
            || lhs.getRangeShape().width() == rhs.getRangeShape().height());
      }

      constexpr
      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
        using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          return getRHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<RHSRange, Scalar>)
        {
          return getLHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Vector> || std::is_same_v<LHSRange, Math::Matrix>)
        {
          return getLHS().getRangeShape().product(getRHS().getRangeShape());
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      Mult& traceOf(Geometry::Attribute attr)
      {
        getLHS().traceOf(attr);
        getRHS().traceOf(attr);
        return *this;
      }

      inline
      constexpr
      const LHS& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return this->object(getLHS().getValue(p)) * this->object(getRHS().getValue(p));
      }

      inline Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
  };

  template <class LHSDerived, class RHSDerived>
  Mult(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Mult<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived>
  inline
  constexpr
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class Number, class RHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator*(Number lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Mult(ScalarFunction(lhs), rhs);
  }

  template <class Number, class LHSDerived, typename = std::enable_if_t<std::is_arithmetic_v<Number>>>
  inline
  constexpr
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, Number rhs)
  {
    return Mult(lhs, ScalarFunction(rhs));
  }

  /**
   * @ingroup MultSpecializations
   * @brief Left Multiplication of a ShapeFunctionBase by a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *   f A(u)
   * @f]
   * where @f$ f @f$ is the function (scalar, vector or matrix valued), and
   * $A(u)$ is the shape function (scalar, vector, or matrix valued).
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>, FES, Space>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = ShapeFunctionBase<RHSDerived, FES, Space>;
      using Parent = ShapeFunctionBase<Mult<LHS, RHS>, FES, Space>;

      constexpr
      Mult(const LHS& lhs, const RHS& rhs)
        : Parent(rhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        // assert(lhs.getRangeType() == RangeType::Scalar
        //     || rhs.getRangeType() == RangeType::Scalar
        //     || lhs.getRangeShape().width() == rhs.getRangeShape().height());
      }

      constexpr
      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getRHS().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
        using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          return getRHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<RHSRange, Scalar>)
        {
          return getLHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Vector> || std::is_same_v<LHSRange, Math::Matrix>)
        {
          return getLHS().getRangeShape().product(getRHS().getRangeShape());
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return getRHS().getDOFs(element);
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getRHS().getFiniteElementSpace();
      }

      inline
      constexpr
      const LHS& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const auto& lhs = this->object(getLHS().getValue(p));
        const auto& rhs = this->object(getRHS().getTensorBasis(p));
        return TensorBasis(rhs.getDOFs(), [&](size_t i){ return lhs * rhs(i); });
      }

      inline Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Mult(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, Space>&)
    -> Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(Scalar lhs, const ShapeFunctionBase<RHSDerived, FES, Space>& rhs)
  {
    return Mult(ScalarFunction(lhs), rhs);
  }

  /**
   * @ingroup MultSpecializations
   * @brief Right multiplication of a ShapeFunctionBase by a FunctionBase
   *
   * Represents the following expression:
   * @f[
   *    A(u) f
   * @f]
   * where @f$ f @f$ is the function (scalar, vector or matrix valued), and
   * $A(u)$ is the shape function (scalar, vector, or matrix valued).
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Mult<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>> final
    : public ShapeFunctionBase<Mult<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>, FES, Space>
  {
    public:
      using LHS = ShapeFunctionBase<RHSDerived, FES, Space>;
      using RHS = FunctionBase<LHSDerived>;
      using Parent = ShapeFunctionBase<Mult<LHS, RHS>, FES, Space>;

      constexpr
      Mult(const LHS& lhs, const RHS& rhs)
        : Parent(lhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeType() == RangeType::Scalar
            || rhs.getRangeType() == RangeType::Scalar
            || lhs.getRangeShape().width() == rhs.getRangeShape().height());
      }

      constexpr
      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      const auto& getLeaf() const
      {
        return getRHS().getLeaf();
      }

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
        using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          return getRHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<RHSRange, Scalar>)
        {
          return getLHS().getRangeShape();
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Vector> || std::is_same_v<LHSRange, Math::Matrix>)
        {
          return getLHS().getRangeShape().product(getRHS().getRangeShape());
        }
        else
        {
          assert(false);
          return { 0, 0 };
        }
      }

      inline
      constexpr
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return getRHS().getDOFs(element);
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getRHS().getFiniteElementSpace();
      }

      inline
      constexpr
      const LHS& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        const auto& lhs = this->object(getLHS().getTensorBasis(p));
        const auto& rhs = this->object(getRHS().getValue(p));
        return TensorBasis(lhs.getDOFs(), [&](size_t i){ return lhs(i) * rhs; });
      }

      inline Mult* copy() const noexcept override
      {
        return new Mult(*this);
      }

    private:
      std::unique_ptr<LHS> m_lhs;
      std::unique_ptr<RHS> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Mult(const ShapeFunctionBase<LHSDerived, FES, Space>&, const FunctionBase<RHSDerived>&)
    -> Mult<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>;

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs, const FunctionBase<RHSDerived>& rhs)
  {
    return Mult(lhs, rhs);
  }

  template <class LHSDerived, class FES, ShapeFunctionSpaceType Space>
  inline
  constexpr
  auto
  operator*(const ShapeFunctionBase<LHSDerived, FES, Space>& lhs, Scalar rhs)
  {
    return Mult(lhs, ScalarFunction(rhs));
  }

  // /* ||--- OPTIMIZATIONS ----------------------------------------------------
  //  * Mult<FunctionBase, ShapeFunctionBase<Space>>
  //  * ---------------------------------------------------------------------->>
  //  */

  // /**
  //  * @ingroup MultSpecializations
  //  * @brief Left Multiplication of a ShapeFunction by a FunctionBase
  //  *
  //  * Represents the following expression:
  //  * @f[
  //  *   f u
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or or matrix valued).
  //  */
  // template <class FES, ShapeFunctionSpaceType Space>
  // class Mult<FunctionBase, ShapeFunction<FES, Space>>
  //   : public Mult<FunctionBase, ShapeFunctionBase<Space>>
  // {
  //   public:
  //     using Parent = Mult<FunctionBase, ShapeFunctionBase<Space>>;
  //     using LHS = FunctionBase;
  //     using RHS = ShapeFunction<FES, Space>;

  //     constexpr
  //     Mult(const FunctionBase& lhs, const ShapeFunction<FES, Space>& rhs)
  //       : Parent(lhs, rhs)
  //     {}

  //     constexpr
  //     Mult(const Mult& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Mult(Mult&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual FunctionBase& getLHS() override
  //     {
  //       return static_cast<FunctionBase&>(Parent::getLHS());
  //     }

  //     virtual const FunctionBase& getLHS() const override
  //     {
  //       return static_cast<const FunctionBase&>(Parent::getLHS());
  //     }

  //     virtual ShapeFunction<FES, Space>& getRHS() override
  //     {
  //       return static_cast<ShapeFunction<FES, Space>&>(Parent::getRHS());
  //     }

  //     virtual const ShapeFunction<FES, Space>& getRHS() const override
  //     {
  //       return static_cast<const ShapeFunction<FES, Space>&>(Parent::getRHS());
  //     }

  //     virtual Mult* copy() const noexcept override
  //     {
  //       return new Mult(*this);
  //     }
  // };
  // template <class FES, ShapeFunctionSpaceType Space>
  // Mult(const FunctionBase&, const ShapeFunction<FES, Space>& rhs)
  //   -> Mult<FunctionBase, ShapeFunction<FES, Space>>;

  // template <class FES, ShapeFunctionSpaceType Space>
  // Mult<FunctionBase, ShapeFunction<FES, Space>>
  // operator*(const FunctionBase& lhs, const ShapeFunction<FES, Space>& rhs)
  // {
  //   return Mult(lhs, rhs);
  // }

  // template <class T, class FES, ShapeFunctionSpaceType Space>
  // std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, ShapeFunction<FES, Space>>>
  // operator*(T v, const ShapeFunction<FES, Space>& rhs)
  // {
  //   return Mult(ScalarFunction(v), rhs);
  // }

  // /**
  //  * @ingroup MultSpecializations
  //  * @brief Left Multiplication of the gradient of a ShapeFunction by a FunctionBase
  //  *
  //  * Represents the following expression:
  //  * @f[
  //  *   f \nabla u
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES, ShapeFunctionSpaceType Space>
  // class Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>
  //   : public Mult<FunctionBase, ShapeFunctionBase<Space>>
  // {
  //   public:
  //     using Parent = Mult<FunctionBase, ShapeFunctionBase<Space>>;
  //     using LHS = FunctionBase;
  //     using RHS = Grad<ShapeFunction<FES, Space>>;

  //     constexpr
  //     Mult(const FunctionBase& lhs, const Grad<ShapeFunction<FES, Space>>& rhs)
  //       : Parent(lhs, rhs)
  //     {}

  //     constexpr
  //     Mult(const Mult& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Mult(Mult&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual FunctionBase& getLHS() override
  //     {
  //       return static_cast<FunctionBase&>(Parent::getLHS());
  //     }

  //     virtual const FunctionBase& getLHS() const override
  //     {
  //       return static_cast<const FunctionBase&>(Parent::getLHS());
  //     }

  //     virtual Grad<ShapeFunction<FES, Space>>& getRHS() override
  //     {
  //       return static_cast<Grad<ShapeFunction<FES, Space>>&>(Parent::getRHS());
  //     }

  //     virtual const Grad<ShapeFunction<FES, Space>>& getRHS() const override
  //     {
  //       return static_cast<const Grad<ShapeFunction<FES, Space>>&>(Parent::getRHS());
  //     }

  //     virtual Mult* copy() const noexcept override
  //     {
  //       return new Mult(*this);
  //     }
  // };
  // template <class FES, ShapeFunctionSpaceType Space>
  // Mult(const FunctionBase&, const Grad<ShapeFunction<FES, Space>>&)
  //   -> Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>;

  // template <class FES, ShapeFunctionSpaceType Space>
  // Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>
  // operator*(const FunctionBase& lhs, const Grad<ShapeFunction<FES, Space>>& rhs)
  // {
  //   return Mult(lhs, rhs);
  // }

  // template <class T, class FES, ShapeFunctionSpaceType Space>
  // std::enable_if_t<
  //   std::is_arithmetic_v<T>, Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>>
  // operator*(T v, const Grad<ShapeFunction<FES, Space>>& rhs)
  // {
  //   return Mult(ScalarFunction(v), rhs);
  // }

  // /**
  //  * @ingroup MultSpecializations
  //  * @brief Left Multiplication of the gradient of a ShapeFunction by the
  //  * Jacobian of a FunctionBase
  //  *
  //  * Represents the following expression:
  //  * @f[
  //  *   f \mathbf{J} u
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES, ShapeFunctionSpaceType Space>
  // class Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>
  //   : public Mult<FunctionBase, ShapeFunctionBase<Space>>
  // {
  //   public:
  //     using Parent = Mult<FunctionBase, ShapeFunctionBase<Space>>;
  //     using LHS = FunctionBase;
  //     using RHS = Jacobian<ShapeFunction<FES, Space>>;

  //     constexpr
  //     Mult(const FunctionBase& lhs, const Jacobian<ShapeFunction<FES, Space>>& rhs)
  //       : Parent(lhs, rhs)
  //     {}

  //     constexpr
  //     Mult(const Mult& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Mult(Mult&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual FunctionBase& getLHS() override
  //     {
  //       return static_cast<FunctionBase&>(Parent::getLHS());
  //     }

  //     virtual const FunctionBase& getLHS() const override
  //     {
  //       return static_cast<const FunctionBase&>(Parent::getLHS());
  //     }

  //     virtual Jacobian<ShapeFunction<FES, Space>>& getRHS() override
  //     {
  //       return static_cast<Jacobian<ShapeFunction<FES, Space>>&>(Parent::getRHS());
  //     }

  //     virtual const Jacobian<ShapeFunction<FES, Space>>& getRHS() const override
  //     {
  //       return static_cast<const Jacobian<ShapeFunction<FES, Space>>&>(Parent::getRHS());
  //     }

  //     virtual Mult* copy() const noexcept override
  //     {
  //       return new Mult(*this);
  //     }
  // };
  // template <class FES, ShapeFunctionSpaceType Space>
  // Mult(const FunctionBase&, const Jacobian<ShapeFunction<FES, Space>>&)
  //   -> Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>;

  // template <class FES, ShapeFunctionSpaceType Space>
  // Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>
  // operator*(const FunctionBase& lhs, const Jacobian<ShapeFunction<FES, Space>>& rhs)
  // {
  //   return Mult(lhs, rhs);
  // }

  // template <class T, class FES, ShapeFunctionSpaceType Space>
  // std::enable_if_t<
  //   std::is_arithmetic_v<T>, Mult<FunctionBase, Jacobian<ShapeFunction<FES, Space>>>>
  // operator*(T v, const Jacobian<ShapeFunction<FES, Space>>& rhs)
  // {
  //   return Mult(ScalarFunction(v), rhs);
  // }

  // /* <<-- OPTIMIZATIONS -----------------------------------------------------
  //  * Mult<FunctionBase, ShapeFunctionBase<Space>>
  //  * ----------------------------------------------------------------------||
  //  */
}

#endif
