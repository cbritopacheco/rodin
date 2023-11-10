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

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Mult<
      Variational::FunctionBase<LHSDerived>,
      Variational::ShapeFunctionBase<RHSDerived, FESType, SpaceType>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;

    using LHS =
      Variational::FunctionBase<LHSDerived>;

    using RHS =
      Variational::ShapeFunctionBase<RHSDerived, FESType, SpaceType>;

    using LHSRange =
      typename FormLanguage::Traits<LHS>::RangeType;

    using RHSRange =
      typename FormLanguage::Traits<RHS>::RangeType;

    using RangeType =
      std::conditional_t<
        // If
        std::is_same_v<LHSRange, Scalar>,
        // Then
        RHSRange,
        // Else
        std::conditional_t<
          // If
          std::is_same_v<LHSRange, Math::Vector>,
          // Then
          std::conditional_t<
            // If
            std::is_same_v<RHSRange, Scalar>,
            // Then
            Math::Vector,
            // Else
            std::conditional_t<
              // If
              std::is_same_v<RHSRange, Math::Vector>,
              // Then
              Math::Matrix,
              // Else
              std::conditional_t<
                // If
                std::is_same_v<RHSRange, Math::Matrix>,
                // Then
                Math::Matrix,
                // Else
                void
              >
            >
          >,
          // Else
          std::conditional_t<
            // If
            std::is_same_v<LHSRange, Math::Matrix>,
            // Then
            Math::Matrix,
            // Else
            void
          >
        >
      >;
  };

  template <class LHSDerived, class RHSDerived, class FESType, Variational::ShapeFunctionSpaceType SpaceType>
  struct Traits<
    Variational::Mult<
      Variational::ShapeFunctionBase<LHSDerived, FESType, SpaceType>,
      Variational::FunctionBase<RHSDerived>>>
  {
    using FES = FESType;
    static constexpr Variational::ShapeFunctionSpaceType Space = SpaceType;
  };
}

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

      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      Mult(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeType() == RangeType::Scalar
            || rhs.getRangeType() == RangeType::Scalar
            || lhs.getRangeShape().width() == rhs.getRangeShape().height());
      }

      Mult(const Mult& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Mult(Mult&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
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
  template <class LHSDerived, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FESType, SpaceType>> final
    : public ShapeFunctionBase<Mult<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FESType, SpaceType>>>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      using LHS = FunctionBase<LHSDerived>;

      using RHS = ShapeFunctionBase<RHSDerived, FES, Space>;

      using Parent = ShapeFunctionBase<Mult<LHS, RHS>, FES, Space>;

      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      constexpr
      Mult(const LHS& lhs, const RHS& rhs)
        : Parent(rhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

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
      size_t getDOFs(const Geometry::Polytope& element) const
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
        const auto lhs = getLHS().getValue(p);
        const auto rhs = getRHS().getTensorBasis(p);
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
  template <class LHSDerived, class RHSDerived, class FESType, ShapeFunctionSpaceType SpaceType>
  class Mult<ShapeFunctionBase<LHSDerived, FESType, SpaceType>, FunctionBase<RHSDerived>> final
    : public ShapeFunctionBase<Mult<ShapeFunctionBase<LHSDerived, FESType, SpaceType>, FunctionBase<RHSDerived>>>
  {
    public:
      using FES = FESType;
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

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
      size_t getDOFs(const Geometry::Polytope& element) const
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
}

#endif
