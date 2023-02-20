/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_H
#define RODIN_VARIATIONAL_DOT_H

#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup DotSpecializations Dot Template Specializations
   * @brief Template specializations of the Dot class.
   * @see Dot
   */

  /**
   * @ingroup DotSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class Dot<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public ScalarFunctionBase<Dot<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = FunctionBase<RHSDerived>;
      using Parent = ScalarFunctionBase<Dot<LHS, RHS>>;

      constexpr
      Dot(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

      constexpr
      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      Dot& traceOf(Geometry::Attribute attrs)
      {
        m_lhs.traceOf(attrs);
        m_rhs.traceOf(attrs);
        return *this;
      }

      inline
      constexpr
      LHS& getLHS()
      {
        return m_lhs;
      }

      inline
      constexpr
      RHS& getRHS()
      {
        return m_rhs;
      }

      inline
      constexpr
      const LHS& getLHS() const
      {
        return m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        return m_rhs;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        assert(m_lhs.getRangeShape() == m_rhs.getRangeShape());
        using LHSRange = FormLanguage::RangeOf<typename FormLanguage::Traits<LHS>::ResultType>;
        using RHSRange = FormLanguage::RangeOf<typename FormLanguage::Traits<RHS>::ResultType>;
        static_assert(std::is_same_v<LHSRange, RHSRange>);
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          return m_lhs.getValue(p) * m_rhs.getValue(p);
        } else if constexpr (std::is_same_v<LHSRange, Math::Vector>)
        {
          return m_lhs.getValue(p).dot(m_rhs.getValue(p));
        }
        else if constexpr (std::is_same_v<RHSRange, Math::Matrix>)
        {
          return (m_lhs.getValue(p) * m_rhs.getValue(p).transpose()).trace();
        }
        else
        {
          assert(false);
          return void();
        }
      }

    private:
      FunctionBase<LHSDerived> m_lhs;
      FunctionBase<RHSDerived> m_rhs;
  };
  template <class LHSDerived, class RHSDerived>
  Dot(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Dot<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @ingroup DotSpecializations
   * @brief Dot product between a FunctionBase and a ShapeFunctionBase.
   *
   * Represents the mathematical expression:
   * @f[
   *   \Lambda : A(u)
   * @f]
   * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$, @f$ \Lambda \in
   * \mathbb{R}^{p \times q} @f$.
   */
  template <class LHSDerived, class RHSDerived, ShapeFunctionSpaceType Space>
  class Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, Space>> final
    : public ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, Space>>, Space>
  {
    public:
      using LHS = FunctionBase<LHSDerived>;
      using RHS = ShapeFunctionBase<RHSDerived, Space>;
      using Parent =
        ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, Space>>, Space>;

      constexpr
      Dot(const LHS& lhs, const RHS& rhs)
        : m_lhs(lhs), m_rhs(rhs)
      {
        assert(m_lhs.getRangeShape() == m_rhs.getRangeShape());
      }

      constexpr
      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs), m_rhs(other.m_rhs)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      inline
      constexpr
      const LHS& getLHS() const
      {
        return m_lhs;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        return m_rhs;
      }

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
        return { 1, 1 };
      }

      inline
      size_t getDOFs(const Geometry::Simplex& element) const
      {
        return getRHS().getDOFs(element);
      }

      inline
      constexpr
      auto& getFiniteElementSpace()
      {
        return m_rhs.getFiniteElementSpace();
      }

      inline
      constexpr
      const auto& getFiniteElementSpace() const
      {
        return getRHS().getFiniteElementSpace();
      }

      inline
      constexpr
      auto getOperator(ShapeComputator& compute, const Geometry::Point& p) const
      {
        assert(m_lhs.getRangeShape() == m_rhs.getRangeShape());
        using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
        using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
        static_assert(std::is_same_v<LHSRange, RHSRange>);
        const auto lhs = getLHS().getValue(p);
        const auto rhs = getRHS().getOperator(compute, p);
        if constexpr (std::is_same_v<LHSRange, Scalar>)
        {
          return lhs * rhs;
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Vector>)
        {
          return TensorBasis([&](size_t i) { return lhs.dot(rhs(i)); });
        }
        else if constexpr (std::is_same_v<LHSRange, Math::Matrix>)
        {
          return TensorBasis([&](size_t i) { return (lhs * rhs(i).transpose()).trace(); });
        }
        else
        {
          assert(false);
          return void();
        }
      }

    private:
      FunctionBase<LHSDerived> m_lhs;
      ShapeFunctionBase<RHSDerived, Space> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, ShapeFunctionSpaceType Space>
  Dot(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, Space>&)
    -> Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, Space>>;

  /**
   * @ingroup DotSpecializations
   */
  template <class LHSDerived, class RHSDerived>
  class Dot<ShapeFunctionBase<LHSDerived, TrialSpace>, ShapeFunctionBase<RHSDerived, TestSpace>> final
    : public FormLanguage::Base
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialSpace>;
      using RHS = ShapeFunctionBase<RHSDerived, TestSpace>;
      using Parent = FormLanguage::Base;

      constexpr
      Dot(const LHS& lhs, const RHS& rhs)
        : m_trial(lhs), m_test(rhs)
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

      constexpr
      Dot(const Dot& other)
        : Base(other),
          m_trial(other.m_trial), m_test(other.m_test)
      {}

      constexpr
      Dot(Dot&& other)
        : Base(std::move(other)),
          m_trial(std::move(other.m_trial)), m_test(std::move(other.m_test))
      {}

      inline
      constexpr
      const LHS& getLHS() const
      {
        return m_trial;
      }

      inline
      constexpr
      const RHS& getRHS() const
      {
        return m_test;
      }

      Math::Matrix getMatrix(ShapeComputator& compute, const Geometry::Point& p) const
      {
        assert(m_trial.getRangeShape() == m_test.getRangeShape());
        const auto trial = m_trial.getOperator(compute, p);
        const auto test = m_test.getOperator(compute, p);
        Math::Matrix res(test.getDOFs(), trial.getDOFs());
        for (Math::Matrix::Index i = 0; i < res.rows(); i++)
          for (Math::Matrix::Index j = 0; j < res.cols(); j++)
            res(i, j) = (trial(i) * test(j).transpose()).trace();
        return res;
      }

      inline Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }

    private:
      LHS m_trial;
      RHS m_test;
  };

  template <class LHSDerived, class RHSDerived>
  Dot(const ShapeFunctionBase<LHSDerived, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestSpace>&)
    -> Dot<ShapeFunctionBase<LHSDerived, TrialSpace>, ShapeFunctionBase<RHSDerived, TestSpace>>;

  // /* ||-- OPTIMIZATIONS -----------------------------------------------------
  //  * Dot<FunctionBase, ShapeFunctionBase<Space>>
  //  * ---------------------------------------------------------------------->>
  //  */

  // /**
  //  * @ingroup DotSpecializations
  //  * @f[
  //  *   f \cdot u
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or vector valued).
  //  */
  // template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  // class Dot<FunctionBase<LHSDerived>, ShapeFunction<RHSDerived, FES, Space>> final
  //   : public Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<ShapeFunction<RHSDerived, FES, Space>, Space>>
  // {
  //   public:
  //     using LHS = FunctionBase<LHSDerived>;
  //     using RHS = ShapeFunction<RHSDerived, FES, Space>;
  //     using Parent = Dot<FunctionBase<LHS>, ShapeFunctionBase<RHS, Space>>;

  //     constexpr
  //     Dot(const LHS& f, const RHS& u)
  //       : Parent(f, u)
  //     {}

  //     constexpr
  //     Dot(const Dot& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Dot(Dot&& other)
  //       : Parent(other)
  //     {}

  //     inline
  //     constexpr
  //     RHS& getRHS()
  //     {
  //       return static_cast<RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     constexpr
  //     const LHS& getLHS() const
  //     {
  //       return static_cast<const LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     const RHS& getRHS() const
  //     {
  //       return static_cast<const RHS&>(Parent::getRHS());
  //     }
  // };

  // template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  // Dot(const FunctionBase<LHSDerived>&, const ShapeFunction<RHSDerived, FES, Space>&)
  //   -> Dot<FunctionBase<LHSDerived>, ShapeFunction<RHSDerived, FES, Space>>;

  /* <<-- OPTIMIZATIONS -----------------------------------------------------
   * Dot<FunctionBase, ShapeFunctionBase<TestSpace>>
   * ----------------------------------------------------------------------||
   */

  // /* ||-- OPTIMIZATIONS -----------------------------------------------------
  //  * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  //  * ---------------------------------------------------------------------->>
  //  */

  // /**
  //  * @ingroup DotSpecializations
  //  *
  //  * @f[
  //  *   (f u) \cdot v
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or vector valued).
  //  */
  // template <class LHSDerived1, class LHSDerived2, class RHSDerived, class FES>
  // class Dot<
  //         Mult<FunctionBase<LHSDerived1>, ShapeFunction<LHSDerived2, FES, TrialSpace>>,
  //         ShapeFunction<RHSDerived, FES, TestSpace>> final
  //   : public Dot<
  //         ShapeFunctionBase<Mult<FunctionBase<LHSDerived1>, ShapeFunction<LHSDerived2, FES, TrialSpace>>, TrialSpace>,
  //         ShapeFunctionBase<ShapeFunction<RHSDerived, FES, TestSpace>, TestSpace>>
  // {
  //   public:
  //     using LHS = Mult<FunctionBase<LHSDerived1>, ShapeFunction<LHSDerived2, FES, TrialSpace>>;
  //     using RHS = ShapeFunction<RHSDerived, FES, TestSpace>;
  //     using Parent = Dot<ShapeFunctionBase<LHS, TrialSpace>, ShapeFunctionBase<RHS, TestSpace>>;

  //     constexpr
  //     Dot(const LHS& fu, const RHS& v)
  //       : Parent(fu, v)
  //     {}

  //     constexpr
  //     Dot(const Dot& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Dot(Dot&& other)
  //       : Parent(other)
  //     {}

  //     inline
  //     constexpr
  //     LHS& getLHS()
  //     {
  //       return static_cast<LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     RHS& getRHS()
  //     {
  //       return static_cast<RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     constexpr
  //     const LHS& getLHS() const
  //     {
  //       return static_cast<const LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     const RHS& getRHS() const
  //     {
  //       return static_cast<const RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     Dot* copy() const noexcept
  //     override
  //     {
  //       return new Dot(*this);
  //     }
  // };
  // template <class LHSDerived1, class LHSDerived2, class RHSDerived, class FES>
  // Dot(const Mult<FunctionBase<LHSDerived1>, ShapeFunction<LHSDerived2, FES, TrialSpace>>&,
  //     const ShapeFunctionBase<RHSDerived, TestSpace>&)
  //   -> Dot<
  //       Mult<FunctionBase<LHSDerived1>, ShapeFunction<LHSDerived2, FES, TrialSpace>>,
  //       ShapeFunction<RHSDerived, FES, TestSpace>>;

  // /**
  //  * @ingroup DotSpecializations
  //  *
  //  * @f[
  //  *   (f \nabla u) \cdot \nabla v
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class LHSDerived1, class LHSDerived2, class RHSDerived, class FES>
  // class Dot<
  //       Mult<FunctionBase<LHSDerived1>, Grad<ShapeFunction<LHSDerived2, FES, TrialSpace>>>,
  //       Grad<ShapeFunction<RHSDerived, FES, TestSpace>>> final
  //   : public Dot<
  //       ShapeFunctionBase<Mult<FunctionBase<LHSDerived1>, Grad<ShapeFunction<LHSDerived2, FES, TrialSpace>>>, TrialSpace>,
  //       ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, FES, TestSpace>>, TestSpace>>
  // {
  //   public:
  //     using LHS = Mult<FunctionBase<LHSDerived1>, Grad<ShapeFunction<LHSDerived2, FES, TrialSpace>>>;
  //     using RHS = Grad<ShapeFunction<RHSDerived, FES, TestSpace>>;
  //     using Parent = Dot<ShapeFunctionBase<LHS, TrialSpace>, ShapeFunctionBase<RHS, TestSpace>>;

  //     constexpr
  //     Dot(const LHS& fgu, const RHS& gv)
  //       : Parent(fgu, gv)
  //     {}

  //     constexpr
  //     Dot(const Dot& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Dot(Dot&& other)
  //       : Parent(other)
  //     {}

  //     inline
  //     constexpr
  //     LHS& getLHS()
  //     {
  //       return static_cast<LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     RHS& getRHS()
  //     {
  //       return static_cast<RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     constexpr
  //     const LHS& getLHS() const
  //     {
  //       return static_cast<const LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     const RHS& getRHS() const
  //     {
  //       return static_cast<const RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     Dot* copy() const noexcept
  //     override
  //     {
  //       return new Dot(*this);
  //     }
  // };
  // template <class LHSDerived1, class LHSDerived2, class RHSDerived, class FES>
  // Dot(const Mult<FunctionBase<LHSDerived1>, Grad<ShapeFunction<LHSDerived2, FES, TrialSpace>>>&,
  //     const Grad<ShapeFunction<RHSDerived, FES, TestSpace>>&)
  //   -> Dot<
  //       Mult<FunctionBase<LHSDerived1>, Grad<ShapeFunction<LHSDerived2, FES, TrialSpace>>>,
  //       Grad<ShapeFunction<RHSDerived, FES, TestSpace>>>;

  // /**
  //  * @ingroup DotSpecializations
  //  *
  //  * Represents the following expression:
  //  * @f[
  //  *   \nabla u \cdot \nabla v
  //  * @f]
  //  */
  // template <class LHSDerived, class RHSDerived, class FES>
  // class Dot<
  //     Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>,
  //     Grad<ShapeFunction<RHSDerived, FES, TestSpace>>> final
  //   : public Dot<
  //       ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>, TrialSpace>,
  //       ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, FES, TestSpace>>, TestSpace>>
  // {
  //   public:
  //     using LHS = Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>;
  //     using RHS = Grad<ShapeFunction<RHSDerived, FES, TestSpace>>;
  //     using Parent = Dot<ShapeFunctionBase<LHS, TrialSpace>, ShapeFunctionBase<RHS, TestSpace>>;

  //     constexpr
  //     Dot(const LHS& nu, const RHS& nv)
  //       : Parent(nu, nv)
  //     {}

  //     constexpr
  //     Dot(const Dot& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Dot(Dot&& other)
  //       : Parent(std::move(other))
  //     {}

  //     inline
  //     constexpr
  //     LHS& getLHS()
  //     {
  //       return static_cast<LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     RHS& getRHS()
  //     {
  //       return static_cast<RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     constexpr
  //     const LHS& getLHS() const
  //     {
  //       return static_cast<const LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     const RHS& getRHS() const
  //     {
  //       return static_cast<const RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     Dot* copy() const noexcept
  //     override
  //     {
  //       return new Dot(*this);
  //     }
  // };
  // template <class LHSDerived, class RHSDerived, class FES>
  // Dot(const Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>&, const Grad<ShapeFunction<RHSDerived, FES, TestSpace>>&)
  //   -> Dot<Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>, Grad<ShapeFunction<RHSDerived, FES, TestSpace>>>;

  // /**
  //  * @ingroup DotSpecializations
  //  *
  //  * Represents the following expression:
  //  * @f[
  //  *   \mathbf{J} u \cdot \mathbf{J} v
  //  * @f]
  //  */
  // template <class LHSDerived, class RHSDerived, class FES>
  // class Dot<
  //     Jacobian<ShapeFunction<LHSDerived, FES, TrialSpace>>,
  //     Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>> final
  //   : public Dot<
  //       ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, FES, TrialSpace>>, TrialSpace>,
  //       ShapeFunctionBase<Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>, TestSpace>>
  // {
  //   public:
  //     using LHS = Jacobian<ShapeFunction<LHSDerived, FES, TrialSpace>>;
  //     using RHS = Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>;
  //     using Parent = Dot<ShapeFunctionBase<LHS, TrialSpace>, ShapeFunctionBase<RHS, TestSpace>>;

  //     constexpr
  //     Dot(const LHS& nu, const RHS& nv)
  //       : Parent(nu, nv)
  //     {}

  //     constexpr
  //     Dot(const Dot& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Dot(Dot&& other)
  //       : Parent(std::move(other))
  //     {}

  //     inline
  //     constexpr
  //     LHS& getLHS()
  //     {
  //       return static_cast<LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     RHS& getRHS()
  //     {
  //       return static_cast<RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     constexpr
  //     const LHS& getLHS() const
  //     {
  //       return static_cast<const LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     const RHS& getRHS() const
  //     {
  //       return static_cast<const RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     Dot* copy() const noexcept
  //     override
  //     {
  //       return new Dot(*this);
  //     }
  // };
  // template <class LHSDerived, class RHSDerived, class FES>
  // Dot(const Jacobian<ShapeFunction<LHSDerived, FES, TrialSpace>>&,
  //     const Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>&)
  //   -> Dot<
  //       Jacobian<ShapeFunction<LHSDerived, FES, TrialSpace>>,
  //       Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>>;

  // /**
  //  * @ingroup DotSpecializations
  //  *
  //  * @f[
  //  *   (f \mathbf{J} u) \cdot \mathbf{J} v
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class LHSDerived1, class LHSDerived2, class RHSDerived, class FES>
  // class Dot<
  //     Mult<FunctionBase<LHSDerived1>, Jacobian<ShapeFunction<LHSDerived2, FES, TrialSpace>>>,
  //     Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>> final
  //   : public Dot<
  //       ShapeFunctionBase<Mult<FunctionBase<LHSDerived1>, Jacobian<ShapeFunction<LHSDerived2, FES, TrialSpace>>>, TrialSpace>,
  //       ShapeFunctionBase<Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>, TestSpace>>
  // {
  //   public:
  //     using LHS = Mult<FunctionBase<LHSDerived1>, Jacobian<ShapeFunction<LHSDerived2, FES, TrialSpace>>>;
  //     using RHS = Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>;
  //     using Parent = Dot<ShapeFunctionBase<LHS, TrialSpace>, ShapeFunctionBase<RHS, TestSpace>>;

  //     constexpr
  //     Dot(const LHS& fgu, const RHS& gv)
  //       : Parent(fgu, gv)
  //     {}

  //     constexpr
  //     Dot(const Dot& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Dot(Dot&& other)
  //       : Parent(other)
  //     {}

  //     inline
  //     constexpr
  //     LHS& getLHS()
  //     {
  //       return static_cast<LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     RHS& getRHS()
  //     {
  //       return static_cast<RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     constexpr
  //     const LHS& getLHS() const
  //     {
  //       return static_cast<const LHS&>(Parent::getLHS());
  //     }

  //     inline
  //     constexpr
  //     const RHS& getRHS() const
  //     {
  //       return static_cast<const RHS&>(Parent::getRHS());
  //     }

  //     inline
  //     Dot* copy() const noexcept
  //     override
  //     {
  //       return new Dot(*this);
  //     }
  // };
  // template <class LHSDerived1, class LHSDerived2, class RHSDerived, class FES>
  // Dot(const Mult<FunctionBase<LHSDerived1>, Jacobian<ShapeFunction<LHSDerived2, FES, TrialSpace>>>&,
  //     const Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>&)
  //   -> Dot<
  //       Mult<FunctionBase<LHSDerived1>, Jacobian<ShapeFunction<LHSDerived2, FES, TrialSpace>>>,
  //       Jacobian<ShapeFunction<RHSDerived, FES, TestSpace>>>;

  // /* <<-- OPTIMIZATIONS -----------------------------------------------------
  //  * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  //  * ----------------------------------------------------------------------||
  //  */
}

#endif
