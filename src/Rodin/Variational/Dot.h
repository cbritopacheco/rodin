/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_H
#define RODIN_VARIATIONAL_DOT_H

/**
 * @file
 * @brief Dot product (inner product) operations for functions and shape functions.
 *
 * This file provides dot product operations for vector and matrix-valued functions:
 * - Vector dot product: @f$ \mathbf{u} \cdot \mathbf{v} = \sum_i u_i v_i @f$
 * - Matrix Frobenius inner product: @f$ A : B = \sum_{ij} A_{ij} B_{ij} @f$
 */

#include <Eigen/Core>

#include "Rodin/Types.h"
#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Traits.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"
#include "RealFunction.h"

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived>
  struct Traits<
    Variational::Dot<
      Variational::FunctionBase<LHSDerived>,
      Variational::FunctionBase<RHSDerived>>>
  {
    using LHSType = Variational::FunctionBase<LHSDerived>;

    using RHSType = Variational::FunctionBase<RHSDerived>;

    using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

    using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

    using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

    using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

    using RangeType = ScalarType;
  };

  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Dot<
      Variational::FunctionBase<LHSDerived>,
      Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using LHSType = Variational::FunctionBase<LHSDerived>;

    using RHSType = Variational::ShapeFunctionBase<RHSDerived, FESType, Space>;

    using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

    using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

    using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

    using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

    using RangeType = ScalarType;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  struct Traits<
    Variational::Dot<
      Variational::ShapeFunctionBase<LHSDerived, TrialFES, Variational::TrialSpace>,
      Variational::ShapeFunctionBase<RHSDerived, TestFES, Variational::TestSpace>>>
  {
    using LHSType = Variational::ShapeFunctionBase<LHSDerived, TrialFES, Variational::TrialSpace>;

    using RHSType = Variational::ShapeFunctionBase<RHSDerived, TestFES, Variational::TestSpace>;

    using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

    using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

    using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

    using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

    using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

    using RangeType = ScalarType;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup DotSpecializations Dot Template Specializations
   * @brief Template specializations of the Dot class.
   *
   * Provides dot product (inner product) operations for:
   * - Function × Function → Real-valued function
   * - Function × ShapeFunction → Real-valued shape function  
   * - ShapeFunction × Function → Real-valued shape function
   * - ShapeFunction × ShapeFunction → Bilinear form entry
   *
   * @see Dot
   */

  /**
   * @ingroup DotSpecializations
   * @brief Dot product between two functions.
   *
   * Represents the inner product:
   * @f[
   *    (f \cdot g)(x) = f(x) \cdot g(x)
   * @f]
   * For vectors: @f$ \mathbf{u} \cdot \mathbf{v} = \sum_i u_i v_i @f$
   *
   * For matrices: @f$ A : B = \sum_{ij} A_{ij} B_{ij} @f$ (Frobenius inner product)
   *
   * @tparam LHSDerived Type of the first function
   * @tparam RHSDerived Type of the second function
   */
  template <class LHSDerived, class RHSDerived>
  class Dot<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>> final
    : public RealFunctionBase<Dot<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>>
  {
    public:
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

      using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

      using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

      using RangeType = ScalarType;

      using Parent = RealFunctionBase<Dot<LHSType, RHSType>>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      Dot(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Dot(Dot&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      template <class ... Args>
      Dot& traceOf(const Args& ... args)
      {
        m_lhs->traceOf(args...);
        m_rhs->traceOf(args...);
        return *this;
      }

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return Math::dot(this->object(getLHS().getValue(p)), this->object(getRHS().getValue(p)));
      }

      Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<FunctionBase<LHSDerived>> m_lhs;
      std::unique_ptr<FunctionBase<RHSDerived>> m_rhs;
  };

  /**
   * @brief Deduction guide for Dot product of two functions.
   */
  template <class LHSDerived, class RHSDerived>
  Dot(const FunctionBase<LHSDerived>&, const FunctionBase<RHSDerived>&)
    -> Dot<FunctionBase<LHSDerived>, FunctionBase<RHSDerived>>;

  /**
   * @ingroup DotSpecializations
   * @brief Dot product between a function and a shape function.
   *
   * Represents the mathematical expression:
   * @f[
   *   (\Lambda \cdot u)(x) = \Lambda(x) \cdot u(x)
   * @f]
   * with @f$ \Lambda @f$ a function and @f$ u @f$ a shape function.
   * For matrices: @f$ A : B @f$ denotes the Frobenius inner product.
   *
   * @tparam LHSDerived Type of the function
   * @tparam RHSDerived Type of the shape function
   * @tparam FES Finite element space type
   * @tparam Space Shape function space (Trial or Test)
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = ShapeFunctionBase<RHSDerived, FESType, Space>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

      using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

      using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

      using RangeType = ScalarType;

      using Parent = ShapeFunctionBase<Dot<LHSType, RHSType>, FESType, Space>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      Dot(const LHSType& lhs, const RHSType& rhs)
        : Parent(rhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      constexpr
      LHSType& getLHS()
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      RHSType& getRHS()
      {
        assert(m_rhs);
        return *m_rhs;
      }

      constexpr
      const auto& getLeaf() const
      {
        return getRHS().getLeaf();
      }

      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getRHS().getDOFs(element);
      }

      const FESType& getFiniteElementSpace() const
      {
        return getRHS().getFiniteElementSpace();
      }

      const Geometry::Point& getPoint() const
      {
        return getRHS().getPoint();
      }

      Dot& setPoint(const Geometry::Point& p)
      {
        m_rhs->setPoint(p);
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        const auto& p = getRHS().getPoint();
        return Math::dot(this->object(getLHS().getValue(p)), this->object(getRHS().getBasis(local)));
      }

      Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief Deduction guide for Dot product of function and shape function.
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Dot(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, Space>&)
    -> Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>;

  /**
   * @ingroup DotSpecializations
   * @brief Dot product between a shape function and a function.
   *
   * Represents the mathematical expression:
   * @f[
   *    (u \cdot \Lambda)(x) = u(x) \cdot \Lambda(x)
   * @f]
   * with @f$ u @f$ a shape function and @f$ \Lambda @f$ a function.
   *
   * @tparam LHSDerived Type of the shape function
   * @tparam RHSDerived Type of the function
   * @tparam FES Finite element space type
   * @tparam Space Shape function space (Trial or Test)
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Dot<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>> final
    : public ShapeFunctionBase<Dot<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>, FES, Space>
  {
    public:
      using FESType = FES;
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

      using LHSType = ShapeFunctionBase<LHSDerived, FESType, Space>;

      using RHSType = FunctionBase<RHSDerived>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

      using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

      using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

      using RangeType = ScalarType;

      using Parent = ShapeFunctionBase<Dot<LHSType, RHSType>, FESType, Space>;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      Dot(const LHSType& lhs, const RHSType& rhs)
        : Parent(lhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      constexpr
      const auto& getLeaf() const
      {
        return getLHS().getLeaf();
      }

      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getLHS().getDOFs(element);
      }

      const FESType& getFiniteElementSpace() const
      {
        return getLHS().getFiniteElementSpace();
      }

      const Geometry::Point& getPoint() const
      {
        return getLHS().getPoint();
      }

      Dot& setPoint(const Geometry::Point& p)
      {
        m_lhs->setPoint(p);
        return *this;
      }

      constexpr
      auto getBasis(size_t local) const
      {
        const auto& p = getLHS().getPoint();
        return Math::dot(this->object(getLHS().getBasis(local)), this->object(getRHS().getValue(p)));
      }

      Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief Deduction guide for Dot product of shape function and function.
   */
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Dot(const ShapeFunctionBase<LHSDerived, FES, Space>&, const FunctionBase<RHSDerived>&)
    -> Dot<ShapeFunctionBase<LHSDerived, FES, Space>, FunctionBase<RHSDerived>>;

  /**
   * @ingroup DotSpecializations
   * @brief Dot product of trial and test shape functions for bilinear forms.
   *
   * Represents the inner product in a bilinear form:
   * @f[
   *    (u, v) \mapsto u(x) : v(x)
   * @f]
   * Used in assembling bilinear forms @f$ a(u, v) @f$.
   *
   * @tparam LHSDerived Type of the trial shape function
   * @tparam TrialFES Trial finite element space type
   * @tparam RHSDerived Type of the test shape function
   * @tparam TestFES Test finite element space type
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Dot<
    ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
    ShapeFunctionBase<RHSDerived, TestFES, TestSpace>> final
    : public FormLanguage::Base
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;

      using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;

      using LHSScalarType = typename FormLanguage::Traits<LHSRangeType>::ScalarType;

      using RHSScalarType = typename FormLanguage::Traits<RHSRangeType>::ScalarType;

      using ScalarType = typename FormLanguage::Mult<LHSScalarType, RHSScalarType>::Type;

      using RangeType = ScalarType;

      using Parent = FormLanguage::Base;

      static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);

      constexpr
      Dot(const LHSType& lhs, const RHSType& rhs)
        : m_trial(lhs.copy()), m_test(rhs.copy())
      {}

      constexpr
      Dot(const Dot& other)
        : Base(other),
          m_trial(other.m_trial->copy()), m_test(other.m_test->copy())
      {}

      constexpr
      Dot(Dot&& other)
        : Base(std::move(other)),
          m_trial(std::move(other.m_trial)), m_test(std::move(other.m_test))
      {}

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_trial);
        return *m_trial;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_test);
        return *m_test;
      }

      const Geometry::Point& getPoint() const
      {
        return m_trial->getPoint();
      }

      Dot& setPoint(const Geometry::Point& p)
      {
        m_trial->setPoint(p);
        m_test->setPoint(p);
        return *this;
      }

      constexpr
      auto operator()(size_t tr, size_t te)
      {
        return Math::dot(this->object(getLHS().getBasis(tr)), this->object(getRHS().getBasis(te)));
      }

      Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_trial;
      std::unique_ptr<RHSType> m_test;
  };

  /**
   * @brief Deduction guide for Dot product of trial and test shape functions.
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Dot(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&,
      const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
  -> Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>;

  /**
   * @ingroup DotSpecializations
   * @brief Dot product of a potential with a test shape function.
   *
   * Used for integral operators involving potentials in boundary element methods.
   *
   * @tparam KernelType Type of the kernel function
   * @tparam LHSDerived Type of the potential's shape function
   * @tparam TrialFES Trial finite element space type
   * @tparam RHSDerived Type of the test shape function
   * @tparam TestFES Test finite element space type
   */
  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Dot<
      Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>
    : public FormLanguage::Base
  {
    public:
      using LHSType = Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>;

      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;

      using Parent = FormLanguage::Base;

      constexpr
      Dot(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  /**
   * @brief Deduction guide for Dot product of potential and test shape function.
   */
  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Dot(const Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>&,
      const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
  -> Dot<
      Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>;

}

#endif
