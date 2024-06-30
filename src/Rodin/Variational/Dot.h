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

namespace Rodin::FormLanguage
{
  template <class LHSDerived, class RHSDerived>
  struct Traits<Variational::Dot<
          Variational::FunctionBase<LHSDerived>,
          Variational::FunctionBase<RHSDerived>>>
  {
    using LHSType = Variational::FunctionBase<LHSDerived>;
    using RHSType = Variational::FunctionBase<RHSDerived>;
  };

  template <class LHSDerived, class RHSDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Dot<
      Variational::FunctionBase<LHSDerived>,
      Variational::ShapeFunctionBase<RHSDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
  };
}

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
      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = FunctionBase<RHSDerived>;

      using Parent = ScalarFunctionBase<Dot<LHSType, RHSType>>;

      Dot(const LHSType& lhs, const RHSType& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

      Dot(const Dot& other)
        : Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

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
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        assert(getLHS().getRangeShape() == getRHS().getRangeShape());
        using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
        using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
        static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);
        const auto& lhs = this->object(getLHS().getValue(p));
        const auto& rhs = this->object(getRHS().getValue(p));
        if constexpr (std::is_same_v<LHSRangeType, Scalar>)
        {
          return lhs * rhs;
        } else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
        {
          return lhs.dot(rhs);
        }
        else if constexpr (std::is_same_v<RHSRangeType, Math::Matrix<Scalar>>)
        {
          return (lhs.array() * rhs.array()).rowwise().sum().colwise().sum().value();
        }
        else
        {
          assert(false);
          return void();
        }
      }

      inline Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<FunctionBase<LHSDerived>> m_lhs;
      std::unique_ptr<FunctionBase<RHSDerived>> m_rhs;
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
  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  class Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>> final
    : public ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>, FES, Space>
  {
    public:
      using FESType = FES;

      using LHSType = FunctionBase<LHSDerived>;

      using RHSType = ShapeFunctionBase<RHSDerived, FESType, Space>;

      using Parent = ShapeFunctionBase<Dot<LHSType, RHSType>, FESType, Space>;

      constexpr
      Dot(const LHSType& lhs, const RHSType& rhs)
        : Parent(rhs.getFiniteElementSpace()),
          m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

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

      inline
      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
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
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getRHS().getDOFs(element);
      }

      inline
      constexpr
      auto getTensorBasis(const Geometry::Point& p) const
      {
        assert(m_lhs->getRangeShape() == m_rhs->getRangeShape());
        using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
        using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
        static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);
        const auto& lhs = this->object(getLHS().getValue(p));
        const auto rhs = getRHS().getTensorBasis(p);
        if constexpr (std::is_same_v<LHSRangeType, Scalar>)
        {
          return lhs * rhs;
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
        {
          return TensorBasis(rhs.getDOFs(),
              [&](size_t i){ return lhs.dot(rhs(i)); });
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<Scalar>>)
        {
          return TensorBasis(rhs.getDOFs(),
              [&](size_t i){ return (lhs(i).array() * rhs(i).array()).rowwise().sum().colwise().sum().value(); });
        }
        else
        {
          assert(false);
          return void();
        }
      }

      inline
      Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class LHSDerived, class RHSDerived, class FES, ShapeFunctionSpaceType Space>
  Dot(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, Space>&)
    -> Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, Space>>;

  /**
   * @ingroup DotSpecializations
   * @brief Represents the dot product of trial and test operators.
   *
   * Constructs an instance representing the following expression:
   * @f[
   *    A(u) : B(v)
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>> final
    : public FormLanguage::Base
  {
    public:
      using LHSType = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;
      using RHSType = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;
      using Parent = FormLanguage::Base;

      constexpr
      Dot(const LHSType& lhs, const RHSType& rhs)
        : m_trial(lhs.copy()), m_test(rhs.copy())
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

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

      inline
      constexpr
      const LHSType& getLHS() const
      {
        assert(m_trial);
        return *m_trial;
      }

      inline
      constexpr
      const RHSType& getRHS() const
      {
        assert(m_test);
        return *m_test;
      }

      /**
       * @brief Gets the element matrix at a point.
       *
       * Computes the @f$ m \times n @f$ element matrix @f$ M @f$ defined by:
       * @f[
       *    M = V U^T
       * @f]
       * where @f$ n @f$ is the number of trial degrees of freedom, and @f$ m
       * @f$ is the number of test degrees of freedom.
       */
      void assemble(const Geometry::Point& p)
      {
        assert(getLHS().getRangeShape() == getRHS().getRangeShape());
        using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
        using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
        static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);
        const auto& trial = getLHS().getTensorBasis(p);
        const auto& test = getRHS().getTensorBasis(p);
        m_matrix.resize(test.getDOFs(), trial.getDOFs());
        if constexpr (std::is_same_v<LHSRangeType, Scalar>)
        {
          for (size_t i = 0; i < test.getDOFs(); i++)
            for (size_t j = 0; j < trial.getDOFs(); j++)
              m_matrix(i, j) = test(i) * trial(j);
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Vector<Scalar>>)
        {
          for (size_t i = 0; i < test.getDOFs(); i++)
            for (size_t j = 0; j < trial.getDOFs(); j++)
              m_matrix(i, j) = test(i).dot(trial(j));
        }
        else if constexpr (std::is_same_v<LHSRangeType, Math::Matrix<Scalar>>)
        {
          for (size_t i = 0; i < test.getDOFs(); i++)
            for (size_t j = 0; j < trial.getDOFs(); j++)
              m_matrix(i, j) = (test(i).array() * trial(j).array()).rowwise().sum().colwise().sum().value();
        }
        else
        {
          assert(false);
          m_matrix.setConstant(NAN);
        }
      }

      const Math::Matrix<Scalar>& getMatrix() const
      {
        return m_matrix;
      }

      inline Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_trial;
      std::unique_ptr<RHSType> m_test;

      Math::Matrix<Scalar> m_matrix;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Dot(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&,
      const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
  -> Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>;

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
      {
        assert(lhs.getRangeShape() == rhs.getRangeShape());
      }

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

      inline
      constexpr
      const LHSType& getLHS() const
      {
        assert(m_lhs);
        return *m_lhs;
      }

      inline
      constexpr
      const RHSType& getRHS() const
      {
        assert(m_rhs);
        return *m_rhs;
      }

      inline
      Dot* copy() const noexcept final override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<LHSType> m_lhs;
      std::unique_ptr<RHSType> m_rhs;
  };

  template <class KernelType, class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Dot(const Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>&,
      const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
  -> Dot<
      Potential<KernelType, ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>>,
      ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>;

}

#endif
