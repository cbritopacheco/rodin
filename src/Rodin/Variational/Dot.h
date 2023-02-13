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
#include "Rodin/Math/DenseMatrix.h"

#include "ForwardDecls.h"

#include "Function.h"
#include "ShapeFunction.h"
#include "ScalarFunction.h"

#include "Exceptions.h"

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
  template <>
  class Dot<FunctionBase, FunctionBase> : public ScalarFunctionBase
  {
    public:
      using LHS = FunctionBase;
      using RHS = FunctionBase;

      Dot(const FunctionBase& a, const FunctionBase& b);

      Dot(const Dot& other);

      Dot(Dot&& other);

      Dot& traceOf(Geometry::Attribute attrs) override;

      virtual LHS& getLHS()
      {
        return *m_a;
      }

      virtual RHS& getRHS()
      {
        return *m_b;
      }

      virtual const LHS& getLHS() const
      {
        return *m_a;
      }

      virtual const RHS& getRHS() const
      {
        return *m_b;
      }

      virtual FunctionValue getValue(const Geometry::Point& p) const override
      {
        assert(m_a->getRangeShape() == m_b->getRangeShape());
        assert(m_a->getRangeType() == m_b->getRangeType());
        switch (m_a->getRangeType())
        {
          case RangeType::Scalar:
          {
            return m_a->getValue(p).scalar() * m_b->getValue(p).scalar();
          }
          case RangeType::Vector:
          {
            Math::Vector a = m_a->getValue(p).vector();
            Math::Vector b = m_b->getValue(p).vector();
            return a.dot(b);
          }
          case RangeType::Matrix:
          {
            Math::Matrix a = m_a->getValue(p).matrix();
            Math::Matrix b = m_b->getValue(p).matrix();
            return (a * b.transpose()).trace();
          }
          default:
          {
            assert(false);
            return 0.0;
          }
        }
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
    private:
      std::unique_ptr<FunctionBase> m_a, m_b;
  };
  Dot(const FunctionBase&, const FunctionBase&) -> Dot<FunctionBase, FunctionBase>;

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
  template <ShapeFunctionSpaceType Space>
  class Dot<FunctionBase, ShapeFunctionBase<Space>> : public ShapeFunctionBase<Space>
  {
    public:
      using Parent = ShapeFunctionBase<Space>;
      using LHS = FunctionBase;
      using RHS = ShapeFunctionBase<Space>;

      Dot(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
        : m_lhs(lhs.copy()), m_rhs(rhs.copy())
      {
        if (lhs.getRangeShape() != rhs.getRangeShape())
          RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
      }

      Dot(const Dot& other)
        :  Parent(other),
          m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
      {}

      Dot(Dot&& other)
        :  Parent(std::move(other)),
          m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
      {}

      virtual LHS& getLHS()
      {
        return *m_lhs;
      }

      virtual RHS& getRHS()
      {
        return *m_rhs;
      }

      virtual const LHS& getLHS() const
      {
        return *m_lhs;
      }

      virtual const RHS& getRHS() const
      {
        return *m_rhs;
      }

      const ShapeFunctionBase<Space>& getLeaf() const override
      {
        return getRHS().getLeaf();
      }

      int getRows() const override
      {
        return 1;
      }

      int getDOFs(const Geometry::Simplex& element) const override
      {
        return getRHS().getDOFs(element);
      }

      int getColumns() const override
      {
        return 1;
      }

      FiniteElementSpaceBase& getFiniteElementSpace() override
      {
        return m_rhs->getFiniteElementSpace();
      }

      const FiniteElementSpaceBase& getFiniteElementSpace() const override
      {
        return getRHS().getFiniteElementSpace();
      }

      virtual TensorBasis getOperator(
          ShapeComputator& compute, const Geometry::Point& p) const override
      {
        const size_t dofs = getRHS().getDOFs(p.getSimplex());
        const TensorBasis basis = getRHS().getOperator(compute, p);
        switch (getLHS().getRangeType())
        {
          case RangeType::Scalar:
          {
            return getLHS().getValue(p) * basis;
          }
          case RangeType::Vector:
          {
            TensorBasis res(static_cast<int>(dofs), 1, 1);
            const Math::Vector vec = getLHS().getValue(p);
            for (size_t i = 0; i < dofs; i++)
            {
              const auto dof = Math::slice(basis, 0, i);
              res(static_cast<int>(i), 0, 0) = (dof * vec.transpose()).trace();
            }
            return res;
          }
          case RangeType::Matrix:
          {
            TensorBasis res(static_cast<int>(dofs), 1, 1);
            const Math::Matrix mat = getLHS().getValue(p);
            for (size_t i = 0; i < dofs; i++)
            {
              const auto dof = Math::slice(basis, 0, i);
              res(static_cast<int>(i), 0, 0) = (dof * mat.transpose()).trace();
            }
            return res;
          }
        }
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<FunctionBase> m_lhs;
      std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
  };
  template <ShapeFunctionSpaceType Space>
  Dot(const FunctionBase&, const ShapeFunctionBase<Space>&)
    -> Dot<FunctionBase, ShapeFunctionBase<Space>>;

  /**
   * @ingroup DotSpecializations
   */
  template <>
  class Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
    : public FormLanguage::Base
  {
    public:
      using Parent = FormLanguage::Base;
      using LHS = ShapeFunctionBase<TrialSpace>;
      using RHS = ShapeFunctionBase<TestSpace>;

      Dot(const ShapeFunctionBase<TrialSpace>& lhs, const ShapeFunctionBase<TestSpace>& rhs)
        : m_trial(lhs.copy()), m_test(rhs.copy())
      {
        if (lhs.getRangeShape() != rhs.getRangeShape())
          RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
      }

      Dot(const Dot& other)
        :  Base(other),
          m_trial(other.m_trial->copy()), m_test(other.m_test->copy())
      {}

      Dot(Dot&& other)
        :  Base(std::move(other)),
          m_trial(std::move(other.m_trial)), m_test(std::move(other.m_test))
      {}

      inline virtual LHS& getLHS()
      {
        assert(m_trial);
        return *m_trial;
      }

      inline virtual RHS& getRHS()
      {
        assert(m_test);
        return *m_test;
      }

      inline virtual const LHS& getLHS() const
      {
        assert(m_trial);
        return *m_trial;
      }

      inline virtual const RHS& getRHS() const
      {
        assert(m_test);
        return *m_test;
      }

      inline virtual Math::Matrix getMatrix(
          ShapeComputator& compute, const Geometry::Point& p) const
      {
        assert(m_trial->getRangeShape() == m_test->getRangeShape());
        const TensorBasis trial = m_trial->getOperator(compute, p);
        const TensorBasis test = m_test->getOperator(compute, p);
        assert(trial.rank() == test.rank());
        Math::Matrix res(test.dimension(0), trial.dimension(0));
        for (Math::Matrix::Index i = 0; i < res.rows(); i++)
        {
          const auto testdof = Math::slice(test, 0, i);
          for (Math::Matrix::Index j = 0; j < res.cols(); j++)
          {
            const auto trialdof = Math::slice(trial, 0, j);
            res(i, j) = (testdof * trialdof.transpose()).trace();
          }
        }
        return res;
      }

      inline virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }

    private:
      std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_trial;
      std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_test;
  };
  Dot(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
    -> Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;

  /* ||-- OPTIMIZATIONS -----------------------------------------------------
   * Dot<FunctionBase, ShapeFunctionBase<Space>>
   * ---------------------------------------------------------------------->>
   */

  /**
   * @ingroup DotSpecializations
   * @f[
   *   f \cdot u
   * @f]
   * where @f$ f @f$ is a function (scalar or vector valued).
   */
  template <class FES, ShapeFunctionSpaceType Space>
  class Dot<FunctionBase, ShapeFunction<FES, Space>>
    : public Dot<FunctionBase, ShapeFunctionBase<Space>>
  {
    public:
      using Parent = Dot<FunctionBase, ShapeFunctionBase<Space>>;
      using LHS = FunctionBase;
      using RHS = ShapeFunction<FES, Space>;

      constexpr
      Dot(const FunctionBase& f, const ShapeFunction<FES, Space>& u)
        : Parent(f, u)
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(other)
      {}

      virtual LHS& getLHS() override
      {
        return static_cast<LHS&>(Parent::getLHS());
      }

      virtual RHS& getRHS() override
      {
        return static_cast<RHS&>(Parent::getRHS());
      }

      virtual const LHS& getLHS() const override
      {
        return static_cast<const LHS&>(Parent::getLHS());
      }

      virtual const RHS& getRHS() const override
      {
        return static_cast<const RHS&>(Parent::getRHS());
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
  };
  template <class FES, ShapeFunctionSpaceType Space>
  Dot(const FunctionBase&, const ShapeFunction<FES, Space>&)
    -> Dot<FunctionBase, ShapeFunction<FES, Space>>;

  /* <<-- OPTIMIZATIONS -----------------------------------------------------
   * Dot<FunctionBase, ShapeFunctionBase<TestSpace>>
   * ----------------------------------------------------------------------||
   */

  /* ||-- OPTIMIZATIONS -----------------------------------------------------
   * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
   * ---------------------------------------------------------------------->>
   */

  /**
   * @ingroup DotSpecializations
   *
   * @f[
   *   (f u) \cdot v
   * @f]
   * where @f$ f @f$ is a function (scalar or vector valued).
   */
  template <class FES>
  class Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>
    : public Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  {
    public:
      using Parent = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
      using LHS = Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>;
      using RHS = ShapeFunction<FES, TestSpace>;

      constexpr
      Dot(const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>& fu,
          const ShapeFunction<FES, TestSpace>& v)
        : Parent(fu, v)
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(other)
      {}

      virtual LHS& getLHS() override
      {
        return static_cast<LHS&>(Parent::getLHS());
      }

      virtual RHS& getRHS() override
      {
        return static_cast<RHS&>(Parent::getRHS());
      }

      virtual const LHS& getLHS() const override
      {
        return static_cast<const LHS&>(Parent::getLHS());
      }

      virtual const RHS& getRHS() const override
      {
        return static_cast<const RHS&>(Parent::getRHS());
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
  };
  template <class FES>
  Dot(const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&, const ShapeFunction<FES, TestSpace>&)
    -> Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;

  /**
   * @ingroup DotSpecializations
   *
   * @f[
   *   (f \nabla u) \cdot \nabla v
   * @f]
   * where @f$ f @f$ is a function (scalar or matrix valued).
   */
  template <class FES>
  class Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>
    : public Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  {
    public:
      using Parent = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
      using LHS = Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>;
      using RHS = Grad<ShapeFunction<FES, TestSpace>>;

      constexpr
      Dot(const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>& fgu,
          const Grad<ShapeFunction<FES, TestSpace>>& gv)
        : Parent(fgu, gv)
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(other)
      {}

      virtual LHS& getLHS() override
      {
        return static_cast<LHS&>(Parent::getLHS());
      }

      virtual RHS& getRHS() override
      {
        return static_cast<RHS&>(Parent::getRHS());
      }

      virtual const LHS& getLHS() const override
      {
        return static_cast<const LHS&>(Parent::getLHS());
      }

      virtual const RHS& getRHS() const override
      {
        return static_cast<const RHS&>(Parent::getRHS());
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
  };
  template <class FES>
  Dot(const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&, const Grad<ShapeFunction<FES, TestSpace>>&)
    -> Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;

  /**
   * @ingroup DotSpecializations
   *
   * Represents the following expression:
   * @f[
   *   \nabla u \cdot \nabla v
   * @f]
   */
  template <class FES>
  class Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>
    : public Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  {
    public:
      using Parent = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
      using LHS = Grad<ShapeFunction<FES, TrialSpace>>;
      using RHS = Grad<ShapeFunction<FES, TestSpace>>;

      constexpr
      Dot(const Grad<ShapeFunction<FES, TrialSpace>>& nu, const Grad<ShapeFunction<FES, TestSpace>>& nv)
        : Parent(nu, nv)
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other))
      {}

      virtual LHS& getLHS() override
      {
        return static_cast<LHS&>(Parent::getLHS());
      }

      virtual RHS& getRHS() override
      {
        return static_cast<RHS&>(Parent::getRHS());
      }

      virtual const LHS& getLHS() const override
      {
        return static_cast<const LHS&>(Parent::getLHS());
      }

      virtual const RHS& getRHS() const override
      {
        return static_cast<const RHS&>(Parent::getRHS());
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
  };
  template <class FES>
  Dot(const Grad<ShapeFunction<FES, TrialSpace>>&, const Grad<ShapeFunction<FES, TestSpace>>&)
    -> Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>;

  /**
   * @ingroup DotSpecializations
   *
   * Represents the following expression:
   * @f[
   *   \mathbf{J} u \cdot \mathbf{J} v
   * @f]
   */
  template <class FES>
  class Dot<Jacobian<ShapeFunction<FES, TrialSpace>>, Jacobian<ShapeFunction<FES, TestSpace>>>
    : public Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  {
    public:
      using Parent = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
      using LHS = Jacobian<ShapeFunction<FES, TrialSpace>>;
      using RHS = Jacobian<ShapeFunction<FES, TestSpace>>;

      constexpr
      Dot(const Jacobian<ShapeFunction<FES, TrialSpace>>& nu, const Jacobian<ShapeFunction<FES, TestSpace>>& nv)
        : Parent(nu, nv)
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(std::move(other))
      {}

      virtual LHS& getLHS() override
      {
        return static_cast<LHS&>(Parent::getLHS());
      }

      virtual RHS& getRHS() override
      {
        return static_cast<RHS&>(Parent::getRHS());
      }

      virtual const LHS& getLHS() const override
      {
        return static_cast<const LHS&>(Parent::getLHS());
      }

      virtual const RHS& getRHS() const override
      {
        return static_cast<const RHS&>(Parent::getRHS());
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
  };
  template <class FES>
  Dot(const Jacobian<ShapeFunction<FES, TrialSpace>>&, const Jacobian<ShapeFunction<FES, TestSpace>>&)
    -> Dot<Jacobian<ShapeFunction<FES, TrialSpace>>, Jacobian<ShapeFunction<FES, TestSpace>>>;

  /**
   * @ingroup DotSpecializations
   *
   * @f[
   *   (f \mathbf{J} u) \cdot \mathbf{J} v
   * @f]
   * where @f$ f @f$ is a function (scalar or matrix valued).
   */
  template <class FES>
  class Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>
    : public Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
  {
    public:
      using Parent = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
      using LHS = Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>;
      using RHS = Jacobian<ShapeFunction<FES, TestSpace>>;

      constexpr
      Dot(const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>& fgu,
          const Jacobian<ShapeFunction<FES, TestSpace>>& gv)
        : Parent(fgu, gv)
      {}

      constexpr
      Dot(const Dot& other)
        : Parent(other)
      {}

      constexpr
      Dot(Dot&& other)
        : Parent(other)
      {}

      virtual LHS& getLHS() override
      {
        return static_cast<LHS&>(Parent::getLHS());
      }

      virtual RHS& getRHS() override
      {
        return static_cast<RHS&>(Parent::getRHS());
      }

      virtual const LHS& getLHS() const override
      {
        return static_cast<const LHS&>(Parent::getLHS());
      }

      virtual const RHS& getRHS() const override
      {
        return static_cast<const RHS&>(Parent::getRHS());
      }

      virtual Dot* copy() const noexcept override
      {
        return new Dot(*this);
      }
  };
  template <class FES>
  Dot(const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>&, const Jacobian<ShapeFunction<FES, TestSpace>>&)
    -> Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>;

  /* <<-- OPTIMIZATIONS -----------------------------------------------------
   * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
   * ----------------------------------------------------------------------||
   */
}

#endif
