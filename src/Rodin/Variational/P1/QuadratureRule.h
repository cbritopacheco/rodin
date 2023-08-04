#ifndef RODIN_VARIATIONAL_P1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_P1_QUADRATURERULE_H

#include "Rodin/Variational/QuadratureRule.h"

#include "P1.h"
#include "P1Element.h"

namespace Rodin::Variational
{
  // /**
  //  * @ingroup QuadratureRuleSpecializations
  //  *
  //  * @f[
  //  * \int A u \cdot v \ dx
  //  * @f]
  //  */
  // template <class FunctionDerived, class LHSDerived, class RHSDerived, class Context>
  // class QuadratureRule<
  //   Dot<
  //     Mult<
  //       FunctionBase<FunctionDerived>,
  //       ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>,
  //     ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>
  //   : public BilinearFormIntegratorBase
  // {
  //   public:
  //     using Parent = BilinearFormIntegratorBase;
  //     using LHS = Mult<
  //       FunctionBase<FunctionDerived>,
  //       ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>;
  //     using RHS = ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>;
  //     using Integrand = Dot<LHS, RHS>;
  // };

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int \nabla u \cdot \nabla v \ dx
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;
      using LHS = ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>, TrialSpace>;
      using RHS = ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>;
      using Integrand = Dot<LHS, RHS>;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : BilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      void assemble(const Geometry::Polytope& polytope) override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& trial = integrand.getLHS();
        const auto& trans = polytope.getTransformation();
        const auto& fes = trial.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        assert(dofs == trial.getDOFs(polytope));
        const QF::QF1P1 qf(polytope.getGeometry());
        auto& res = getMatrix();
        res = Math::Matrix::Zero(dofs, dofs);
        m_gradient.resize(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
          const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
          const auto jacInvT = p.getJacobianInverse().transpose();
          for (size_t local = 0; local < dofs; local++)
            m_gradient[local] = jacInvT * fe.getGradient(local)(rc);
          for (size_t i = 0; i < dofs; i++)
            res(i, i) += qf.getWeight(k) * p.getDistortion() * m_gradient[i].squaredNorm();
          for (size_t i = 0; i < dofs; i++)
            for (size_t j = 0; j < i; j++)
              res(i, j) += qf.getWeight(k) * p.getDistortion() * m_gradient[i].dot(m_gradient[j]);
        }
        res.template triangularView<Eigen::Upper>() = res.transpose();
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
      std::vector<Math::Vector> m_gradient;
  };

  template <class LHSDerived, class RHSDerived, class Context>
  QuadratureRule(const Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>>&)
    -> QuadratureRule<Dot<
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>, TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int A \nabla u \cdot \nabla v \ dx
   * @f]
   */
  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<
  Dot<
    ShapeFunctionBase<
      Mult<
        FunctionBase<LHSFunctionDerived>,
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>, TrialSpace>>,
      P1<Scalar, Context>, TrialSpace>,
    ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;

      using LHS =
        ShapeFunctionBase<Mult<FunctionBase<LHSFunctionDerived>,
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>,
        TrialSpace>>, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>,
        TrialSpace>;

      using RHS = ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>;

      using Integrand = Dot<LHS, RHS>;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : BilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      void assemble(const Geometry::Polytope& polytope) override
      {
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& integrand = getIntegrand();
        const auto& trial = integrand.getLHS();
        const auto& f = static_cast<const Mult<FunctionBase<LHSFunctionDerived>,
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>,
        TrialSpace>>, P1<Scalar, Context>, TrialSpace>>&>(integrand.getLHS()).getLHS();
        const auto& trans = polytope.getTransformation();
        const auto& fes = trial.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        assert(dofs == trial.getDOFs(polytope));
        const QF::QF1P1 qf(polytope.getGeometry());
        auto& res = getMatrix();
        res = Math::Matrix::Zero(dofs, dofs);
        m_gradient.resize(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
          const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
          const auto jacInvT = p.getJacobianInverse().transpose();
          for (size_t local = 0; local < dofs; local++)
            m_gradient[local] = jacInvT * fe.getGradient(local)(rc);
          const auto fv = f.getValue(p);
          for (size_t i = 0; i < dofs; i++)
            res(i, i) += qf.getWeight(k) * p.getDistortion() * (fv * m_gradient[i]).dot(m_gradient[i]);
          for (size_t i = 0; i < dofs; i++)
            for (size_t j = 0; j < i; j++)
              res(i, j) += qf.getWeight(k) * p.getDistortion() * (fv * m_gradient[i]).dot(m_gradient[j]);
        }
        res.template triangularView<Eigen::Upper>() = res.transpose();
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
      std::vector<Math::Vector> m_gradient;
  };

  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Context>
  QuadratureRule(const Dot<
        ShapeFunctionBase<Mult<FunctionBase<LHSFunctionDerived>,
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>,
        TrialSpace>>, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>,
        TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>>&)
    -> QuadratureRule<Dot<
        ShapeFunctionBase<Mult<FunctionBase<LHSFunctionDerived>,
        ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>,
        TrialSpace>>, P1<Scalar, Context>, TrialSpace>>, P1<Scalar, Context>,
        TrialSpace>,
        ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>, P1<Scalar, Context>, TestSpace>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int f v \ dx
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<
      ShapeFunctionBase<Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>>,
      P1<Scalar, Context>, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using LHS = FunctionBase<LHSDerived>;

      using RHS = ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>;

      using FES = P1<Scalar, Context>;

      using Integrand = ShapeFunctionBase<Dot<LHS, RHS>, P1<Scalar, Context>, TestSpace>;

      using Parent = LinearFormIntegratorBase;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      void assemble(const Geometry::Polytope& polytope) final override
      {
        using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;
        using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;
        static_assert(std::is_same_v<LHSRange, RHSRange>);
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& integrand = static_cast<const Dot<LHS, RHS>&>(getIntegrand());
        const auto& f = integrand.getLHS();
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        assert(integrand.getRangeType() == RangeType::Scalar);
        const QF::QF1P1 qf(polytope.getGeometry());
        auto& res = getVector();
        res = Math::Vector::Zero(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
          const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
          if constexpr (std::is_same_v<Scalar, LHSRange>)
          {
            for (size_t local = 0; local < dofs; local++)
              res.coeffRef(local) += qf.getWeight(k) * p.getDistortion() * f(p) * fe.getBasis(local)(rc);
          }
          else if constexpr (std::is_same_v<Math::Vector, LHSRange>)
          {
            for (size_t local = 0; local < dofs; local++)
              res.coeffRef(local) += qf.getWeight(k) * p.getDistortion() * f(p).dot(fe.getBasis(local)(rc));
          }
          else if constexpr (std::is_same_v<Math::Matrix, LHSRange>)
          {
            for (size_t local = 0; local < dofs; local++)
              res.coeffRef(local) +=
                qf.getWeight(k) * p.getDistortion() * (f(p).array() * fe.getBasis(local)(rc).array()).rowwise().sum().colwise().sum().value();
          }
          else
          {
            assert(false);
            res.setConstant(NAN);
          }
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class LHSDerived, class RHSDerived, class Context>
  QuadratureRule(const ShapeFunctionBase<
    Dot<
      FunctionBase<LHSDerived>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>>,
    P1<Scalar, Context>, TestSpace>&)
  -> QuadratureRule<ShapeFunctionBase<
      Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>>,
      P1<Scalar, Context>, TestSpace>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   *  \int f v \ dx
   * @f]
   */
  template <class NestedDerived, class Context>
  class QuadratureRule<
    ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using Nested = ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>;

      using FES = P1<Scalar, Context>;

      using Integrand = ShapeFunctionBase<Nested, P1<Scalar, Context>, TestSpace>;

      using Parent = LinearFormIntegratorBase;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : Parent(integrand.getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy())
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand))
      {}

      inline
      constexpr
      const Integrand& getIntegrand() const
      {
        assert(m_integrand);
        return *m_integrand;
      }

      void assemble(const Geometry::Polytope& polytope) final override
      {
        using Range = typename FormLanguage::Traits<Integrand>::RangeType;
        static_assert(std::is_same_v<Range, Scalar>);
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& integrand = static_cast<const Nested&>(getIntegrand());
        const auto& fes = integrand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        assert(integrand.getRangeType() == RangeType::Scalar);
        const QF::QF1P1 qf(polytope.getGeometry());
        auto& res = getVector();
        res = Math::Vector::Zero(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
          const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
          for (size_t local = 0; local < dofs; local++)
            res.coeffRef(local) += qf.getWeight(k) * p.getDistortion() * fe.getBasis(local)(rc);
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class NestedDerived, class Context>
  QuadratureRule(const ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>&)
  -> QuadratureRule<ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>, P1<Scalar, Context>, TestSpace>>;
}

#endif

