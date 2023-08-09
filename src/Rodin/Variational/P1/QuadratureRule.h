#ifndef RODIN_VARIATIONAL_P1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_P1_QUADRATURERULE_H

#include "Rodin/Variational/QuadratureRule.h"

#include "P1.h"
#include "P1Element.h"

namespace Rodin::Variational
{
  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int (A u) \cdot v \ dx
   * @f]
   */
  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>,
            P1<Scalar, Context>, TrialSpace>>,
            P1<Scalar, Context>, TrialSpace>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>,
          P1<Scalar, Context>, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;

      using FES = P1<Scalar, Context>;

      using Coefficient =
        FunctionBase<CoefficientDerived>;

      using Multiplicand =
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, FES, TrialSpace>>>;

      using LHS =
        ShapeFunctionBase<
          Mult<
            FunctionBase<CoefficientDerived>,
            ShapeFunctionBase<ShapeFunction<LHSDerived, FES, TrialSpace>>>>;

      using RHS =
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, FES, TestSpace>>;

      using Integrand = Dot<LHS, RHS>;

      using CoefficientRange = typename FormLanguage::Traits<Coefficient>::RangeType;

      using MultiplicandRange = typename FormLanguage::Traits<Multiplicand>::RangeType;

      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      static_assert(std::is_same_v<LHSRange, RHSRange>);

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
        const auto& coeff = integrand.getLHS().getDerived().getLHS();
        const auto& multiplicand = integrand.getLHS();
        const auto& trans = polytope.getTransformation();
        const auto& fes = multiplicand.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        const QF::QF1P1 qf(polytope.getGeometry());
        auto& res = getMatrix();
        res = Math::Matrix::Zero(dofs, dofs);
        if constexpr (std::is_same_v<CoefficientRange, Scalar>)
        {
          static_assert(std::is_same_v<MultiplicandRange, RHSRange>);
          if constexpr (std::is_same_v<MultiplicandRange, Scalar>)
          {
            for (size_t k = 0; k < qf.getSize(); k++)
            {
              Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
              const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
              const Scalar c = coeff.getValue(p);
              for (size_t i = 0; i < dofs; i++)
              {
                const Scalar basis = fe.getBasis(i)(rc);
                res(i, i) += qf.getWeight(k) * p.getDistortion() * c * basis * basis;
              }
              for (size_t i = 0; i < dofs; i++)
                for (size_t j = 0; j < i; j++)
                  res(i, j) += qf.getWeight(k) * p.getDistortion() * c * fe.getBasis(i)(rc) * fe.getBasis(j)(rc);
            }
            res.template triangularView<Eigen::Upper>() = res.transpose();
          }
          else
          {
            assert(false); // Not handled yet
            res.setConstant(NAN);
            return;
          }
        }
        else
        {
          assert(false); // Not handled yet
          res.setConstant(NAN);
          return;
        }
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
      std::vector<Math::Vector> m_vvalues;
      std::vector<Math::Matrix> m_mvalues;
  };

  template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Context>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<CoefficientDerived>,
          ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>>,
      ShapeFunctionBase<
        ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int (A \nabla u) \cdot \nabla v \ dx
   * @f]
   */
  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<
  Dot<
    ShapeFunctionBase<
      Mult<
        FunctionBase<LHSFunctionDerived>,
        ShapeFunctionBase<
          Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>,
        P1<Scalar, Context>, TrialSpace>>,
      P1<Scalar, Context>, TrialSpace>,
    ShapeFunctionBase<
      Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>,
    P1<Scalar, Context>, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;

      using FES = P1<Scalar, Context>;

      using LHS =
        ShapeFunctionBase<
          Mult<
            FunctionBase<LHSFunctionDerived>,
            ShapeFunctionBase<
              Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>>>>;

      using RHS =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, FES, TestSpace>>>;

      using Integrand = Dot<LHS, RHS>;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : BilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_rgradient(other.m_rgradient),
          m_pgradient(other.m_pgradient)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_rgradient(std::move(other.m_rgradient)),
          m_pgradient(std::move(other.m_pgradient))
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
        const auto& f = integrand.getLHS().getDerived().getLHS();
        const auto& trans = polytope.getTransformation();
        const auto& fes = trial.getFiniteElementSpace();
        const auto& fe = fes.getFiniteElement(d, idx);
        const size_t dofs = fe.getCount();
        assert(dofs == trial.getDOFs(polytope));
        const QF::QF1P1 qf(polytope.getGeometry());
        auto& res = getMatrix();
        res = Math::Matrix::Zero(dofs, dofs);
        m_rgradient.resize(dofs);
        m_pgradient.resize(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
          const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
          const auto jacInvT = p.getJacobianInverse().transpose();
          for (size_t local = 0; local < dofs; local++)
          {
            fe.getGradient(local)(m_rgradient[local], rc);
            m_pgradient[local] = jacInvT * m_rgradient[local];
          }
          const auto fv = f.getValue(p);
          for (size_t i = 0; i < dofs; i++)
            res(i, i) += qf.getWeight(k) * p.getDistortion() * (fv * m_pgradient[i]).dot(m_pgradient[i]);
          for (size_t i = 0; i < dofs; i++)
            for (size_t j = 0; j < i; j++)
              res(i, j) += qf.getWeight(k) * p.getDistortion() * (fv * m_pgradient[i]).dot(m_pgradient[j]);
        }
        res.template triangularView<Eigen::Upper>() = res.transpose();
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
      std::vector<Math::SpatialVector> m_pgradient;
      std::vector<Math::SpatialVector> m_rgradient;
  };

  template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Context>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Mult<
          FunctionBase<LHSFunctionDerived>,
          ShapeFunctionBase<
            Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>>>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int \nabla u \cdot \nabla v \ dx
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>,
          P1<Scalar, Context>, TrialSpace>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>,
          P1<Scalar, Context>, TestSpace>>>
    : public BilinearFormIntegratorBase
  {
    public:
      using Parent = BilinearFormIntegratorBase;

      using FES = P1<Scalar, Context>;

      using LHS =
        ShapeFunctionBase<
          Grad<ShapeFunction<LHSDerived, FES, TrialSpace>>>;

      using RHS =
        ShapeFunctionBase<
          Grad<ShapeFunction<RHSDerived, FES, TestSpace>>>;

      using Integrand = Dot<LHS, RHS>;

      constexpr
      QuadratureRule(const Integrand& integrand)
        : BilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
          m_integrand(integrand.copy())
      {}

      constexpr
      QuadratureRule(const QuadratureRule& other)
        : Parent(other),
          m_integrand(other.m_integrand->copy()),
          m_rgradient(other.m_rgradient),
          m_pgradient(other.m_pgradient)
      {}

      constexpr
      QuadratureRule(QuadratureRule&& other)
        : Parent(std::move(other)),
          m_integrand(std::move(other.m_integrand)),
          m_rgradient(std::move(other.m_rgradient)),
          m_pgradient(std::move(other.m_pgradient))
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
        m_pgradient.resize(dofs);
        for (size_t k = 0; k < qf.getSize(); k++)
        {
          Geometry::Point p(polytope, trans, std::ref(qf.getPoint(k)));
          const Math::Vector& rc = p.getCoordinates(Geometry::Point::Coordinates::Reference);
          const auto jacInvT = p.getJacobianInverse().transpose();
          for (size_t local = 0; local < dofs; local++)
            m_pgradient[local] = jacInvT * fe.getGradient(local)(rc);
          for (size_t i = 0; i < dofs; i++)
            res(i, i) += qf.getWeight(k) * p.getDistortion() * m_pgradient[i].squaredNorm();
          for (size_t i = 0; i < dofs; i++)
            for (size_t j = 0; j < i; j++)
              res(i, j) += qf.getWeight(k) * p.getDistortion() * m_pgradient[i].dot(m_pgradient[j]);
        }
        res.template triangularView<Eigen::Upper>() = res.transpose();
      }

      virtual Region getRegion() const override = 0;

      virtual QuadratureRule* copy() const noexcept override = 0;

    private:
      std::unique_ptr<Integrand> m_integrand;
      std::vector<Math::SpatialVector> m_rgradient;
      std::vector<Math::SpatialVector> m_pgradient;
  };

  template <class LHSDerived, class RHSDerived, class Context>
  QuadratureRule(const
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>&)
  ->
  QuadratureRule<
    Dot<
      ShapeFunctionBase<
        Grad<ShapeFunction<LHSDerived, P1<Scalar, Context>, TrialSpace>>>,
      ShapeFunctionBase<
        Grad<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   * \int f \cdot v \ dx
   * @f]
   */
  template <class LHSDerived, class RHSDerived, class Context>
  class QuadratureRule<
    ShapeFunctionBase<
      Dot<
        FunctionBase<LHSDerived>,
        ShapeFunctionBase<
          ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>,
            P1<Scalar, Context>, TestSpace>>,
      P1<Scalar, Context>, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using FES = P1<Scalar, Context>;

      using LHS = FunctionBase<LHSDerived>;

      using RHS = ShapeFunctionBase<ShapeFunction<RHSDerived, FES, TestSpace>>;

      using Integrand = ShapeFunctionBase<Dot<LHS, RHS>>;

      using Parent = LinearFormIntegratorBase;

      using LHSRange = typename FormLanguage::Traits<LHS>::RangeType;

      using RHSRange = typename FormLanguage::Traits<RHS>::RangeType;

      static_assert(std::is_same_v<LHSRange, RHSRange>);

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
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& integrand = getIntegrand().getDerived();
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
            return;
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
      ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>&)
  ->
  QuadratureRule<ShapeFunctionBase<
    Dot<
      FunctionBase<LHSDerived>,
      ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Scalar, Context>, TestSpace>>>>>;

  /**
   * @ingroup QuadratureRuleSpecializations
   *
   * @f[
   *  \int v \ dx
   * @f]
   */
  template <class NestedDerived, class Context>
  class QuadratureRule<
    ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>,
        P1<Scalar, Context>, TestSpace>>
    : public LinearFormIntegratorBase
  {
    public:
      using Parent = LinearFormIntegratorBase;

      using FES = P1<Scalar, Context>;

      using Integrand =
        ShapeFunctionBase<
          ShapeFunction<NestedDerived, FES, TestSpace>>;

      using IntegrandRange = typename FormLanguage::Traits<Integrand>::RangeType;

      static_assert(std::is_same_v<IntegrandRange, Scalar>);

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
        const size_t d = polytope.getDimension();
        const Index idx = polytope.getIndex();
        const auto& trans = polytope.getTransformation();
        const auto& integrand = getIntegrand().getDerived();
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
  QuadratureRule(const
    ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>>&)
  ->
  QuadratureRule<
    ShapeFunctionBase<
      ShapeFunction<NestedDerived, P1<Scalar, Context>, TestSpace>>>;
}

#endif

