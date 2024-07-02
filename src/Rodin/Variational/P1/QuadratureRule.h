#ifndef RODIN_VARIATIONAL_P1_QUADRATURERULE_H
#define RODIN_VARIATIONAL_P1_QUADRATURERULE_H

#include "Rodin/Variational/QuadratureRule.h"
#include "Rodin/QF/QF1P1.h"
#include "Rodin/QF/GrundmannMoller.h"

#include "P1.h"
#include "P1Element.h"

// namespace Rodin::Variational
// {
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of a P1 ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int v \ dx \: ,
//    * @f]
//    * where @f$ v \in \mathbb{P}_1 @f$.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class NestedDerived, class Range, class Mesh>
//   class QuadratureRule<
//     ShapeFunctionBase<
//       ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>
//     : public LinearFormIntegratorBase<
//         typename FormLanguage::Traits<
//           typename FormLanguage::Traits<
//             ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>>::RangeType>::NumberType>
//   {
//     public:
//       using FESType = P1<Range, Mesh>;
// 
//       using NumberType = typename FormLanguage::Traits<FESType>::NumberType;
// 
//       using IntegrandType =
//         ShapeFunctionBase<ShapeFunction<NestedDerived, FESType, TestSpace>>;
// 
//       using IntegrandRangeType = typename FormLanguage::Traits<FESType>::RangeType;
// 
//       using Parent = LinearFormIntegratorBase<NumberType>;
// 
//       static_assert(std::is_same_v<IntegrandRangeType, NumberType>);
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : Parent(integrand.getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy())
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand))
//       {}
// 
//       inline
//       constexpr
//       const IntegrandType& getIntegrand() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) final override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& trans = polytope.getTransformation();
//         const auto& integrand = getIntegrand().getDerived();
//         const auto& fes = integrand.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         const Geometry::Point p(polytope, trans, std::cref(rc));
//         const Real distortion = p.getDistortion();
//         auto& res = this->getVector();
//         res.resize(dofs);
//         res.setZero();
//         for (size_t local = 0; local < dofs; local++)
//           res.coeffRef(local) += w * distortion * fe.getBasis(local)(rc);
//       }
// 
//       virtual Integrator::Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
//   };
// 
//   /**
//    * @ingroup RodinCTAD
//    */
//   template <class NestedDerived, class Range, class Mesh>
//   QuadratureRule(const ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>&)
//     -> QuadratureRule<ShapeFunctionBase<ShapeFunction<NestedDerived, P1<Range, Mesh>, TestSpace>>>;
// 
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of the Dot product of some coefficient function and a
//    * P1 ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int f \cdot v \ dx \: ,
//    * @f]
//    * where @f$ v \in \mathbb{P}_1 @f$.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int f \cdot v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class LHSDerived, class RHSDerived, class Range, class Mesh>
//   class QuadratureRule<
//     ShapeFunctionBase<
//       Dot<
//         FunctionBase<LHSDerived>,
//         ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>,
//           P1<Range, Mesh>, TestSpace>>
//     : public LinearFormIntegratorBase<
//         typename FormLanguage::Traits<
//           typename FormLanguage::Traits<
//             ShapeFunctionBase<
//               Dot<
//                 FunctionBase<LHSDerived>,
//                 ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>
//           ::RangeType>
//         ::NumberType>
//   {
//     public:
//       using FESType = P1<Range, Mesh>;
// 
//       using LHSType = FunctionBase<LHSDerived>;
// 
//       using RHSType =
//         ShapeFunctionBase<ShapeFunction<RHSDerived, FESType, TestSpace>, FESType, TestSpace>;
// 
//       using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
// 
//       using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
// 
//       using IntegrandType = ShapeFunctionBase<Dot<LHSType, RHSType>>;
// 
//       using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
// 
//       using NumberType = typename FormLanguage::Traits<IntegrandType>::NumberType;
// 
//       using Parent = LinearFormIntegratorBase<NumberType>;
// 
//       static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : Parent(integrand.getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy())
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand))
//       {}
// 
//       inline
//       constexpr
//       const IntegrandType& getIntegrand() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) final override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& trans = polytope.getTransformation();
//         const auto& integrand = getIntegrand().getDerived();
//         const auto& f = integrand.getLHS();
//         const auto& fes = integrand.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         const Geometry::Point p(polytope, trans, std::cref(rc));
//         const Real distortion = p.getDistortion();
//         auto& res = this->getVector();
//         res.resize(dofs);
//         res.setZero();
//         if constexpr (std::is_same_v<NumberType, LHSRangeType>)
//         {
//           for (size_t local = 0; local < dofs; local++)
//             res.coeffRef(local) += w * distortion * f(p) * fe.getBasis(local)(rc);
//         }
//         else if constexpr (std::is_same_v<Math::Vector<NumberType>, LHSRangeType>)
//         {
//           for (size_t local = 0; local < dofs; local++)
//             res.coeffRef(local) += w * distortion * f(p).dot(fe.getBasis(local)(rc));
//         }
//         else if constexpr (std::is_same_v<Math::Matrix<NumberType>, LHSRangeType>)
//         {
//           for (size_t local = 0; local < dofs; local++)
//             res.coeffRef(local) +=
//               w * distortion * (f(p).array() * fe.getBasis(local)(rc).array()).rowwise().sum().colwise().sum().value();
//         }
//         else
//         {
//           assert(false);
//           res.setConstant(NAN);
//           return;
//         }
//       }
// 
//       virtual Integrator::Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
//   };
// 
//   /**
//    * @ingroup RodinCTAD
//    */
//   template <class LHSDerived, class RHSDerived, class Range, class Mesh>
//   QuadratureRule(
//       const ShapeFunctionBase<
//         Dot<FunctionBase<LHSDerived>,
//         ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
//     -> QuadratureRule<ShapeFunctionBase<
//         Dot<
//           FunctionBase<LHSDerived>,
//           ShapeFunctionBase<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;
// 
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of the anisotropic Dot product of two instances of the
//    * P1 ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int (A u) \cdot v \ dx \: ,
//    * @f]
//    * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$
//    * A @f$ is a coefficient function.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int (A u) \cdot v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash u, v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Range, class Mesh>
//   class QuadratureRule<
//     Dot<
//       ShapeFunctionBase<
//         Mult<
//           FunctionBase<CoefficientDerived>,
//           ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>,
//             P1<Range, Mesh>, TrialSpace>>, P1<Range, Mesh>, TrialSpace>,
//       ShapeFunctionBase<
//         ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>, P1<Range, Mesh>, TestSpace>>>
//     : public LocalBilinearFormIntegratorBase<
//         typename FormLanguage::Traits<
//           Dot<
//             ShapeFunctionBase<
//               Mult<
//                 FunctionBase<CoefficientDerived>,
//                 ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>,
//             ShapeFunctionBase<
//               ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>::NumberType>
//   {
//     public:
//       using FESType = P1<Range, Mesh>;
// 
//       using CoefficientType = FunctionBase<CoefficientDerived>;
// 
//       using MultiplicandType =
//         Mult<
//           FunctionBase<CoefficientDerived>,
//           ShapeFunctionBase<ShapeFunction<LHSDerived, FESType, TrialSpace>>>;
// 
//       using LHSType =
//         ShapeFunctionBase<
//           Mult<
//             FunctionBase<CoefficientDerived>,
//             ShapeFunctionBase<ShapeFunction<LHSDerived, FESType, TrialSpace>>>>;
// 
//       using RHSType =
//         ShapeFunctionBase<
//           ShapeFunction<RHSDerived, FESType, TestSpace>>;
// 
//       using IntegrandType = Dot<LHSType, RHSType>;
// 
//       using CoefficientRangeType = typename FormLanguage::Traits<CoefficientType>::RangeType;
// 
//       using MultiplicandRangeType = typename FormLanguage::Traits<MultiplicandType>::RangeType;
// 
//       using LHSRangeType = typename FormLanguage::Traits<LHSType>::RangeType;
// 
//       using RHSRangeType = typename FormLanguage::Traits<RHSType>::RangeType;
// 
//       using NumberType = typename FormLanguage::Traits<IntegrandType>::NumberType;
// 
//       using Parent = LocalBilinearFormIntegratorBase<NumberType>;
// 
// 
//       static_assert(std::is_same_v<LHSRangeType, RHSRangeType>);
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : LocalBilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy())
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand))
//       {}
// 
//       inline
//       constexpr
//       const IntegrandType& getIntegrandType() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& integrand = getIntegrandType();
//         const auto& coeff = integrand.getLHS().getDerived().getLHS();
//         const auto& multiplicand = integrand.getLHS();
//         const auto& trans = polytope.getTransformation();
//         const auto& fes = multiplicand.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         auto& res = getMatrix();
//         res.resize(dofs, dofs);
//         res.setZero();
//         if constexpr (std::is_same_v<CoefficientRangeType, NumberType>)
//         {
//           static_assert(std::is_same_v<MultiplicandRangeType, RHSRangeType>);
//           if constexpr (std::is_same_v<MultiplicandRangeType, NumberType>)
//           {
//             const Geometry::Point p(polytope, trans, std::cref(rc));
//             const Real distortion = p.getDistortion();
//             const Real c = coeff.getValue(p);
//             for (size_t i = 0; i < dofs; i++)
//             {
//               const NumberType basis = fe.getBasis(i)(rc);
//               res(i, i) += w * distortion * c * basis * basis;
//             }
// 
//             for (size_t i = 0; i < dofs; i++)
//               for (size_t j = 0; j < i; j++)
//                 res(i, j) += w * distortion * c * fe.getBasis(i)(rc) * fe.getBasis(j)(rc);
//             res.template triangularView<Eigen::Upper>() = res.transpose();
//           }
//           else
//           {
//             assert(false); // Not handled yet
//             res.setConstant(NAN);
//             return;
//           }
//         }
//         else
//         {
//           assert(false); // Not handled yet
//           res.setConstant(NAN);
//           return;
//         }
//       }
// 
//       virtual Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
//       std::vector<Math::Vector<NumberType>> m_vvalues;
//       std::vector<Math::Matrix<NumberType>> m_mvalues;
//   };
// 
//   template <class CoefficientDerived, class LHSDerived, class RHSDerived, class Number, class Mesh>
//   QuadratureRule(const
//     Dot<
//       ShapeFunctionBase<
//         Mult<
//           FunctionBase<CoefficientDerived>,
//           ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Number, Mesh>, TrialSpace>>>>,
//       ShapeFunctionBase<
//         ShapeFunction<RHSDerived, P1<Number, Mesh>, TestSpace>>>&)
//   ->
//   QuadratureRule<
//     Dot<
//       ShapeFunctionBase<
//         Mult<
//           FunctionBase<CoefficientDerived>,
//           ShapeFunctionBase<ShapeFunction<LHSDerived, P1<Number, Mesh>, TrialSpace>>>>,
//       ShapeFunctionBase<
// 
//         ShapeFunction<RHSDerived, P1<Number, Mesh>, TestSpace>>>>;
// 
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of the isotropic Dot product of two instances of the P1
//    * Grad of ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int \nabla u \cdot \nabla v \ dx \: ,
//    * @f]
//    * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$ A
//    * @f$ is a coefficient function.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int \nabla u \cdot \nabla v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash u, v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class LHSDerived, class RHSDerived, class Range, class Mesh>
//   class QuadratureRule<
//     Dot<
//       ShapeFunctionBase<
//         Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>, P1<Range, Mesh>, TrialSpace>,
//       ShapeFunctionBase<
//         Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>, P1<Range, Mesh>, TestSpace>>>
//     : public LocalBilinearFormIntegratorBase
//   {
//     public:
//       using FESType = P1<Range, Mesh>;
// 
//       using LHSType =
//         ShapeFunctionBase<
//           Grad<ShapeFunction<LHSDerived, FESType, TrialSpace>>>;
// 
//       using RHSType =
//         ShapeFunctionBase<
//           Grad<ShapeFunction<RHSDerived, FESType, TestSpace>>>;
// 
//       using IntegrandType = Dot<LHSType, RHSType>;
// 
//       using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
// 
//       using NumberType = typename FormLanguage::Traits<IntegrandRangeType>::NumberType;
// 
//       using Parent = LocalBilinearFormIntegratorBase;
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : LocalBilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy()),
//           m_rgradient(other.m_rgradient),
//           m_pgradient(other.m_pgradient)
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand)),
//           m_rgradient(std::move(other.m_rgradient)),
//           m_pgradient(std::move(other.m_pgradient))
//       {}
// 
//       inline
//       constexpr
//       const IntegrandType& getIntegrandType() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& integrand = getIntegrandType();
//         const auto& trial = integrand.getLHS();
//         const auto& trans = polytope.getTransformation();
//         const auto& fes = trial.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         assert(dofs == trial.getDOFs(polytope));
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         const Geometry::Point p(polytope, trans, std::cref(rc));
//         const Real distortion = p.getDistortion();
//         const auto jacInvT = p.getJacobianInverse().transpose();
//         auto& res = getMatrix();
//         res.resize(dofs, dofs);
//         res.setZero();
//         m_pgradient.resize(dofs);
//         for (size_t local = 0; local < dofs; local++)
//           m_pgradient[local] = jacInvT * fe.getGradient(local)(rc);
//         for (size_t i = 0; i < dofs; i++)
//           res(i, i) += w * distortion * m_pgradient[i].squaredNorm();
//         for (size_t i = 0; i < dofs; i++)
//           for (size_t j = 0; j < i; j++)
//             res(i, j) += w * distortion * m_pgradient[i].dot(m_pgradient[j]);
//         res.template triangularView<Eigen::Upper>() = res.transpose();
//       }
// 
//       virtual Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
//       std::vector<Math::SpatialVector<NumberType>> m_rgradient;
//       std::vector<Math::SpatialVector<NumberType>> m_pgradient;
//   };
// 
//   template <class LHSDerived, class RHSDerived, class Range, class Mesh>
//   QuadratureRule(
//       const Dot<
//         ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
//         ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
//     -> QuadratureRule<Dot<
//           ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
//           ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;
// 
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of the anisotropic Dot product of two instances of the
//    * P1 Grad of ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int (A \nabla u) \cdot \nabla v \ dx \: ,
//    * @f]
//    * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in \mathbb{P}_1 @f$, and @f$ A
//    * @f$ is a coefficient function.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int (A \nabla u) \cdot \nabla v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash u, v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Range, class Mesh>
//   class QuadratureRule<
//   Dot<
//     ShapeFunctionBase<
//       Mult<
//         FunctionBase<LHSFunctionDerived>,
//         ShapeFunctionBase<
//           Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>, P1<Range, Mesh>, TrialSpace>>,
//       P1<Range, Mesh>, TrialSpace>,
//     ShapeFunctionBase<
//       Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>, P1<Range, Mesh>, TestSpace>>>
//     : public LocalBilinearFormIntegratorBase
//   {
//     public:
//       using FESType = P1<Range, Mesh>;
// 
//       using LHSType =
//         ShapeFunctionBase<
//           Mult<
//             FunctionBase<LHSFunctionDerived>,
//             ShapeFunctionBase<
//               Grad<ShapeFunction<LHSDerived, FESType, TrialSpace>>>>>;
// 
//       using RHSType =
//         ShapeFunctionBase<
//           Grad<ShapeFunction<RHSDerived, FESType, TestSpace>>>;
// 
//       using IntegrandType = Dot<LHSType, RHSType>;
// 
//       using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
// 
//       using NumberType = typename FormLanguage::Traits<IntegrandRangeType>::NumberType;
// 
//       using Parent = LocalBilinearFormIntegratorBase;
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : LocalBilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy()),
//           m_rgradient(other.m_rgradient),
//           m_pgradient(other.m_pgradient)
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand)),
//           m_rgradient(std::move(other.m_rgradient)),
//           m_pgradient(std::move(other.m_pgradient))
//       {}
// 
//       inline
//       constexpr
//       const IntegrandType& getIntegrandType() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& integrand = getIntegrandType();
//         const auto& trial = integrand.getLHS();
//         const auto& f = integrand.getLHS().getDerived().getLHS();
//         const auto& trans = polytope.getTransformation();
//         const auto& fes = trial.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         assert(dofs == trial.getDOFs(polytope));
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         const Geometry::Point p(polytope, trans, std::cref(rc));
//         const Real distortion = p.getDistortion();
//         const auto jacInvT = p.getJacobianInverse().transpose();
//         auto& res = getMatrix();
//         res.resize(dofs, dofs);
//         res.setZero();
//         m_rgradient.resize(dofs);
//         m_pgradient.resize(dofs);
//         for (size_t local = 0; local < dofs; local++)
//         {
//           fe.getGradient(local)(m_rgradient[local], rc);
//           m_pgradient[local] = jacInvT * m_rgradient[local];
//         }
//         const auto fv = f.getValue(p);
//         for (size_t i = 0; i < dofs; i++)
//           res(i, i) += w * distortion * (fv * m_pgradient[i]).dot(m_pgradient[i]);
//         for (size_t i = 0; i < dofs; i++)
//           for (size_t j = 0; j < i; j++)
//             res(i, j) += w * distortion * (fv * m_pgradient[i]).dot(m_pgradient[j]);
//         res.template triangularView<Eigen::Upper>() = res.transpose();
//       }
// 
//       virtual Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
//       std::vector<Math::SpatialVector<NumberType>> m_rgradient;
//       std::vector<Math::SpatialVector<NumberType>> m_pgradient;
//   };
// 
//   /**
//    * @ingroup RodinCTAD
//    */
//   template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Range, class Mesh>
//   QuadratureRule(
//       const Dot<
//         ShapeFunctionBase<
//           Mult<
//             FunctionBase<LHSFunctionDerived>,
//             ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>>,
//         ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
//     -> QuadratureRule<
//           Dot<ShapeFunctionBase<
//             Mult<
//               FunctionBase<LHSFunctionDerived>,
//               ShapeFunctionBase<Grad<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>>>,
//           ShapeFunctionBase<Grad<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;
// 
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of the isotropic Frobenius inner product two instances
//    * of the P1 Jacobian of ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int \mathbf{J} \: u : \mathbf{J} \: v \ dx \: ,
//    * @f]
//    * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in
//    * \mathbb{P}_1 @f$, and @f$ A @f$ is coefficient function.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int \mathbf{J} \: u : \mathbf{J} \: v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash u, v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class LHSDerived, class RHSDerived, class Range, class Mesh>
//   class QuadratureRule<
//     Dot<
//       ShapeFunctionBase<
//         Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>,
//           P1<Range, Mesh>, TrialSpace>,
//       ShapeFunctionBase<
//         Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>,
//           P1<Range, Mesh>, TestSpace>>>
//     : public LocalBilinearFormIntegratorBase
//   {
//     public:
//       using FESType = P1<Range, Mesh>;
// 
//       using LHSType =
//         ShapeFunctionBase<
//           Jacobian<ShapeFunction<LHSDerived, FESType, TrialSpace>>>;
// 
//       using RHSType =
//         ShapeFunctionBase<
//           Jacobian<ShapeFunction<RHSDerived, FESType, TestSpace>>>;
// 
//       using IntegrandType = Dot<LHSType, RHSType>;
// 
//       using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
// 
//       using NumberType = typename FormLanguage::Traits<IntegrandRangeType>::NumberType;
// 
//       using Parent = LocalBilinearFormIntegratorBase;
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : LocalBilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy()),
//           m_rjac(other.m_rjac),
//           m_pjac(other.m_pjac)
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand)),
//           m_rjac(std::move(other.m_rjac)),
//           m_pjac(std::move(other.m_pjac))
//       {}
// 
//       inline
//       constexpr
//       const IntegrandType& getIntegrandType() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& integrand = getIntegrandType();
//         const auto& trial = integrand.getLHS();
//         const auto& trans = polytope.getTransformation();
//         const auto& fes = trial.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         assert(dofs == trial.getDOFs(polytope));
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         const Geometry::Point p(polytope, trans, std::cref(rc));
//         const Real distortion = p.getDistortion();
//         auto& res = getMatrix();
//         res.resize(dofs, dofs);
//         res.setZero();
//         m_pjac.resize(dofs);
//         for (size_t local = 0; local < dofs; local++)
//           m_pjac[local] = fe.getJacobian(local)(rc) * p.getJacobianInverse();
//         for (size_t i = 0; i < dofs; i++)
//           res(i, i) += w * distortion * m_pjac[i].squaredNorm();
//         for (size_t i = 0; i < dofs; i++)
//           for (size_t j = 0; j < i; j++)
//             res(i, j) += w * distortion * m_pjac[i].dot(m_pjac[j]);
//         res.template triangularView<Eigen::Upper>() = res.transpose();
//       }
// 
//       virtual Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
//       std::vector<Math::SpatialVector<NumberType>> m_rjac;
//       std::vector<Math::SpatialVector<NumberType>> m_pjac;
//   };
// 
//   /**
//    * @ingroup RodinCTAD
//    */
//   template <class LHSDerived, class RHSDerived, class Range, class Mesh>
//   QuadratureRule(const
//     Dot<
//       ShapeFunctionBase<
//         Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
//       ShapeFunctionBase<
//         Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>&)
//     ->
//       QuadratureRule<
//         Dot<
//           ShapeFunctionBase<
//             Jacobian<ShapeFunction<LHSDerived, P1<Range, Mesh>, TrialSpace>>>,
//           ShapeFunctionBase<
//             Jacobian<ShapeFunction<RHSDerived, P1<Range, Mesh>, TestSpace>>>>>;
// 
//   /**
//    * @ingroup QuadratureRuleSpecializations
//    * @brief Integration of the anisotropic Frobenius inner product two
//    * instances of the P1 Jacobian of ShapeFunction.
//    *
//    * This class represents the CTAD for the expression:
//    * @f[
//    * \int (A \: \mathbf{J} \: u) : \mathbf{J} \: v \ dx \: ,
//    * @f]
//    * where @f$ u \in \mathbb{P}_1 @f$ and @f$ v \in
//    * \mathbb{P}_1 @f$, and @f$ A @f$ is a coefficient function.
//    *
//    * Judgement
//    * ---------
//    *
//    * The following judgement specifies that the expression is a well formed type
//    * of QuadratureRule.
//    * @f[
//    * \dfrac
//    * {\vdash \int (A \: \mathbf{J} \: u) : \mathbf{J} \: v \ dx :
//    * \texttt{QuadratureRule}}
//    * {\vdash u, v : \mathbb{P}_1}
//    * @f]
//    */
//   template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Mesh>
//   class QuadratureRule<
//   Dot<
//     ShapeFunctionBase<
//       Mult<
//         FunctionBase<LHSFunctionDerived>,
//         ShapeFunctionBase<
//           Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>,
//         P1<Math::Vector<Real>, Mesh>, TrialSpace>>,
//       P1<Math::Vector<Real>, Mesh>, TrialSpace>,
//     ShapeFunctionBase<
//       Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>,
//     P1<Math::Vector<Real>, Mesh>, TestSpace>>>
//     : public LocalBilinearFormIntegratorBase
//   {
//     public:
//       /// Type of finite element space
//       using FESType = P1<Math::Vector<Real>, Mesh>;
// 
//       /// Type of left hand side
//       using LHSType =
//         ShapeFunctionBase<
//           Mult<
//             FunctionBase<LHSFunctionDerived>,
//             ShapeFunctionBase<
//               Jacobian<ShapeFunction<LHSDerived, FESType, TrialSpace>>>>>;
// 
//       /// Type of right hand side
//       using RHSType =
//         ShapeFunctionBase<
//           Jacobian<ShapeFunction<RHSDerived, FESType, TestSpace>>>;
// 
//       /// Type of integrand
//       using IntegrandType = Dot<LHSType, RHSType>;
// 
//       using IntegrandRangeType = typename FormLanguage::Traits<IntegrandType>::RangeType;
// 
//       using NumberType = typename FormLanguage::Traits<IntegrandRangeType>::NumberType;
// 
//       using Parent = LocalBilinearFormIntegratorBase;
// 
//       constexpr
//       QuadratureRule(const IntegrandType& integrand)
//         : LocalBilinearFormIntegratorBase(integrand.getLHS().getLeaf(), integrand.getRHS().getLeaf()),
//           m_integrand(integrand.copy())
//       {}
// 
//       constexpr
//       QuadratureRule(const QuadratureRule& other)
//         : Parent(other),
//           m_integrand(other.m_integrand->copy()),
//           m_rjac(other.m_rjac),
//           m_pjac(other.m_pjac)
//       {}
// 
//       constexpr
//       QuadratureRule(QuadratureRule&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand)),
//           m_rjac(std::move(other.m_rjac)),
//           m_pjac(std::move(other.m_pjac))
//       {}
// 
//       /**
//        * @brief Gets the integrand.
//        */
//       inline
//       constexpr
//       const IntegrandType& getIntegrandType() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       void assemble(const Geometry::Polytope& polytope) override
//       {
//         const size_t d = polytope.getDimension();
//         const Index idx = polytope.getIndex();
//         const auto& integrand = getIntegrandType();
//         const auto& trial = integrand.getLHS();
//         const auto& f = integrand.getLHS().getDerived().getLHS();
//         const auto& trans = polytope.getTransformation();
//         const auto& fes = trial.getFiniteElementSpace();
//         const auto& fe = fes.getFiniteElement(d, idx);
//         const size_t dofs = fe.getCount();
//         assert(dofs == trial.getDOFs(polytope));
//         const QF::QF1P1 qf(polytope.getGeometry());
//         assert(qf.getSize() == 1);
//         const Real w = qf.getWeight(0);
//         const auto& rc = qf.getPoint(0);
//         const Geometry::Point p(polytope, trans, std::cref(rc));
//         const Real distortion = p.getDistortion();
//         auto& res = getMatrix();
//         res = Math::Matrix<Real>::Zero(dofs, dofs);
//         m_rjac.resize(dofs);
//         m_pjac.resize(dofs);
//         for (size_t local = 0; local < dofs; local++)
//         {
//           fe.getJacobian(local)(m_rjac[local], rc);
//           m_pjac[local] = m_rjac[local] * p.getJacobianInverse();
//         }
//         const auto fv = f.getValue(p);
//         for (size_t i = 0; i < dofs; i++)
//         {
//           const auto lhs = fv * m_pjac[i];
//           const auto rhs = m_pjac[i];
//           res(i, i) += w * distortion * (
//               lhs.array() * rhs.array()).rowwise().sum().colwise().sum().value();
//         }
//         for (size_t i = 0; i < dofs; i++)
//         {
//           const auto lhs = fv * m_pjac[i];
//           for (size_t j = 0; j < i; j++)
//           {
//             const auto rhs = m_pjac[j];
//             res(i, j) += w * distortion * (
//                 lhs.array() * rhs.array()).rowwise().sum().colwise().sum().value();
//           }
//         }
//         res.template triangularView<Eigen::Upper>() = res.transpose();
//       }
// 
//       virtual Region getRegion() const override = 0;
// 
//       virtual QuadratureRule* copy() const noexcept override = 0;
// 
//     private:
//       std::unique_ptr<IntegrandType> m_integrand;
// 
//       std::vector<Math::SpatialMatrix<Real>> m_rjac;
//       std::vector<Math::SpatialMatrix<Real>> m_pjac;
//   };
// 
//   /**
//    * @ingroup RodinCTAD
//    */
//   template <class LHSFunctionDerived, class LHSDerived, class RHSDerived, class Mesh>
//   QuadratureRule(const
//     Dot<
//       ShapeFunctionBase<
//         Mult<
//           FunctionBase<LHSFunctionDerived>,
//           ShapeFunctionBase<Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>>>>,
//       ShapeFunctionBase<
//         Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>>>&)
//   ->
//   QuadratureRule<
//     Dot<
//       ShapeFunctionBase<
//         Mult<
//           FunctionBase<LHSFunctionDerived>,
//           ShapeFunctionBase<
//             Jacobian<ShapeFunction<LHSDerived, P1<Math::Vector<Real>, Mesh>, TrialSpace>>>>>,
//       ShapeFunctionBase<
//         Jacobian<ShapeFunction<RHSDerived, P1<Math::Vector<Real>, Mesh>, TestSpace>>>>>;
// }

#endif

