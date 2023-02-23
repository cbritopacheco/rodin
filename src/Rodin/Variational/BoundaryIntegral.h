/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H
#define RODIN_VARIATIONAL_BOUNDARYINTEGRAL_H

#include <cassert>
#include <set>
#include <utility>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "Dot.h"
#include "LinearForm.h"
#include "ForwardDecls.h"
#include "GridFunction.h"
#include "Function.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "FiniteElement.h"
#include "MatrixFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

// namespace Rodin::Variational
// {
//   /**
//    * @defgroup BoundaryIntegralSpecializations BoundaryIntegral Template Specializations
//    * @brief Template specializations of the BoundaryIntegral class.
//    *
//    * @see BoundaryIntegral
//    */
// 
//   /**
//    * @ingroup BoundaryIntegralSpecializations
//    * @brief Boundary integration of the dot product of a trial and test operators.
//    */
//   template <class LHSDerived, class RHSDerived>
//   class BoundaryIntegral<Dot<ShapeFunctionBase<LHSDerived, TrialSpace>, ShapeFunctionBase<RHSDerived, TestSpace>>> final
//     : public BilinearFormIntegratorBase
//   {
//     public:
//       using LHS = ShapeFunctionBase<LHSDerived, TrialSpace>;
//       using RHS = ShapeFunctionBase<LHSDerived, TestSpace>;
//       using Integrand = Dot<LHS, RHS>;
//       using Parent = BilinearFormIntegratorBase;
// 
//       using IntegrationOrder =
//         std::function<size_t(
//             const Geometry::Simplex&, const Geometry::SimplexTransformation&,
//             const FiniteElement&, const FiniteElement&)>;
// 
//       inline
//       const Integrand& getIntegrand() const
//       {
//         return m_prod;
//       }
// 
//       inline
//       Region getRegion() const override
//       {
//         return Region::Boundary;
//       }
// 
//       inline
//       Math::Matrix getMatrix(const Geometry::Simplex& element) const override
//       {
//         assert(false);
//       }
// 
//       inline
//       BoundaryIntegral* copy() const noexcept override
//       {
//         return new BoundaryIntegral(*this);
//       }
//     private:
//       Integrand m_prod;
//       IntegrationOrder m_intOrder;
//   };
// 
//   template <class LHSDerived, class RHSDerived>
//   BoundaryIntegral(
//       const ShapeFunctionBase<LHSDerived, TrialSpace>&,
//       const ShapeFunctionBase<RHSDerived, TestSpace>&)
//     -> BoundaryIntegral<Dot<
//         ShapeFunctionBase<LHSDerived, TrialSpace>, ShapeFunctionBase<RHSDerived, TestSpace>>>;
// 
//   template <class LHSDerived, class RHSDerived>
//   BoundaryIntegral(
//       const Dot<ShapeFunctionBase<LHSDerived, TrialSpace>,
//       ShapeFunctionBase<RHSDerived, TestSpace>>&)
//     -> BoundaryIntegral<Dot<
//         ShapeFunctionBase<LHSDerived, TrialSpace>, ShapeFunctionBase<RHSDerived, TestSpace>>>;
// 
//   template <class NestedDerived>
//   class BoundaryIntegral<ShapeFunctionBase<NestedDerived, TestSpace>> final
//     : public LinearFormIntegratorBase
//   {
//     public:
//       using Integrand = ShapeFunctionBase<NestedDerived, TestSpace>;
//       using Parent = LinearFormIntegratorBase;
// 
//       using IntegrationOrder =
//         std::function<size_t(
//             const Geometry::Simplex&, const Geometry::SimplexTransformation&,
//             const FiniteElement&)>;
// 
//       template <class LHSDerived, class RHSDerived>
//       constexpr
//       BoundaryIntegral(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, TestSpace>& rhs)
//         : BoundaryIntegral(Dot(lhs, rhs))
//       {}
// 
//       /**
//        * @brief Integral of a scalar valued test operator
//        *
//        * Given
//        * @f[
//        *   A : V_h \rightarrow \mathbb{R}
//        * @f]
//        * constructs an instance representing the following integral
//        * @f[
//        *   \int_\Omega A(v) \ dx \ .
//        * @f]
//        */
//       constexpr
//       BoundaryIntegral(const Integrand& integrand)
//         : Parent(integrand.getLeaf()),
//           m_integrand(integrand),
//           m_intOrder(
//               [](const Geometry::Simplex&, const Geometry::SimplexTransformation& trans,
//                  const FiniteElement& fe) -> size_t
//               {
//                 return fe.getOrder() + trans.getHandle().OrderW();
//               })
//       {}
// 
//       constexpr
//       BoundaryIntegral(const BoundaryIntegral& other)
//         : Parent(other),
//           m_integrand(other.m_integrand),
//           m_intOrder(other.m_intOrder)
//       {}
// 
//       constexpr
//       BoundaryIntegral(BoundaryIntegral&& other)
//         : Parent(std::move(other)),
//           m_integrand(std::move(other.m_integrand))
//       {}
// 
//       inline
//       BoundaryIntegral& setIntegrationOrder(IntegrationOrder order)
//       {
//         m_intOrder = order;
//         return *this;
//       }
// 
//       inline
//       size_t getIntegrationOrder(
//           const Geometry::Simplex& simplex, const Geometry::SimplexTransformation& trans,
//           const FiniteElement& fe) const
//       {
//         return m_intOrder(simplex, trans, fe);
//       }
// 
//       inline
//       constexpr
//       const Integrand& getIntegrand() const
//       {
//         assert(m_integrand);
//         return *m_integrand;
//       }
// 
//       inline
//       Region getRegion() const final override
//       {
//         return Region::Domain;
//       }
// 
//       inline
//       Math::Vector getVector(const Geometry::Simplex& simplex) const final override
//       {
// 
//         // const auto& test = *m_integrand;
//         // assert(test.getRangeType() == RangeType::Scalar);
//         // auto& trans = simplex.getTransformation();
//         // assert(false);
//         // const size_t order = 0;
//         // // const size_t order = getIntegrationOrder(test.getFiniteElementSpace(), simplex);
//         // ShapeComputator compute;
// 
//         // Math::Vector res = Math::Vector::Zero(test.getDOFs(simplex));
//         // for (const auto& p : simplex.getIntegrationRule(order))
//         // {
//         //   const TensorBasis basis =
//         //     trans.Weight() * trans.GetIntPoint().weight * test.getOperator(compute, p);
//         //   res += Eigen::Map<const Math::Vector>(basis.data(), basis.size());
//         // }
//         // return res;
//       }
// 
//       virtual BoundaryIntegral* copy() const noexcept override
//       {
//         return new BoundaryIntegral(*this);
//       }
// 
//     private:
//       Integrand m_integrand;
//       IntegrationOrder m_intOrder;
//   };
// 
//   template <class NestedDerived>
//   BoundaryIntegral(const ShapeFunctionBase<NestedDerived, TestSpace>&)
//     -> BoundaryIntegral<ShapeFunctionBase<NestedDerived, TestSpace>>;
// 
//   template <class LHSDerived, class RHSDerived>
//   BoundaryIntegral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, TestSpace>&)
//     -> BoundaryIntegral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, TestSpace>>, TestSpace>>;
// 
//   // template <class FES>
//   // class BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
//   //   : public Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
//   // {
//   //   public:
//   //     using Parent =
//   //       Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
//   //     using Integrand =
//   //       Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;
//   //     using Parent::Parent;
//   //     Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
//   //     BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
//   // };
//   // template <class FES>
//   // BoundaryIntegral(
//   //     const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&,
//   //     const ShapeFunction<FES, TestSpace>&)
//   //   -> BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
//   // template <class FES>
//   // BoundaryIntegral(
//   //     const Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>&)
//   //   -> BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
// 
//   // template <class FES>
//   // class BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>
//   //   : public Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>
//   // {
//   //   public:
//   //     using Parent =
//   //       Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
//   //     using Integrand =
//   //       Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;
//   //     using Parent::Parent;
//   //     Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
//   //     BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
//   // };
//   // template <class FES>
//   // BoundaryIntegral(
//   //     const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&,
//   //     const Grad<ShapeFunction<FES, TestSpace>>&)
//   //   -> BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
//   // template <class FES>
//   // BoundaryIntegral(
//   //     const Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>&)
//   //   -> BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
// 
//   // template <class FES>
//   // class BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
//   //   : public Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
//   // {
//   //   public:
//   //     using Parent    = Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
//   //     using Integrand  = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;
//   //     using Parent::Parent;
//   //     Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
//   //     BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
//   // };
//   // template <class FES>
//   // BoundaryIntegral(const FunctionBase&, const ShapeFunction<FES, TestSpace>&)
//   //   -> BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
//   // template <class FES>
//   // BoundaryIntegral(const Dot<FunctionBase, ShapeFunction<FES, TestSpace>>&)
//   //   -> BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
// }
// 
#endif

