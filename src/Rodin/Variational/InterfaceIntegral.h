/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTERFACEINTEGRAL_H
#define RODIN_VARIATIONAL_INTERFACEINTEGRAL_H

#include <cassert>
#include <set>
#include <utility>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Utility/MFEM.h"

#include "Dot.h"
#include "LinearForm.h"
#include "ForwardDecls.h"
#include "GridFunction.h"
#include "Function.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "MatrixFunction.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

#include "Utility.h"

namespace Rodin::Variational
{
  /**
   * @defgroup InterfaceIntegralSpecializations InterfaceIntegral Template Specializations
   * @brief Template specializations of the InterfaceIntegral class.
   *
   * @see InterfaceIntegral
   */

  /**
   * @ingroup InterfaceIntegralSpecializations
   * @brief Interface integration of the dot product of a trial and test operators.
   *
   * Given two operators defined over trial and test spaces @f$ U_h
   * @f$ and @f$ V_h @f$,
   * @f[
   *   A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
   * @f]
   * this class represents the integral of their dot product:
   * @f[
   *   \int_{\mathcal{I}_h} A(u) : B(v) \ dx
   * @f]
   * over the interface @f$ \mathcal{I}_h @f$ of the triangulation @f$
   * \mathcal{T}_h @f$.
   */
  // template <>
  // class InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  //   : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent   = Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
  //     using Parent::Parent;
  //     Region getRegion() const override { return Region::Interface; }
  //     InterfaceIntegral* copy() const noexcept override { return new InterfaceIntegral(*this); }
  // };
  // InterfaceIntegral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
  //   -> InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  // InterfaceIntegral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
  //   -> InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

  // template <>
  // class FaceIntegral<ShapeFunctionBase<TestSpace>> : public Integral<ShapeFunctionBase<TestSpace>>
  // {
  //   public:
  //     using Parent   = Integral<ShapeFunctionBase<TestSpace>>;
  //     using Integrand = ShapeFunctionBase<TestSpace>;
  //     using Parent::Parent;
  //     Region getRegion() const override { return Region::Faces; }
  //     FaceIntegral* copy() const noexcept override { return new FaceIntegral(*this); }
  // };
  // FaceIntegral(const FunctionBase&, const ShapeFunctionBase<TestSpace>&)
  //   -> FaceIntegral<ShapeFunctionBase<TestSpace>>;
  // FaceIntegral(const ShapeFunctionBase<TestSpace>&)
  //   -> FaceIntegral<ShapeFunctionBase<TestSpace>>;
}

#endif


