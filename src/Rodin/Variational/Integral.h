/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <cassert>
#include <set>
#include <utility>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "Dot.h"
#include "Function.h"
#include "LinearForm.h"
#include "ForwardDecls.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "FiniteElement.h"
#include "MatrixFunction.h"
#include "GaussianQuadrature.h"
#include "LinearFormIntegrator.h"
#include "BilinearFormIntegrator.h"

namespace Rodin::Variational
{

  /**
   * @defgroup IntegralSpecializations Integral Template Specializations
   * @brief Template specializations of the Integral class.
   *
   * @see Integral
   */

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of the dot product of a trial and test operators.
   *
   * Given two operators defined over trial and test spaces @f$ U_h
   * @f$ and @f$ V_h @f$,
   * @f[
   *   A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
   * @f]
   * this class represents the integral of their dot product:
   * @f[
   *   \int_{\mathcal{T}_h} A(u) : B(v) \ dx
   * @f]
   */
  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  class Integral<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>> final
    : public GaussianQuadrature<Dot<
          ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>,
          ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>
  {
    public:
      using LHS = ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>;
      using RHS = ShapeFunctionBase<RHSDerived, TestFES, TestSpace>;
      using Integrand = Dot<LHS, RHS>;
      using Parent = GaussianQuadrature<Integrand>;

      /**
       * @brief Integral of the dot product of trial and test operators
       *
       * Constructs an instance representing the following integral:
       * @f[
       *   \int_\Omega A(u) : B(v) \ dx
       * @f]
       *
       * @param[in] lhs Trial operator @f$ A(u) @f$
       * @param[in] rhs Test operator @f$ B(v) @f$
       */
      constexpr
      Integral(const LHS& lhs, const RHS& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      /**
       * @brief Integral of the dot product of trial and test operators
       *
       * Constructs the following object representing the following
       * integral:
       * @f[
       *   \int_\Omega A(u) : B(v) \ dx
       * @f]
       *
       * @param[in] prod Dot product instance
       */
      constexpr
      Integral(const Integrand& prod)
        : Parent(prod)
      {}

      constexpr
      Integral(const Integral& other)
        : Parent(other)
      {}

      constexpr
      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Domain;
      }

      inline Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      std::unique_ptr<Integrand> m_prod;
  };

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>&)
    -> Integral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  template <class LHSDerived, class TrialFES, class RHSDerived, class TestFES>
  Integral(const ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>&, const ShapeFunctionBase<RHSDerived, TestFES, TestSpace>&)
    -> Integral<Dot<ShapeFunctionBase<LHSDerived, TrialFES, TrialSpace>, ShapeFunctionBase<RHSDerived, TestFES, TestSpace>>>;

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a test operator.
   *
   * Given an operator defined over a test space @f$ V_h @f$
   * @f[
   *   A : V_h \rightarrow \mathbb{R},
   * @f]
   * this class will represent its integral
   * @f[
   *   \int_\Omega A(v) \ dx \ .
   * @f]
   */
  template <class NestedDerived, class FES>
  class Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>> final
    : public GaussianQuadrature<ShapeFunctionBase<NestedDerived, FES, TestSpace>>
  {
    public:
      using Integrand = ShapeFunctionBase<NestedDerived, FES, TestSpace>;
      using Parent = GaussianQuadrature<Integrand>;

      template <class LHSDerived, class RHSDerived>
      constexpr
      Integral(const FunctionBase<LHSDerived>& lhs, const ShapeFunctionBase<RHSDerived, FES, TestSpace>& rhs)
        : Integral(Dot(lhs, rhs))
      {}

      constexpr
      Integral(const Integrand& integrand)
        : Parent(integrand)
      {}

      constexpr
      Integral(const Integral& other)
        : Parent(other)
      {}

      constexpr
      Integral(Integral&& other)
        : Parent(std::move(other))
      {}

      inline Integrator::Region getRegion() const override
      {
        return Integrator::Region::Domain;
      }

      inline Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      std::unique_ptr<Integrand> m_integrand;
  };

  template <class NestedDerived, class FES>
  Integral(const ShapeFunctionBase<NestedDerived, FES, TestSpace>&)
    -> Integral<ShapeFunctionBase<NestedDerived, FES, TestSpace>>;

  template <class LHSDerived, class RHSDerived, class FES>
  Integral(const FunctionBase<LHSDerived>&, const ShapeFunctionBase<RHSDerived, FES, TestSpace>&)
    -> Integral<ShapeFunctionBase<Dot<FunctionBase<LHSDerived>, ShapeFunctionBase<RHSDerived, FES, TestSpace>>, FES, TestSpace>>;

  /**
   * @ingroup IntegralSpecializations
   * @brief Integration of a GridFunction object.
   */
  template <class FES>
  class Integral<GridFunction<FES>> final : public FormLanguage::Base
  {
    public:
      using Integrand = GridFunction<FES>;
      using Parent = FormLanguage::Base;

      /**
       * @brief Constructs the integral object
       */
      Integral(Integrand& u)
        : m_u(u),
          m_v(u.getFiniteElementSpace()),
          m_one(u.getFiniteElementSpace()),
          m_lf(m_v),
          m_assembled(false)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
        m_one = ScalarFunction(Scalar(1));
        m_lf.from(Integral(u * m_v));
      }

      Integral(const Integral& other)
        : Parent(other),
          m_u(other.m_u),
          m_v(other.m_u.get().getFiniteElementSpace()),
          m_one(other.m_u.get().getFiniteElementSpace()),
          m_lf(m_v),
          m_assembled(false)
      {}

      Integral(Integral&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_v(std::move(other.m_v)),
          m_one(std::move(other.m_one)),
          m_lf(std::move(other.m_lf)),
          m_assembled(std::move(other.m_assembled))
      {}

      /**
       * @brief Integrates the expression and returns the value
       * @returns Value of integral
       */
      inline
      Scalar compute()
      {
        m_lf.assemble();
        m_assembled = true;
        return m_lf(m_one);
      }

      inline
      Integral* copy() const noexcept override
      {
        return new Integral(*this);
      }

    private:
      std::reference_wrapper<GridFunction<FES>>         m_u;
      TestFunction<FES>                                 m_v;
      GridFunction<FES>                                 m_one;

      LinearForm<FES, Context::Serial, mfem::Vector>    m_lf;
      bool m_assembled;
  };

  template <class FES>
  Integral(GridFunction<FES>&) -> Integral<GridFunction<FES>>;

  // /* ||-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  //  * ---------------------------------------------------------------------->>
  //  */

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * @f[
  //  * \int_\Omega \nabla u \cdot \nabla v \ dx
  //  * @f]
  //  */
  // template <class FES>
  // class Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>
  //   : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>;

  //     constexpr
  //     Integral(const Grad<ShapeFunction<FES, TrialSpace>>& gu, const Grad<ShapeFunction<FES, TestSpace>>& gv)
  //       : Integral(Dot(gu, gv))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     Math::Matrix getMatrix(const Geometry::Simplex& element) const override
  //     {
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       if (&trial == &test)
  //       {
  //         mfem::DenseMatrix mat;
  //         const int order =
  //           getIntegrationOrder(
  //               getIntegrand().getLHS().getFiniteElementSpace(),
  //               getIntegrand().getRHS().getFiniteElementSpace(),
  //               element);
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         mfem::ConstantCoefficient one(1.0);
  //         mfem::DiffusionIntegrator bfi(one);
  //         bfi.SetIntRule(ir);
  //         bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //         return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(const Grad<ShapeFunction<FES, TrialSpace>>&, const Grad<ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>;
  // template <class FES>
  // Integral(const Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>&)
  //   -> Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>;

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega (f u) \cdot v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
  //   : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;

  //     constexpr
  //     Integral(
  //         const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>& fu,
  //         const ShapeFunction<FES, TestSpace>& v)
  //       : Integral(Dot(fu, v))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Math::Matrix getMatrix(
  //         const Geometry::Simplex& element) const override
  //     {
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       const int order =
  //         getIntegrationOrder(
  //           getIntegrand().getLHS().getFiniteElementSpace(),
  //           getIntegrand().getRHS().getFiniteElementSpace(),
  //           element);

  //       if (&trial == &test)
  //       {
  //         mfem::DenseMatrix mat;
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
  //         switch (getIntegrand().getLHS().getLHS().getRangeType())
  //         {
  //           case RangeType::Scalar:
  //           {
  //             switch (getIntegrand().getLHS().getRHS().getRangeType())
  //             {
  //               case RangeType::Scalar:
  //               {
  //                 mfem::MassIntegrator bfi(q.template get<RangeType::Scalar>());
  //                 bfi.SetIntRule(ir);
  //                 bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //                 return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //               }
  //               case RangeType::Vector:
  //               {
  //                 mfem::VectorMassIntegrator bfi(q.template get<RangeType::Scalar>());
  //                 bfi.SetIntRule(ir);
  //                 bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //                 return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //               }
  //               case RangeType::Matrix:
  //               {
  //                 assert(false); // Unsupported
  //                 break;
  //               }
  //             }
  //             break;
  //           }
  //           case RangeType::Vector:
  //           {
  //             assert(false); // Unsupported
  //             break;
  //           }
  //           case RangeType::Matrix:
  //           {
  //             switch (getIntegrand().getLHS().getRHS().getRangeType())
  //             {
  //               case RangeType::Scalar:
  //               {
  //                 assert(false); // Unsupported
  //                 break;
  //               }
  //               case RangeType::Vector:
  //               {
  //                 mfem::VectorMassIntegrator bfi(q.template get<RangeType::Matrix>());
  //                 bfi.SetIntRule(ir);
  //                 bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //                 return Eigen::Map<Math::Matrix>(mat.Data(), mat.NumRows(), mat.NumCols());
  //               }
  //               case RangeType::Matrix:
  //               {
  //                 assert(false); // Unsupported
  //                 break;
  //               }
  //             }
  //             break;
  //           }
  //         }
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(
  //     const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&,
  //     const ShapeFunction<FES, TestSpace>&)
  //   -> Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
  // template <class FES>
  // Integral(
  //     const Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega (f \nabla u) \cdot \nabla v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<
  //   Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>,
  //   Grad<ShapeFunction<FES, TestSpace>>>>
  //     : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;

  //     constexpr
  //     Integral(
  //         const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>& fgu,
  //         const Grad<ShapeFunction<FES, TestSpace>>& gv)
  //       : Integral(Dot(fgu, gv))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual Math::Matrix getMatrix(
  //         const Geometry::Simplex& element) const override
  //     {
  //       mfem::DenseMatrix mat;
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       const int order =
  //         getIntegrationOrder(
  //           getIntegrand().getLHS().getFiniteElementSpace(),
  //           getIntegrand().getRHS().getFiniteElementSpace(),
  //           element);
  //       if (&trial == &test)
  //       {
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
  //         switch (getIntegrand().getLHS().getLHS().getRangeType())
  //         {
  //           case RangeType::Scalar:
  //           {
  //             mfem::DiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
  //             bfi.SetIntRule(ir);
  //             bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //             break;
  //           }
  //           case RangeType::Vector:
  //           {
  //             assert(false); // Unsupported
  //             break;
  //           }
  //           case RangeType::Matrix:
  //           {
  //             mfem::DiffusionIntegrator bfi(q.template get<RangeType::Matrix>());
  //             bfi.SetIntRule(ir);
  //             bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //             break;
  //           }
  //         }
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //       Math::Matrix res = Eigen::Map<Math::Matrix>(mat.GetData(), mat.NumRows(), mat.NumCols());
  //       return res;
  //     }

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(
  //     const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&,
  //     const Grad<ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
  // template <class FES>
  // Integral(
  //     const Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega (f \nabla u) \cdot \nabla v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<
  //   Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>,
  //   Jacobian<ShapeFunction<FES, TestSpace>>>>
  //     : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  // {
  //   public:
  //     using Parent =
  //       Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
  //     using Integrand =
  //       Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>;

  //     constexpr
  //     Integral(
  //         const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>& fgu,
  //         const Jacobian<ShapeFunction<FES, TestSpace>>& gv)
  //       : Integral(Dot(fgu, gv))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
  //             const Geometry::Simplex& element)
  //           {
  //             const auto& trial = trialFes.getFiniteElement(element);
  //             const auto& test = testFes.getFiniteElement(element);
  //             if (trial.Space() == mfem::FunctionSpace::Pk)
  //               return trial.GetOrder() + test.GetOrder() - 2;
  //             else
  //               return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     virtual Math::Matrix getMatrix(const Geometry::Simplex& element) const override
  //     {
  //       mfem::DenseMatrix mat;
  //       const auto& trial = getIntegrand().getLHS()
  //                              .getFiniteElementSpace()
  //                              .getFiniteElement(element);
  //       const auto& test = getIntegrand().getRHS()
  //                             .getFiniteElementSpace()
  //                             .getFiniteElement(element);
  //       const int order =
  //         getIntegrationOrder(
  //             getIntegrand().getLHS().getFiniteElementSpace(),
  //             getIntegrand().getRHS().getFiniteElementSpace(),
  //             element);
  //       if (&trial == &test)
  //       {
  //         const mfem::IntegrationRule* ir =
  //           trial.Space() == mfem::FunctionSpace::rQk ?
  //             &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
  //             &mfem::IntRules.Get(trial.GetGeomType(), order);
  //         auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
  //         switch (getIntegrand().getLHS().getLHS().getRangeType())
  //         {
  //           case RangeType::Scalar:
  //           {
  //             mfem::VectorDiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
  //             bfi.SetIntRule(ir);
  //             bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
  //             break;
  //           }
  //           case RangeType::Vector:
  //           {
  //             assert(false); // Unsupported
  //             break;
  //           }
  //           case RangeType::Matrix:
  //           {
  //             assert(false); // Unimplemented
  //             break;
  //           }
  //         }
  //       }
  //       else
  //       {
  //         assert(false); // Unimplemented
  //       }
  //       Math::Matrix res = Eigen::Map<Math::Matrix>(mat.GetData(), mat.NumRows(), mat.NumCols());
  //       return res;
  //     }

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(
  //     const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>&,
  //     const Jacobian<ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>>;
  // template <class FES>
  // Integral(
  //     const Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>&)
  //   -> Integral<Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>>;

  // /* <<-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
  //  * ----------------------------------------------------------------------||
  //  */

  // /* ||-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<ShapeFunctionBase<TestSpace>>
  //  * ---------------------------------------------------------------------->>
  //  */

  // /**
  //  * @ingroup IntegralSpecializations
  //  *
  //  * Optimized integration of the expression:
  //  * @f[
  //  *   \int_\Omega f \cdot v \ dx
  //  * @f]
  //  * where @f$ f @f$ is a function (scalar, vector or matrix valued).
  //  */
  // template <class FES>
  // class Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
  //   : public Integral<ShapeFunctionBase<TestSpace>>
  // {
  //   public:
  //     using IntegrationOrder =
  //       std::function<int(const FiniteElementSpaceBase&, const Geometry::Simplex&)>;
  //     using Parent    = Integral<ShapeFunctionBase<TestSpace>>;
  //     using Integrand  = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;

  //     constexpr
  //     Integral(const FunctionBase& f, const ShapeFunction<FES, TestSpace>& v)
  //       : Integral(Dot(f, v))
  //     {}

  //     constexpr
  //     Integral(const Integrand& integrand)
  //       : Parent(integrand)
  //     {
  //       setIntegrationOrder(
  //           [](const FiniteElementSpaceBase& fes, const Geometry::Simplex& element)
  //           {
  //             return 2 * fes.getFiniteElement(element).GetOrder();
  //           });
  //     }

  //     constexpr
  //     Integral(const Integral& other)
  //       : Parent(other)
  //     {}

  //     constexpr
  //     Integral(Integral&& other)
  //       : Parent(std::move(other))
  //     {}

  //     Integral& setIntegrationOrder(IntegrationOrder order)
  //     {
  //       return static_cast<Integral&>(Parent::setIntegrationOrder(order));
  //     }

  //     int getIntegrationOrder(
  //         const FiniteElementSpaceBase& fes, const Geometry::Simplex& element) const
  //     {
  //       return Parent::getIntegrationOrder(fes, element);
  //     }

  //     virtual const Integrand& getIntegrand() const override
  //     {
  //       return static_cast<const Integrand&>(Parent::getIntegrand());
  //     }

  //     virtual Math::Vector getVector(const Geometry::Simplex& element) const override
  //     {
  //       const FunctionBase& f = getIntegrand().getLHS();

  //       const auto& fe = getIntegrand().getFiniteElementSpace()
  //                            .getFiniteElement(element);
  //       const mfem::IntegrationRule *ir =
  //         &mfem::IntRules.Get(
  //             fe.GetGeomType(),
  //             getIntegrationOrder(getIntegrand().getFiniteElementSpace(), element));
  //       auto q = f.build(element.getMesh());

  //       mfem::Vector vec;
  //       switch (f.getRangeType())
  //       {
  //         case RangeType::Scalar:
  //         {
  //           mfem::DomainLFIntegrator lfi(q.get<RangeType::Scalar>());
  //           lfi.SetIntRule(ir);
  //           lfi.AssembleRHSElementVect(fe, element.getTransformation(), vec);
  //           Math::Vector res = Eigen::Map<Math::Vector>(vec.GetData(), vec.Size());
  //           return res;
  //         }
  //         case RangeType::Vector:
  //         {
  //           mfem::VectorDomainLFIntegrator lfi(q.get<RangeType::Vector>());
  //           lfi.SetIntRule(ir);
  //           lfi.AssembleRHSElementVect(fe, element.getTransformation(), vec);
  //           Math::Vector res = Eigen::Map<Math::Vector>(vec.GetData(), vec.Size());
  //           return res;
  //         }
  //         case RangeType::Matrix:
  //         {
  //           assert(false); // Unsupported
  //           Math::Vector res = Eigen::Map<Math::Vector>(nullptr, 0);
  //           return res;
  //         }
  //       }
  //     }

  //     virtual Integral* copy() const noexcept override
  //     {
  //       return new Integral(*this);
  //     }
  // };
  // template <class FES>
  // Integral(const FunctionBase&, const ShapeFunction<FES, TestSpace>&)
  //   -> Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
  // template <class FES>
  // Integral(const Dot<FunctionBase, ShapeFunction<FES, TestSpace>>&)
  //   -> Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;

  // /* <<-- OPTIMIZATIONS -----------------------------------------------------
  //  * Integral<ShapeFunctionBase<TestSpace>>
  //  * ----------------------------------------------------------------------||
  //  */
}

#endif
