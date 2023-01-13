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
#include "LinearForm.h"
#include "ForwardDecls.h"
#include "GridFunction.h"
#include "Function.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "MatrixFunction.h"
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
    * @defgroup BoundaryIntegralSpecializations BoundaryIntegral Template Specializations
    * @brief Template specializations of the BoundaryIntegral class.
    *
    * @see BoundaryIntegral
    */

   /**
    * @defgroup InterfaceIntegralSpecializations InterfaceIntegral Template Specializations
    * @brief Template specializations of the InterfaceIntegral class.
    *
    * @see InterfaceIntegral
    */

   /**
    * @defgroup BoundaryFaceIntegralSpecializations BoundaryFaceIntegral Template Specializations
    * @brief Template specializations of the BoundaryFaceIntegral class.
    *
    * @see BoundaryFaceIntegral
    */

   /**
    * @ingroup IntegralSpecializations
    * @brief Integration of the dot product of a trial and test operators.
    *
    * Given two operators defined over trial and test spaces @f$ U_h
    * @f$ and @f$ V_h @f$,
    * @f[
    *    A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
    * @f]
    * this class represents the integral of their dot product:
    * @f[
    *    \int_{\mathcal{T}_h} A(u) : B(v) \ dx
    * @f]
    */
   template <>
   class Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
      : public BilinearFormIntegratorBase
   {
      public:
         using IntegrationOrderFunction =
            std::function<
               int(const FiniteElementSpaceBase&, const FiniteElementSpaceBase&, const Geometry::Simplex&)>;
         using Parent = BilinearFormIntegratorBase;
         using Integrand =
            Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;

         /**
          * @brief Integral of the dot product of trial and test operators
          *
          * Constructs an instance representing the following integral:
          * @f[
          *    \int_\Omega A(u) : B(v) \ dx
          * @f]
          *
          * @param[in] lhs Trial operator @f$ A(u) @f$
          * @param[in] rhs Test operator @f$ B(v) @f$
          */
         Integral(const ShapeFunctionBase<TrialSpace>& lhs, const ShapeFunctionBase<TestSpace>& rhs)
            : Integral(Dot(lhs, rhs))
         {}

         /**
          * @brief Integral of the dot product of trial and test operators
          *
          * Constructs the following object representing the following
          * integral:
          * @f[
          *    \int_\Omega A(u) : B(v) \ dx
          * @f]
          *
          * @param[in] prod Dot product instance
          */
         Integral(const Integrand& prod)
            :  BilinearFormIntegratorBase(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
               m_prod(prod),
               m_intOrder(
                  [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
                     const Geometry::Simplex& element)
                  {
                     const auto& trial = trialFes.getFiniteElement(element);
                     const auto& test = testFes.getFiniteElement(element);
                     return trial.GetOrder() + test.GetOrder() + element.getTransformation().OrderW();
                  })
         {}

         Integral(const Integral& other)
            : BilinearFormIntegratorBase(other),
              m_prod(other.m_prod),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : BilinearFormIntegratorBase(std::move(other)),
              m_prod(std::move(other.m_prod)),
              m_intOrder(std::move(other.m_intOrder))
         {}

         /**
          * @brief Sets the function which calculates the integration order
          * @param[in] order Function which computes the order of integration
          * @returns Reference to self (for method chaining)
          */
         Integral& setIntegrationOrder(IntegrationOrderFunction order)
         {
            m_intOrder = order;
            return *this;
         }

         int getIntegrationOrder(
               const FiniteElementSpaceBase& trialFes,
               const FiniteElementSpaceBase& testFes,
               const Geometry::Simplex& element) const
         {
            return m_intOrder(trialFes, testFes, element);
         }

         virtual const Integrand& getIntegrand() const
         {
            return m_prod;
         }

         virtual Region getRegion() const override
         {
            return Region::Domain;
         }

         virtual mfem::DenseMatrix getMatrix(
               const Geometry::Simplex& element) const override;

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         Integrand m_prod;
         IntegrationOrderFunction m_intOrder;
   };
   Integral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   Integral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
      -> Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

   /**
    * @ingroup BoundaryIntegralSpecializations
    * @brief Boundary integration of the dot product of a trial and test operators.
    */
   template <>
   class BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
      : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent    = Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
         using Parent::Parent;
         Region getRegion() const override { return Region::Boundary; }
         BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
   };
   BoundaryIntegral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   BoundaryIntegral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

   /**
    * @ingroup InterfaceIntegralSpecializations
    * @brief Interface integration of the dot product of a trial and test operators.
    *
    * Given two operators defined over trial and test spaces @f$ U_h
    * @f$ and @f$ V_h @f$,
    * @f[
    *    A : U_h \rightarrow \mathbb{R}^{p \times q}, \quad B : V_h \rightarrow \mathbb{R}^{p \times q},
    * @f]
    * this class represents the integral of their dot product:
    * @f[
    *    \int_{\mathcal{I}_h} A(u) : B(v) \ dx
    * @f]
    * over the interface @f$ \mathcal{I}_h @f$ of the triangulation @f$
    * \mathcal{T}_h @f$.
    */
   template <>
   class InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
      : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent    = Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
         using Parent::Parent;
         Region getRegion() const override { return Region::Interface; }
         InterfaceIntegral* copy() const noexcept override { return new InterfaceIntegral(*this); }
   };
   InterfaceIntegral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   InterfaceIntegral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
      -> InterfaceIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

   /**
    * @ingroup IntegralSpecializations
    * @brief Integration of a test operator.
    *
    * Given an operator defined over a test space @f$ V_h @f$
    * @f[
    *    A : V_h \rightarrow \mathbb{R},
    * @f]
    * this class will represent its integral
    * @f[
    *    \int_\Omega A(v) \ dx \ .
    * @f]
    */
   template <>
   class Integral<ShapeFunctionBase<TestSpace>> : public LinearFormIntegratorBase
   {
      public:
         using IntegrationOrder =
            std::function<int(const FiniteElementSpaceBase&, const Geometry::Simplex&)>;
         using Parent = LinearFormIntegratorBase;
         using Integrand = ShapeFunctionBase<TestSpace>;

         Integral(const FunctionBase& lhs, const ShapeFunctionBase<TestSpace>& rhs)
            : Integral(Dot(lhs, rhs))
         {}

         /**
          * @brief Integral of a scalar valued test operator
          *
          * Given
          * @f[
          *    A : V_h \rightarrow \mathbb{R}
          * @f]
          * constructs an instance representing the following integral
          * @f[
          *    \int_\Omega A(v) \ dx \ .
          * @f]
          */
         Integral(const Integrand& integrand)
            :  LinearFormIntegratorBase(integrand.getLeaf()),
               m_integrand(integrand.copy()),
               m_intOrder(
                     [](const FiniteElementSpaceBase& fes, const Geometry::Simplex& element)
                     {
                        const auto& fe = fes.getFiniteElement(element);
                        return fe.GetOrder() + element.getTransformation().OrderW();
                     })
         {}

         Integral(const Integral& other)
            : LinearFormIntegratorBase(other),
              m_integrand(other.m_integrand->copy()),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : LinearFormIntegratorBase(std::move(other)),
              m_integrand(std::move(other.m_integrand))
         {}

         Integral& setIntegrationOrder(IntegrationOrder order)
         {
            m_intOrder = order;
            return *this;
         }

         int getIntegrationOrder(const FiniteElementSpaceBase& fes, const Geometry::Simplex& element) const
         {
            return m_intOrder(fes, element);
         }

         virtual const Integrand& getIntegrand() const
         {
            assert(m_integrand);
            return *m_integrand;
         }

         virtual Region getRegion() const override
         {
            return Region::Domain;
         }

         virtual mfem::Vector getVector(
               const Geometry::Simplex& element) const override;

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }

      private:
         std::unique_ptr<Integrand> m_integrand;
         IntegrationOrder m_intOrder;
   };
   Integral(const FunctionBase&, const ShapeFunctionBase<TestSpace>&) -> Integral<ShapeFunctionBase<TestSpace>>;
   Integral(const ShapeFunctionBase<TestSpace>&) -> Integral<ShapeFunctionBase<TestSpace>>;

   template <>
   class FaceIntegral<ShapeFunctionBase<TestSpace>> : public Integral<ShapeFunctionBase<TestSpace>>
   {
      public:
         using Parent    = Integral<ShapeFunctionBase<TestSpace>>;
         using Integrand = ShapeFunctionBase<TestSpace>;
         using Parent::Parent;
         Region getRegion() const override { return Region::Faces; }
         FaceIntegral* copy() const noexcept override { return new FaceIntegral(*this); }
   };
   FaceIntegral(const FunctionBase&, const ShapeFunctionBase<TestSpace>&)
      -> FaceIntegral<ShapeFunctionBase<TestSpace>>;
   FaceIntegral(const ShapeFunctionBase<TestSpace>&)
      -> FaceIntegral<ShapeFunctionBase<TestSpace>>;

   template <>
   class BoundaryIntegral<ShapeFunctionBase<TestSpace>> : public Integral<ShapeFunctionBase<TestSpace>>
   {
      public:
         using Parent    = Integral<ShapeFunctionBase<TestSpace>>;
         using Integrand = ShapeFunctionBase<TestSpace>;
         using Parent::Parent;
         Region getRegion() const override { return Region::Boundary; }
         BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
   };
   BoundaryIntegral(const FunctionBase&, const ShapeFunctionBase<TestSpace>&)
      -> BoundaryIntegral<ShapeFunctionBase<TestSpace>>;
   BoundaryIntegral(const ShapeFunctionBase<TestSpace>&)
      -> BoundaryIntegral<ShapeFunctionBase<TestSpace>>;

   /**
    * @ingroup IntegralSpecializations
    * @brief Integration of a GridFunction object.
    */
   template <class FES>
   class Integral<GridFunction<FES>> : public FormLanguage::Base
   {
      public:
         using Parent    = FormLanguage::Base;
         using Integrand = GridFunction<FES>&;

         /**
          * @brief Constructs the integral object
          */
         Integral(GridFunction<FES>& u)
            : m_u(u),
              m_v(u.getFiniteElementSpace()),
              m_one(u.getFiniteElementSpace()),
              m_lf(m_v),
              m_assembled(false)
         {
            assert(u.getFiniteElementSpace().getVectorDimension() == 1);
            m_one = ScalarFunction(1.0);
            m_lf.from(Integral<ShapeFunctionBase<TestSpace>>(ScalarFunction(u) * m_v));
         }

         Integral(const Integral& other)
            : Integral(other.m_u)
         {}

         Integral(Integral&& other) = default;

         /**
          * @brief Integrates the expression and returns the value
          * @returns Value of integral
          *
          * This method does not cache the integrated value.
          */
         double compute()
         {
            m_lf.assemble();
            m_assembled = true;
            return m_lf(m_one);
         }

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         TestFunction<FES>                m_v;
         GridFunction<FES>&               m_u;
         GridFunction<FES>                m_one;
         LinearForm<FES, Context::Serial, mfem::Vector>    m_lf;
         bool m_assembled;
   };
   template <class FES>
   Integral(GridFunction<FES>&) -> Integral<GridFunction<FES>>;

   /* ||-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @ingroup IntegralSpecializations
    *
    * @f[
    * \int_\Omega \nabla u \cdot \nabla v \ dx
    * @f]
    */
   template <class FES>
   class Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>
      : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent =
            Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand =
            Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>;

         constexpr
         Integral(const Grad<ShapeFunction<FES, TrialSpace>>& gu, const Grad<ShapeFunction<FES, TestSpace>>& gv)
            : Integral(Dot(gu, gv))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Parent(integrand)
         {
            setIntegrationOrder(
                  [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
                     const Geometry::Simplex& element)
                  {
                     const auto& trial = trialFes.getFiniteElement(element);
                     const auto& test = testFes.getFiniteElement(element);
                     if (trial.Space() == mfem::FunctionSpace::Pk)
                        return trial.GetOrder() + test.GetOrder() - 2;
                     else
                        return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
                  });
         }

         constexpr
         Integral(const Integral& other)
            : Parent(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Parent(std::move(other))
         {}

         const Integrand& getIntegrand() const override
         {
            return static_cast<const Integrand&>(Parent::getIntegrand());
         }

         mfem::DenseMatrix getMatrix(const Geometry::Simplex& element) const override
         {
            const auto& trial = getIntegrand().getLHS()
                                              .getFiniteElementSpace()
                                              .getFiniteElement(element);
            const auto& test = getIntegrand().getRHS()
                                             .getFiniteElementSpace()
                                             .getFiniteElement(element);
            if (&trial == &test)
            {
               mfem::DenseMatrix mat;
               const int order =
                  getIntegrationOrder(
                        getIntegrand().getLHS().getFiniteElementSpace(),
                        getIntegrand().getRHS().getFiniteElementSpace(),
                        element);
               const mfem::IntegrationRule* ir =
                  trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(trial.GetGeomType(), order);
               mfem::ConstantCoefficient one(1.0);
               mfem::DiffusionIntegrator bfi(one);
               bfi.SetIntRule(ir);
               bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
               return mat;
            }
            else
            {
               assert(false); // Unimplemented
            }
         }

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
   template <class FES>
   Integral(const Grad<ShapeFunction<FES, TrialSpace>>&, const Grad<ShapeFunction<FES, TestSpace>>&)
      -> Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>;
   template <class FES>
   Integral(const Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>&)
      -> Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>;

   /**
    * @ingroup IntegralSpecializations
    *
    * Optimized integration of the expression:
    * @f[
    *    \int_\Omega (f u) \cdot v \ dx
    * @f]
    * where @f$ f @f$ is a function (scalar or matrix valued).
    */
   template <class FES>
   class Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
      : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent =
            Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand =
            Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;

         constexpr
         Integral(
               const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>& fu,
               const ShapeFunction<FES, TestSpace>& v)
            : Integral(Dot(fu, v))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Parent(integrand)
         {
            setIntegrationOrder(
                  [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
                     const Geometry::Simplex& element)
                  {
                     const auto& trial = trialFes.getFiniteElement(element);
                     const auto& test = testFes.getFiniteElement(element);
                     if (trial.Space() == mfem::FunctionSpace::Pk)
                        return trial.GetOrder() + test.GetOrder() - 2;
                     else
                        return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
                  });
         }

         constexpr
         Integral(const Integral& other)
            : Parent(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Parent(std::move(other))
         {}

         virtual const Integrand& getIntegrand() const override
         {
            return static_cast<const Integrand&>(Parent::getIntegrand());
         }

         virtual mfem::DenseMatrix getMatrix(
               const Geometry::Simplex& element) const override
         {
            const auto& trial = getIntegrand().getLHS()
                                              .getFiniteElementSpace()
                                              .getFiniteElement(element);
            const auto& test = getIntegrand().getRHS()
                                             .getFiniteElementSpace()
                                             .getFiniteElement(element);
            const int order =
               getIntegrationOrder(
                  getIntegrand().getLHS().getFiniteElementSpace(),
                  getIntegrand().getRHS().getFiniteElementSpace(),
                  element);

            if (&trial == &test)
            {
               mfem::DenseMatrix mat;
               const mfem::IntegrationRule* ir =
                  trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(trial.GetGeomType(), order);
               auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
               switch (getIntegrand().getLHS().getLHS().getRangeType())
               {
                  case RangeType::Scalar:
                  {
                     switch (getIntegrand().getLHS().getRHS().getRangeType())
                     {
                        case RangeType::Scalar:
                        {
                           mfem::MassIntegrator bfi(q.template get<RangeType::Scalar>());
                           bfi.SetIntRule(ir);
                           bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
                           break;
                        }
                        case RangeType::Vector:
                        {
                           mfem::VectorMassIntegrator bfi(q.template get<RangeType::Scalar>());
                           bfi.SetIntRule(ir);
                           bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
                           break;
                        }
                        case RangeType::Matrix:
                        {
                           assert(false); // Unsupported
                           break;
                        }
                     }
                     break;
                  }
                  case RangeType::Vector:
                  {
                     assert(false); // Unsupported
                     break;
                  }
                  case RangeType::Matrix:
                  {
                     switch (getIntegrand().getLHS().getRHS().getRangeType())
                     {
                        case RangeType::Scalar:
                        {
                           assert(false); // Unsupported
                           break;
                        }
                        case RangeType::Vector:
                        {
                           mfem::VectorMassIntegrator bfi(q.template get<RangeType::Matrix>());
                           bfi.SetIntRule(ir);
                           bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
                           break;
                        }
                        case RangeType::Matrix:
                        {
                           assert(false); // Unsupported
                           break;
                        }
                     }
                     break;
                  }
               }
               return mat;
            }
            else
            {
               assert(false); // Unimplemented
            }
         }

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
   template <class FES>
   Integral(
         const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&,
         const ShapeFunction<FES, TestSpace>&)
      -> Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
   template <class FES>
   Integral(
         const Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>&)
      -> Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;

   template <class FES>
   class BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
      : public Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>
   {
      public:
         using Parent =
            Integral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
         using Integrand =
            Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>;
         using Parent::Parent;
         Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
         BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
   };
   template <class FES>
   BoundaryIntegral(
         const Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>&,
         const ShapeFunction<FES, TestSpace>&)
      -> BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;
   template <class FES>
   BoundaryIntegral(
         const Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>&)
      -> BoundaryIntegral<Dot<Mult<FunctionBase, ShapeFunction<FES, TrialSpace>>, ShapeFunction<FES, TestSpace>>>;


   /**
    * @ingroup IntegralSpecializations
    *
    * Optimized integration of the expression:
    * @f[
    *    \int_\Omega (f \nabla u) \cdot \nabla v \ dx
    * @f]
    * where @f$ f @f$ is a function (scalar or matrix valued).
    */
   template <class FES>
   class Integral<Dot<
      Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>,
      Grad<ShapeFunction<FES, TestSpace>>>>
         : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent =
            Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand =
            Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;

         constexpr
         Integral(
               const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>& fgu,
               const Grad<ShapeFunction<FES, TestSpace>>& gv)
            : Integral(Dot(fgu, gv))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Parent(integrand)
         {
            setIntegrationOrder(
                  [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
                     const Geometry::Simplex& element)
                  {
                     const auto& trial = trialFes.getFiniteElement(element);
                     const auto& test = testFes.getFiniteElement(element);
                     if (trial.Space() == mfem::FunctionSpace::Pk)
                        return trial.GetOrder() + test.GetOrder() - 2;
                     else
                        return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
                  });
         }

         constexpr
         Integral(const Integral& other)
            : Parent(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Parent(std::move(other))
         {}

         virtual mfem::DenseMatrix getMatrix(
               const Geometry::Simplex& element) const override
         {
            mfem::DenseMatrix mat;
            const auto& trial = getIntegrand().getLHS()
                                              .getFiniteElementSpace()
                                              .getFiniteElement(element);
            const auto& test = getIntegrand().getRHS()
                                             .getFiniteElementSpace()
                                             .getFiniteElement(element);
            const int order =
               getIntegrationOrder(
                  getIntegrand().getLHS().getFiniteElementSpace(),
                  getIntegrand().getRHS().getFiniteElementSpace(),
                  element);
            if (&trial == &test)
            {
               const mfem::IntegrationRule* ir =
                  trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(trial.GetGeomType(), order);
               auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
               switch (getIntegrand().getLHS().getLHS().getRangeType())
               {
                  case RangeType::Scalar:
                  {
                     mfem::DiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
                     bfi.SetIntRule(ir);
                     bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
                     break;
                  }
                  case RangeType::Vector:
                  {
                     assert(false); // Unsupported
                     break;
                  }
                  case RangeType::Matrix:
                  {
                     mfem::DiffusionIntegrator bfi(q.template get<RangeType::Matrix>());
                     bfi.SetIntRule(ir);
                     bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
                     break;
                  }
               }
            }
            else
            {
               assert(false); // Unimplemented
            }
         }

         virtual const Integrand& getIntegrand() const override
         {
            return static_cast<const Integrand&>(Parent::getIntegrand());
         }

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
   template <class FES>
   Integral(
         const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&,
         const Grad<ShapeFunction<FES, TestSpace>>&)
      -> Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
   template <class FES>
   Integral(
         const Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>&)
      -> Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;


   template <class FES>
   class BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>
      : public Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>
   {
      public:
         using Parent =
            Integral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
         using Integrand =
            Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>;
         using Parent::Parent;
         Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
         BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
   };
   template <class FES>
   BoundaryIntegral(
         const Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>&,
         const Grad<ShapeFunction<FES, TestSpace>>&)
      -> BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;
   template <class FES>
   BoundaryIntegral(
         const Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>&)
      -> BoundaryIntegral<Dot<Mult<FunctionBase, Grad<ShapeFunction<FES, TrialSpace>>>, Grad<ShapeFunction<FES, TestSpace>>>>;


   /**
    * @ingroup IntegralSpecializations
    *
    * Optimized integration of the expression:
    * @f[
    *    \int_\Omega (f \nabla u) \cdot \nabla v \ dx
    * @f]
    * where @f$ f @f$ is a function (scalar or matrix valued).
    */
   template <class FES>
   class Integral<Dot<
      Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>,
      Jacobian<ShapeFunction<FES, TestSpace>>>>
         : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent =
            Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand =
            Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>;

         constexpr
         Integral(
               const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>& fgu,
               const Jacobian<ShapeFunction<FES, TestSpace>>& gv)
            : Integral(Dot(fgu, gv))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Parent(integrand)
         {
            setIntegrationOrder(
                  [](const FiniteElementSpaceBase& trialFes, const FiniteElementSpaceBase& testFes,
                     const Geometry::Simplex& element)
                  {
                     const auto& trial = trialFes.getFiniteElement(element);
                     const auto& test = testFes.getFiniteElement(element);
                     if (trial.Space() == mfem::FunctionSpace::Pk)
                        return trial.GetOrder() + test.GetOrder() - 2;
                     else
                        return trial.GetOrder() + test.GetOrder() + trial.GetDim() - 1;
                  });
         }

         constexpr
         Integral(const Integral& other)
            : Parent(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Parent(std::move(other))
         {}

         virtual mfem::DenseMatrix getMatrix(
               const Geometry::Simplex& element) const override
         {
            mfem::DenseMatrix mat;
            const auto& trial = getIntegrand().getLHS()
                                              .getFiniteElementSpace()
                                              .getFiniteElement(element);
            const auto& test = getIntegrand().getRHS()
                                             .getFiniteElementSpace()
                                             .getFiniteElement(element);
            const int order =
               getIntegrationOrder(
                     getIntegrand().getLHS().getFiniteElementSpace(),
                     getIntegrand().getRHS().getFiniteElementSpace(),
                     element);
            if (&trial == &test)
            {
               const mfem::IntegrationRule* ir =
                  trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(trial.GetGeomType(), order);
               auto q = getIntegrand().getLHS().getLHS().build(element.getMesh());
               switch (getIntegrand().getLHS().getLHS().getRangeType())
               {
                  case RangeType::Scalar:
                  {
                     mfem::VectorDiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
                     bfi.SetIntRule(ir);
                     bfi.AssembleElementMatrix(trial, element.getTransformation(), mat);
                     break;
                  }
                  case RangeType::Vector:
                  {
                     assert(false); // Unsupported
                     break;
                  }
                  case RangeType::Matrix:
                  {
                     assert(false); // Unimplemented
                     break;
                  }
               }
            }
            else
            {
               assert(false); // Unimplemented
            }
            return mat;
         }

         virtual const Integrand& getIntegrand() const override
         {
            return static_cast<const Integrand&>(Parent::getIntegrand());
         }

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
   template <class FES>
   Integral(
         const Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>&,
         const Jacobian<ShapeFunction<FES, TestSpace>>&)
      -> Integral<Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>>;
   template <class FES>
   Integral(
         const Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>&)
      -> Integral<Dot<Mult<FunctionBase, Jacobian<ShapeFunction<FES, TrialSpace>>>, Jacobian<ShapeFunction<FES, TestSpace>>>>;

   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
    * ----------------------------------------------------------------------||
    */

   /* ||-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<ShapeFunctionBase<TestSpace>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @ingroup IntegralSpecializations
    *
    * Optimized integration of the expression:
    * @f[
    *    \int_\Omega f \cdot v \ dx
    * @f]
    * where @f$ f @f$ is a function (scalar, vector or matrix valued).
    */
   template <class FES>
   class Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
      : public Integral<ShapeFunctionBase<TestSpace>>
   {
      public:
         using IntegrationOrder =
            std::function<int(const FiniteElementSpaceBase&, const Geometry::Simplex&)>;
         using Parent      = Integral<ShapeFunctionBase<TestSpace>>;
         using Integrand   = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;

         constexpr
         Integral(const FunctionBase& f, const ShapeFunction<FES, TestSpace>& v)
            : Integral(Dot(f, v))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Parent(integrand)
         {
            setIntegrationOrder(
                  [](const FiniteElementSpaceBase& fes, const Geometry::Simplex& element)
                  {
                     return 2 * fes.getFiniteElement(element).GetOrder();
                  });
         }

         constexpr
         Integral(const Integral& other)
            : Parent(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Parent(std::move(other))
         {}

         Integral& setIntegrationOrder(IntegrationOrder order)
         {
            return static_cast<Integral&>(Parent::setIntegrationOrder(order));
         }

         int getIntegrationOrder(
               const FiniteElementSpaceBase& fes, const Geometry::Simplex& element) const
         {
            return Parent::getIntegrationOrder(fes, element);
         }

         virtual const Integrand& getIntegrand() const override
         {
            return static_cast<const Integrand&>(Parent::getIntegrand());
         }

         virtual mfem::Vector getVector(const Geometry::Simplex& element) const override
         {
            const FunctionBase& f = getIntegrand().getLHS();

            const auto& fe = getIntegrand().getFiniteElementSpace()
                                           .getFiniteElement(element);
            const mfem::IntegrationRule *ir =
               &mfem::IntRules.Get(
                     fe.GetGeomType(),
                     getIntegrationOrder(getIntegrand().getFiniteElementSpace(), element));
            auto q = f.build(element.getMesh());

            mfem::Vector vec;
            switch (f.getRangeType())
            {
               case RangeType::Scalar:
               {
                  mfem::DomainLFIntegrator lfi(q.get<RangeType::Scalar>());
                  lfi.SetIntRule(ir);
                  lfi.AssembleRHSElementVect(fe, element.getTransformation(), vec);
                  break;
               }
               case RangeType::Vector:
               {
                  mfem::VectorDomainLFIntegrator lfi(q.get<RangeType::Vector>());
                  lfi.SetIntRule(ir);
                  lfi.AssembleRHSElementVect(fe, element.getTransformation(), vec);
                  break;
               }
               case RangeType::Matrix:
               {
                  assert(false); // Unsupported
                  break;
               }
            }
            return vec;
         }

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
   template <class FES>
   Integral(const FunctionBase&, const ShapeFunction<FES, TestSpace>&)
      -> Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
   template <class FES>
   Integral(const Dot<FunctionBase, ShapeFunction<FES, TestSpace>>&)
      -> Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;

   template <class FES>
   class BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
      : public Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>
   {
      public:
         using Parent      = Integral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
         using Integrand   = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;
         using Parent::Parent;
         Integrator::Region getRegion() const override { return Integrator::Region::Boundary; }
         BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
   };
   template <class FES>
   BoundaryIntegral(const FunctionBase&, const ShapeFunction<FES, TestSpace>&)
      -> BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;
   template <class FES>
   BoundaryIntegral(const Dot<FunctionBase, ShapeFunction<FES, TestSpace>>&)
      -> BoundaryIntegral<Dot<FunctionBase, ShapeFunction<FES, TestSpace>>>;

   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<ShapeFunctionBase<TestSpace>>
    * ----------------------------------------------------------------------||
    */
}

#endif
