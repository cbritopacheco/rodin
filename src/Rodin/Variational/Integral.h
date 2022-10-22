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
    * @brief This module contains information about the template
    * specializations of the Integral class.
    *
    * Throughout this section we utilize the symbols @f$ u \in U_h @f$ for the
    * trial function, and @f$ v \in V_h @f$ for the test function.
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
    *    \int_\Omega A(u) : B(v) \ dx
    * @f]
    */
   template <>
   class Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
      : public BilinearFormIntegratorBase
   {
      public:
         using Parent =
            BilinearFormDomainIntegrator;
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
            : BilinearFormIntegratorBase(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
              m_prod(prod),
              m_intOrder(
                    [](const Bilinear::Assembly::Common& as)
                    { return as.trial.GetOrder() + as.test.GetOrder() + as.trans.OrderW(); })
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
         Integral& setIntegrationOrder(
            std::function<int(const Bilinear::Assembly::Common&)> order)
         {
            m_intOrder = order;
            return *this;
         }

         int getIntegrationOrder(const Bilinear::Assembly::Common& as) const
         {
            return m_intOrder(as);
         }

         virtual const Integrand& getIntegrand() const
         {
            return m_prod;
         }

         virtual IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         virtual void getElementMatrix(const Bilinear::Assembly::Common& as) const override;

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         Integrand m_prod;
         std::function<int(const Bilinear::Assembly::Common&)> m_intOrder;
   };
   Integral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   Integral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
      -> Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

   template <>
   class BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
      : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent    = Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
         using Parent::Parent;
         IntegratorRegion getIntegratorRegion() const override { return IntegratorRegion::Boundary; }
         BoundaryIntegral* copy() const noexcept override { return new BoundaryIntegral(*this); }
   };
   BoundaryIntegral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   BoundaryIntegral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

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
         using Parent = LinearFormDomainIntegrator;
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
                     [](const Linear::Assembly::Common& as)
                     { return as.fe.GetOrder() + as.trans.OrderW(); })
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

         Integral& setIntegrationOrder(std::function<int(const Linear::Assembly::Common&)> order)
         {
            m_intOrder = order;
            return *this;
         }

         int getIntegrationOrder(const Linear::Assembly::Common& as) const
         {
            return m_intOrder(as);
         }

         virtual const Integrand& getIntegrand() const
         {
            assert(m_integrand);
            return *m_integrand;
         }

         virtual IntegratorRegion getIntegratorRegion() const override
         {
            return IntegratorRegion::Domain;
         }

         virtual void getElementVector(const Linear::Assembly::Common& as) const override;

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }

      private:
         std::unique_ptr<Integrand> m_integrand;
         std::function<int(const Linear::Assembly::Common&)> m_intOrder;
   };
   Integral(const FunctionBase&, const ShapeFunctionBase<TestSpace>&) -> Integral<ShapeFunctionBase<TestSpace>>;
   Integral(const ShapeFunctionBase<TestSpace>&) -> Integral<ShapeFunctionBase<TestSpace>>;

   template <>
   class BoundaryIntegral<ShapeFunctionBase<TestSpace>> : public Integral<ShapeFunctionBase<TestSpace>>
   {
      public:
         using Parent    = Integral<ShapeFunctionBase<TestSpace>>;
         using Integrand = ShapeFunctionBase<TestSpace>;
         using Parent::Parent;
         IntegratorRegion getIntegratorRegion() const override { return IntegratorRegion::Boundary; }
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
            if (m_assembled)
               m_lf.update();
            else
               m_lf.assemble();
            m_assembled = true;
            return m_lf(m_one);
         }

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         TestFunction<FES>    m_v;
         GridFunction<FES>&   m_u;
         GridFunction<FES>    m_one;
         LinearForm<FES>      m_lf;
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
                  [](const Bilinear::Assembly::Common& as)
                  {
                     if (as.trial.Space() == mfem::FunctionSpace::Pk)
                        return as.trial.GetOrder() + as.test.GetOrder() - 2;
                     else
                        return as.trial.GetOrder() + as.test.GetOrder() + as.trial.GetDim() - 1;
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

         void getElementMatrix(const Bilinear::Assembly::Common& as) const override
         {
            if (&as.trial == &as.test)
            {
               const int order = getIntegrationOrder(as);
               const mfem::IntegrationRule* ir =
                  as.trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(as.trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(as.trial.GetGeomType(), order);
               mfem::ConstantCoefficient one(1.0);
               mfem::DiffusionIntegrator bfi(one);
               bfi.SetIntRule(ir);
               bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
            }
            else
            {
               assert(false); // Unimplemented
            }
         }

         Integral* copy() const noexcept override
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
                  [](const Bilinear::Assembly::Common& as)
                  {
                     if (as.trial.Space() == mfem::FunctionSpace::Pk)
                        return as.trial.GetOrder() + as.test.GetOrder() - 2;
                     else
                        return as.trial.GetOrder() + as.test.GetOrder() + as.trial.GetDim() - 1;
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

         virtual void getElementMatrix(const Bilinear::Assembly::Common& as) const override
         {
            const int order = getIntegrationOrder(as);
            if (&as.trial == &as.test)
            {
               const mfem::IntegrationRule* ir =
                  as.trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(as.trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(as.trial.GetGeomType(), order);
               auto q = getIntegrand().getLHS().getLHS().build();
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
                           bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
                           break;
                        }
                        case RangeType::Vector:
                        {
                           mfem::VectorMassIntegrator bfi(q.template get<RangeType::Scalar>());
                           bfi.SetIntRule(ir);
                           bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
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
                           bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
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
         IntegratorRegion getIntegratorRegion() const override { return IntegratorRegion::Boundary; }
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
                  [](const Bilinear::Assembly::Common& as)
                  {
                     if (as.trial.Space() == mfem::FunctionSpace::Pk)
                        return as.trial.GetOrder() + as.test.GetOrder() - 2;
                     else
                        return as.trial.GetOrder() + as.test.GetOrder() + as.trial.GetDim() - 1;
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

         virtual void getElementMatrix(const Bilinear::Assembly::Common& as) const override
         {
            const int order = getIntegrationOrder(as);
            if (&as.trial == &as.test)
            {
               const mfem::IntegrationRule* ir =
                  as.trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(as.trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(as.trial.GetGeomType(), order);
               auto q = getIntegrand().getLHS().getLHS().build();
               switch (getIntegrand().getLHS().getLHS().getRangeType())
               {
                  case RangeType::Scalar:
                  {
                     mfem::DiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
                     bfi.SetIntRule(ir);
                     bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
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
                     bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
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
         IntegratorRegion getIntegratorRegion() const override { return IntegratorRegion::Boundary; }
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
                  [](const Bilinear::Assembly::Common& as)
                  {
                     if (as.trial.Space() == mfem::FunctionSpace::Pk)
                        return as.trial.GetOrder() + as.test.GetOrder() - 2;
                     else
                        return as.trial.GetOrder() + as.test.GetOrder() + as.trial.GetDim() - 1;
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

         virtual void getElementMatrix(const Bilinear::Assembly::Common& as) const override
         {
            const int order = getIntegrationOrder(as);
            if (&as.trial == &as.test)
            {
               const mfem::IntegrationRule* ir =
                  as.trial.Space() == mfem::FunctionSpace::rQk ?
                     &mfem::RefinedIntRules.Get(as.trial.GetGeomType(), order) :
                     &mfem::IntRules.Get(as.trial.GetGeomType(), order);
               auto q = getIntegrand().getLHS().getLHS().build();
               switch (getIntegrand().getLHS().getLHS().getRangeType())
               {
                  case RangeType::Scalar:
                  {
                     mfem::VectorDiffusionIntegrator bfi(q.template get<RangeType::Scalar>());
                     bfi.SetIntRule(ir);
                     bfi.AssembleElementMatrix(as.trial, as.trans, as.mat);
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
                  [](const Linear::Assembly::Common& as)
                  {
                     return 2 * as.fe.GetOrder();
                  });

            const int order = integrand.getRHS()
                                       .getFiniteElementSpace()
                                       .getHandle().GetFE(0)->GetOrder();
            setIntegrationOrder(
                  [order](const Linear::Assembly::Device&)
                  {
                     return 2 * order;
                  });
         }

         constexpr
         Integral(const Integral& other)
            : Parent(other),
              m_devIntOrder(other.m_devIntOrder)
         {}

         constexpr
         Integral(Integral&& other)
            : Parent(std::move(other)),
              m_devIntOrder(std::move(other.m_devIntOrder))
         {}

         Integral& setIntegrationOrder(std::function<int(const Linear::Assembly::Common&)> order)
         {
            return static_cast<Integral&>(Parent::setIntegrationOrder(order));
         }

         int getIntegrationOrder(const Linear::Assembly::Common& as) const
         {
            return Parent::getIntegrationOrder(as);
         }

         Integral& setIntegrationOrder(std::function<int(const Linear::Assembly::Device&)> order)
         {
            m_devIntOrder = order;
            return *this;
         }

         int getIntegrationOrder(const Linear::Assembly::Device& as) const
         {
            assert(m_devIntOrder);
            return m_devIntOrder(as);
         }

         virtual const Integrand& getIntegrand() const override
         {
            return static_cast<const Integrand&>(Parent::getIntegrand());
         }

         virtual bool isSupported(Linear::Assembly::Type t) const override
         {
            switch (t)
            {
               case Linear::Assembly::Type::Common:
                  return true;
               case Linear::Assembly::Type::Device:
                  return true;
            }
            return false;
         }

         virtual void getElementVector(const Linear::Assembly::Device& as) const override
         {
            const FunctionBase& f = getIntegrand().getLHS();
            const mfem::IntegrationRule *ir =
               &mfem::IntRules.Get(as.fes.GetFE(0)->GetGeomType(), getIntegrationOrder(as));
            auto q = f.build();

            switch (f.getRangeType())
            {
               case RangeType::Scalar:
               {
                  mfem::DomainLFIntegrator lfi(q.get<RangeType::Scalar>());
                  lfi.SetIntRule(ir);
                  lfi.AssembleDevice(as.fes, as.markers, as.vec);
                  break;
               }
               case RangeType::Vector:
               {
                  mfem::VectorDomainLFIntegrator lfi(q.get<RangeType::Vector>());
                  lfi.SetIntRule(ir);
                  lfi.AssembleDevice(as.fes, as.markers, as.vec);
                  break;
               }
               case RangeType::Matrix:
               {
                  assert(false); // Unsupported
                  break;
               }
            }
         }

         virtual void getElementVector(const Linear::Assembly::Common& as) const override
         {
            const FunctionBase& f = getIntegrand().getLHS();

            const mfem::IntegrationRule *ir =
               &mfem::IntRules.Get(as.fe.GetGeomType(), getIntegrationOrder(as));
            auto q = f.build();

            switch (f.getRangeType())
            {
               case RangeType::Scalar:
               {
                  mfem::DomainLFIntegrator lfi(q.get<RangeType::Scalar>());
                  lfi.SetIntRule(ir);
                  lfi.AssembleRHSElementVect(as.fe, as.trans, as.vec);
                  break;
               }
               case RangeType::Vector:
               {
                  mfem::VectorDomainLFIntegrator lfi(q.get<RangeType::Vector>());
                  lfi.SetIntRule(ir);
                  lfi.AssembleRHSElementVect(as.fe, as.trans, as.vec);
                  break;
               }
               case RangeType::Matrix:
               {
                  assert(false); // Unsupported
                  break;
               }
            }
         }

         virtual Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         std::function<int(const Linear::Assembly::Device&)> m_devIntOrder;
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
         IntegratorRegion getIntegratorRegion() const override { return IntegratorRegion::Boundary; }
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
