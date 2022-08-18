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
    * @brief Integral of the dot product of a trial and a test operator
    *
    * Given two operators defined over trial and test spaces @f$ U_h @f$ and @f$ V_h @f$,
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
      : public BilinearFormDomainIntegrator
   {
      public:
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
         Integral(
               const ShapeFunctionBase<TrialSpace>& lhs,
               const ShapeFunctionBase<TestSpace>& rhs)
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
            : BilinearFormDomainIntegrator(prod.getLHS().getLeaf(), prod.getRHS().getLeaf()),
              m_prod(prod),
              m_intOrder(
                    [](const Bilinear::Assembly::Common& as)
                    { return as.trial.GetOrder() + as.test.GetOrder() + as.trans.OrderW(); })
         {}

         Integral(const Integral& other)
            : BilinearFormDomainIntegrator(other),
              m_prod(other.m_prod),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : BilinearFormDomainIntegrator(std::move(other)),
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

         /**
          * @brief Sets the function which calculates the integration order
          * @param[in] order Function which computes the order of integration
          * @returns Reference to self (for method chaining)
          */
         int getIntegrationOrder(const Bilinear::Assembly::Common& as) const
         {
            return m_intOrder(as);
         }

         const Integrand& getIntegrand() const
         {
            return m_prod;
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
   Integral(
         const Dot<ShapeFunctionBase<TrialSpace>,
         ShapeFunctionBase<TestSpace>>&)
      -> Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   Integral(
         const ShapeFunctionBase<TrialSpace>& lhs,
         const ShapeFunctionBase<TestSpace>& rhs)
      -> Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;


   /**
    * @brief Integral of a scalar valued test function operator
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
   class Integral<ShapeFunctionBase<TestSpace>> : public LinearFormDomainIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<TestSpace>;

         Integral(FunctionBase&& lhs, ShapeFunctionBase<TestSpace>&& rhs)
            : Integral(Dot(std::move(lhs), std::move(rhs)))
         {}

         Integral(FunctionBase&& lhs, const ShapeFunctionBase<TestSpace>& rhs)
            : Integral(Dot(std::move(lhs), rhs))
         {}

         Integral(const FunctionBase& lhs, ShapeFunctionBase<TestSpace>&& rhs)
            : Integral(Dot(lhs, std::move(rhs)))
         {}

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
            :  LinearFormDomainIntegrator(integrand.getLeaf()),
               m_integrand(integrand.copy()),
               m_intOrder(
                     [](const Linear::Assembly::Common& as)
                     { return as.fe.GetOrder() + as.trans.OrderW(); })
         {}

         Integral(const Integral& other)
            : LinearFormDomainIntegrator(other),
              m_integrand(other.m_integrand->copy()),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : LinearFormDomainIntegrator(std::move(other)),
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

         const Integrand& getIntegrand() const
         {
            assert(m_integrand);
            return *m_integrand;
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
   Integral(const ShapeFunctionBase<TestSpace>&) -> Integral<ShapeFunctionBase<TestSpace>>;
   Integral(const FunctionBase&, const ShapeFunctionBase<TestSpace>&) -> Integral<ShapeFunctionBase<TestSpace>>;

   /**
    * @brief Integral of a GridFunction
    */
   template <class FES>
   class Integral<GridFunction<FES>> : public FormLanguage::Base
   {
      public:
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

   template <>
   class BoundaryIntegral<ShapeFunctionBase<TestSpace>> : public LinearFormBoundaryIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<TestSpace>;

         BoundaryIntegral(
               const FunctionBase& lhs,
               const ShapeFunctionBase<TestSpace>& rhs)
            :  LinearFormBoundaryIntegrator(rhs.getLeaf()),
               m_integral(lhs, rhs)
         {}

         BoundaryIntegral(const Integrand& integrand)
            :  LinearFormBoundaryIntegrator(integrand.getLeaf()),
               m_integral(integrand)
         {}

         BoundaryIntegral(const BoundaryIntegral& other)
            : LinearFormBoundaryIntegrator(other),
              m_integral(other.m_integral)
         {}

         BoundaryIntegral(BoundaryIntegral&& other)
            : LinearFormBoundaryIntegrator(std::move(other)),
              m_integral(std::move(other.m_integral))
         {}

         void getElementVector(const Linear::Assembly::Common& as) const override
         {
            m_integral.getElementVector(as);
         }

         BoundaryIntegral* copy() const noexcept override
         {
            return new BoundaryIntegral(*this);
         }
      private:
         Integral<Integrand> m_integral;
   };
   BoundaryIntegral(const ShapeFunctionBase<TestSpace>&)
      -> BoundaryIntegral<ShapeFunctionBase<TestSpace>>;
   BoundaryIntegral(const FunctionBase& lhs, const ShapeFunctionBase<TestSpace>& rhs)
      -> BoundaryIntegral<ShapeFunctionBase<TestSpace>>;

   template <>
   class BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
      : public BilinearFormBoundaryIntegrator
   {
      public:
         using Integrand = Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;

         BoundaryIntegral(
               const ShapeFunctionBase<TrialSpace>& lhs,
               const ShapeFunctionBase<TestSpace>& rhs)
            :  BilinearFormBoundaryIntegrator(lhs.getLeaf(), rhs.getLeaf()),
               m_integral(lhs, rhs)
         {}

         BoundaryIntegral(const Integrand& integrand)
            :  BilinearFormBoundaryIntegrator(integrand.getLHS(), integrand.getRHS()),
               m_integral(integrand)
         {}

         BoundaryIntegral(const BoundaryIntegral& other)
            : BilinearFormBoundaryIntegrator(other),
              m_integral(other.m_integral)
         {}

         BoundaryIntegral(BoundaryIntegral&& other)
            : BilinearFormBoundaryIntegrator(std::move(other)),
              m_integral(std::move(other.m_integral))
         {}

         void getElementMatrix(const Bilinear::Assembly::Common& as) const override
         {
            return m_integral.getElementMatrix(as);
         }

         BoundaryIntegral* copy() const noexcept override
         {
            return new BoundaryIntegral(*this);
         }

      private:
         Integral<Integrand> m_integral;
   };
   BoundaryIntegral(const Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
   BoundaryIntegral(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;

   /* ||-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @internal
    *
    * Optimized integration of the expression:
    * @f[
    *    \int_\Omega \nabla u \cdot \nabla v \ dx
    * @f]
    * where $f$ is a function (scalar or matrix valued).
    */
   template <class FES>
   class Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>
      : public Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
   {
      public:
         using Parent = Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>;
         using Integrand = Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>;

         constexpr
         Integral(const Grad<ShapeFunction<FES, TrialSpace>>& gu, const Grad<ShapeFunction<FES, TestSpace>>& gv)
            : Integral(Dot(gu, gv))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>(integrand)
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
            : Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>(std::move(other))
         {}

         void getElementMatrix(const Bilinear::Assembly::Common& as) const override
         {
            assert(as.trial.GetDim() == as.test.GetDim());
            const int dim  = as.trial.GetDim();
            const int spaceDim = as.trans.GetSpaceDim();
            const bool square = (dim == spaceDim);
            const int order = getIntegrationOrder(as);
            const mfem::IntegrationRule* ir =
               as.trial.Space() == mfem::FunctionSpace::rQk ?
                  &mfem::RefinedIntRules.Get(as.trial.GetGeomType(), order) :
                  &mfem::IntRules.Get(as.trial.GetGeomType(), order);

            if (&as.trial == &as.test)
            {
               const int nd = as.trial.GetDof();

               as.mat.SetSize(nd);
               as.mat = 0.0;

               mfem::DenseMatrix dshape(nd, dim);
               mfem::DenseMatrix dshapedxt(nd, spaceDim);
               for (int i = 0; i < ir->GetNPoints(); i++)
               {
                  const mfem::IntegrationPoint &ip = ir->IntPoint(i);
                  as.trial.CalcDShape(ip, dshape);
                  as.trans.SetIntPoint(&ip);

                  const double tw = as.trans.Weight();
                  const double w = ip.weight / (square ? tw : tw * tw * tw);

                  mfem::Mult(dshape, as.trans.AdjugateJacobian(), dshapedxt);
                  mfem::AddMult_a_AAt(w, dshapedxt, as.mat);
               }
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
   Integral(const Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>> integrand)
      -> Integral<Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>>;

   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>>
    * ----------------------------------------------------------------------||
    */

   /* ||-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<ShapeFunctionBase<TestSpace>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @internal
    *
    * Optimized integration of the expression:
    * @f[
    *    \int_\Omega f v \ dx
    * @f]
    * where @f$ f @f$ is a function (scalar or matrix valued).
    */
   template <class FES>
   class Integral<Mult<FunctionBase, ShapeFunction<FES, TestSpace>>>
      : public Integral<ShapeFunctionBase<TestSpace>>
   {
      public:
         using Parent = Integral<ShapeFunctionBase<TestSpace>>;
         using Integrand = Mult<FunctionBase, ShapeFunction<FES, TestSpace>>;

         constexpr
         Integral(const Integrand& integrand)
            : Integral<ShapeFunctionBase<TestSpace>>(integrand)
         {
            if (integrand.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, integrand.getRangeType()).raise();
         }

         constexpr
         Integral(const Integral& other)
            : Integral<ShapeFunctionBase<TestSpace>>(other)
         {}

         constexpr
         Integral(Integral&& other)
            : Integral<ShapeFunctionBase<TestSpace>>(std::move(other))
         {}

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
   };
   template <class FES>
   Integral(const Mult<FunctionBase, ShapeFunction<FES, TestSpace>>&)
      -> Integral<Mult<FunctionBase, ShapeFunction<FES, TestSpace>>>;

   /**
    * @internal
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
         using Parent = Integral<ShapeFunctionBase<TestSpace>>;
         using Integrand = Dot<FunctionBase, ShapeFunction<FES, TestSpace>>;

         constexpr
         Integral(const FunctionBase& f, const ShapeFunction<FES, TestSpace>& u)
            : Integral(Dot(f, u))
         {}

         constexpr
         Integral(const Integrand& integrand)
            : Integral<ShapeFunctionBase<TestSpace>>(integrand)
         {
            setIntegrationOrder(
                  [](const Linear::Assembly::Common& as)
                  {
                     return 2 * as.fe.GetOrder();
                  });

            const int order = integrand.getShapeFunction()
                                       .getFiniteElementSpace()
                                       .getHandle().GetFE(0)->GetOrder();
            setIntegrationOrder(
                  [order](const Linear::Assembly::Device&)
                  {
                     return 2 * order;
                  }
                  );
         }

         constexpr
         Integral(const Integral& other)
            : Integral<ShapeFunctionBase<TestSpace>>(other),
              m_devIntOrder(other.m_devIntOrder)
         {}

         constexpr
         Integral(Integral&& other)
            : Integral<ShapeFunctionBase<TestSpace>>(std::move(other)),
              m_devIntOrder(std::move(other.m_devIntOrder))
         {}

         const Integrand& getIntegrand() const
         {
            return static_cast<const Integrand&>(Integral<ShapeFunctionBase<TestSpace>>::getIntegrand());
         }

         bool isSupported(Linear::Assembly::Type t) const override
         {
            switch (t)
            {
               case Linear::Assembly::Type::Common:
                  return true;
               case Linear::Assembly::Type::Device:
                  return true;
               default:
                  return false;
            }
            return false;
         }

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

         void getElementVector(const Linear::Assembly::Device& as) const override
         {
            const FunctionBase& f = getIntegrand().getFunction();
            const mfem::IntegrationRule *ir =
               &mfem::IntRules.Get(as.fes.GetFE(0)->GetGeomType(), getIntegrationOrder(as));
            auto q = f.build();

            switch (f.getRangeType())
            {
               case RangeType::Scalar:
               {
                  mfem::DomainLFIntegrator lfi(std::get<Internal::ScalarProxyFunction>(q));
                  lfi.SetIntRule(ir);
                  lfi.AssembleDevice(as.fes, as.markers, as.vec);
                  break;
               }
               case RangeType::Vector:
               {
                  mfem::VectorDomainLFIntegrator lfi(std::get<Internal::VectorProxyFunction>(q));
                  lfi.SetIntRule(ir);
                  lfi.AssembleDevice(as.fes, as.markers, as.vec);
               }
               case RangeType::Matrix:
               {
                  assert(false); // Unsupported
                  break;
               }
            }
         }

         void getElementVector(const Linear::Assembly::Common& as) const override
         {
            const FunctionBase& f = getIntegrand().getFunction();

            const mfem::IntegrationRule *ir =
               &mfem::IntRules.Get(as.fe.GetGeomType(), getIntegrationOrder(as));
            auto q = f.build();

            switch (f.getRangeType())
            {
               case RangeType::Scalar:
               {
                  mfem::DomainLFIntegrator lfi(std::get<Internal::ScalarProxyFunction>(q));
                  lfi.SetIntRule(ir);
                  lfi.AssembleRHSElementVect(as.fe, as.trans, as.vec);
                  break;
               }
               case RangeType::Vector:
               {
                  mfem::VectorDomainLFIntegrator lfi(std::get<Internal::VectorProxyFunction>(q));
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

         Integral* copy() const noexcept override
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


   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Integral<ShapeFunctionBase<TestSpace>>
    * ----------------------------------------------------------------------||
    */

}

#endif
