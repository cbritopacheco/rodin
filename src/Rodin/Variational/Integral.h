/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRAL_H
#define RODIN_VARIATIONAL_INTEGRAL_H

#include <boost/uuid/uuid_io.hpp>
#include <set>
#include <utility>
#include <mfem.hpp>

#include "Dot.h"
#include "LinearForm.h"
#include "ForwardDecls.h"
#include "FormLanguage.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
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
   class Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
      : public BilinearFormDomainIntegrator
   {
      public:
         using Integrand = Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;

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
         Integral(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
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
                    [](const mfem::FiniteElement& trial,
                       const mfem::FiniteElement& test,
                       mfem::ElementTransformation& trans)
                    { return trial.GetOrder() + test.GetOrder() + trans.OrderW(); })
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
            std::function<
               int(const mfem::FiniteElement&, const mfem::FiniteElement&, mfem::ElementTransformation&)> order)
         {
            m_intOrder = order;
            return *this;
         }

         /**
          * @brief Sets the function which calculates the integration order
          * @param[in] order Function which computes the order of integration
          * @returns Reference to self (for method chaining)
          */
         int getIntegrationOrder(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans) const
         {
            return m_intOrder(trial, test, trans);
         }

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) const override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         Integrand m_prod;
         std::function<int(
            const mfem::FiniteElement&, const mfem::FiniteElement&, mfem::ElementTransformation&)> m_intOrder;
   };
   Integral(const Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>&)
      -> Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;
   Integral(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
      -> Integral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;


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
   class Integral<ShapeFunctionBase<Test>>
      : public LinearFormDomainIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         /**
          * @brief Integral of the multiplication between a scalar valued function and operator
          *
          * Given
          * @f[
          * \lambda : \Omega \rightarrow \mathbb{R}, \quad
          * A : V_h \rightarrow \mathbb{R}
          * @f]
          * constructs an instance representing the following integral
          * @f[
          *    \int \lambda A(v) \ dx \ .
          * @f]
          */
         Integral(const ScalarFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
            : Integral(Dot(lhs, rhs))
         {}

         /**
          * @brief Integral of the dot product between a vector valued function and operator
          *
          * Given
          * @f[
          * \vec{\lambda} : \Omega \rightarrow \mathbb{R}^d, \quad
          * \vec{A} : V_h \rightarrow \mathbb{R}^d
          * @f]
          * constructs an instance representing the following integral
          * @f[
          *    \int \vec{\lambda} \cdot \vec{A}(v) \ dx \ .
          * @f]
          */
         Integral(const VectorFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
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
               m_test(integrand.copy()),
               m_intOrder(
                    [](const mfem::FiniteElement& fe,
                       mfem::ElementTransformation& trans)
                    { return fe.GetOrder() + trans.OrderW(); })
         {}

         Integral(const Integral& other)
            : LinearFormDomainIntegrator(other),
              m_test(other.m_test->copy()),
              m_intOrder(other.m_intOrder)
         {}

         Integral(Integral&& other)
            : LinearFormDomainIntegrator(std::move(other)),
              m_test(std::move(other.m_test))
         {}

         /**
          * @brief Sets the function which calculates the integration order
          * @param[in] order Function which computes the order of integration
          * @returns Reference to self (for method chaining)
          */
         Integral& setIntegrationOrder(
            std::function<
               int(const mfem::FiniteElement&, mfem::ElementTransformation&)> order)
         {
            m_intOrder = order;
            return *this;
         }

         /**
          * @brief Sets the function which calculates the integration order
          * @param[in] order Function which computes the order of integration
          * @returns Reference to self (for method chaining)
          */
         int getIntegrationOrder(
               const mfem::FiniteElement& fe, mfem::ElementTransformation& trans) const
         {
            return m_intOrder(fe, trans);
         }

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) const override;

         Integral* copy() const noexcept override
         {
            return new Integral(*this);
         }
      private:
         std::unique_ptr<Integrand> m_test;
         std::function<int(
            const mfem::FiniteElement&, mfem::ElementTransformation&)> m_intOrder;
   };
   Integral(const ShapeFunctionBase<Test>&)
      -> Integral<ShapeFunctionBase<Test>>;
   Integral(const ScalarFunctionBase&, const ShapeFunctionBase<Test>&)
      -> Integral<ShapeFunctionBase<Test>>;
   Integral(const VectorFunctionBase&, const ShapeFunctionBase<Test>&)
      -> Integral<ShapeFunctionBase<Test>>;

   /**
    * @brief Integral of a GridFunction
    */
   template <class FEC, class Trait>
   class Integral<GridFunction<FEC, Trait>> : public FormLanguage::Base
   {
      public:
         /**
          * @brief Constructs the integral object
          */
         Integral(GridFunction<FEC, Trait>& u)
            : m_u(u),
              m_v(u.getFiniteElementSpace()),
              m_one(u.getFiniteElementSpace()),
              m_lf(m_v),
              m_assembled(false)
         {
            assert(u.getFiniteElementSpace().getVectorDimension() == 1);
            m_one = ScalarFunction(1.0);
            m_lf.from(Integral<ShapeFunctionBase<Test>>(ScalarFunction(u) * m_v));
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
         TestFunction<FEC, Trait>    m_v;
         GridFunction<FEC, Trait>&   m_u;
         GridFunction<FEC, Trait>    m_one;
         LinearForm<FEC, Trait>      m_lf;
         bool m_assembled;
   };
   template <class FEC, class Trait>
   Integral(GridFunction<FEC, Trait>&) -> Integral<GridFunction<FEC, Trait>>;

   template <>
   class BoundaryIntegral<ShapeFunctionBase<Test>> : public LinearFormBoundaryIntegrator
   {
      public:
         using Integrand = ShapeFunctionBase<Test>;

         BoundaryIntegral(const ScalarFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
            :  LinearFormBoundaryIntegrator(rhs.getLeaf()),
               m_integral(lhs, rhs)
         {}

         BoundaryIntegral(const VectorFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
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

         void getElementVector(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans, mfem::Vector& vec) const override
         {
            m_integral.getElementVector(fe, trans, vec);
         }

         BoundaryIntegral* copy() const noexcept override
         {
            return new BoundaryIntegral(*this);
         }
      private:
         Integral<Integrand> m_integral;
   };
   BoundaryIntegral(const ShapeFunctionBase<Test>&)
      -> BoundaryIntegral<ShapeFunctionBase<Test>>;
   BoundaryIntegral(const ScalarFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
      -> BoundaryIntegral<ShapeFunctionBase<Test>>;
   BoundaryIntegral(const VectorFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
      -> BoundaryIntegral<ShapeFunctionBase<Test>>;
   BoundaryIntegral(const MatrixFunctionBase& lhs, const ShapeFunctionBase<Test>& rhs)
      -> BoundaryIntegral<ShapeFunctionBase<Test>>;

   template <>
   class BoundaryIntegral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>
      : public BilinearFormBoundaryIntegrator
   {
      public:
         using Integrand = Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>;

         BoundaryIntegral(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs)
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

         void getElementMatrix(
               const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
               mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) const override
         {
            return m_integral.getElementMatrix(trial, test, trans, mat);
         }

         BoundaryIntegral* copy() const noexcept override
         {
            return new BoundaryIntegral(*this);
         }

      private:
         Integral<Integrand> m_integral;
   };
   BoundaryIntegral(const Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;
   BoundaryIntegral(const ShapeFunctionBase<Trial>&, const ShapeFunctionBase<Test>&)
      -> BoundaryIntegral<Dot<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>>;

}

#endif
