/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "ForwardDecls.h"

#include "TrialFunction.h"
#include "TestFunction.h"
#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the gradient @f$ \nabla u @f$ of the scalar function
    * @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^n \rightarrow \mathbb{R} @f$, the gradient
    * @f$ \nabla u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ at the point
    * @f$ x = (x_1, \ldots, x_n) @f$ is defined by
    * @f[
    *    \nabla u (x) =
    *    \left[
    *       \dfrac{\partial u}{\partial x_1}(x), \ldots,
    *       \dfrac{\partial u}{\partial x_n}(x)
    *    \right]^T
    * @f]
    */
   template <>
   class Gradient<GridFunction<H1>> : public VectorCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Gradient of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Gradient(const GridFunction<H1>& u);

         size_t getDimension() const override;

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            m_mfemVectorCoefficient.Eval(value, trans, ip);
         }

         VectorCoefficientBase* copy() const noexcept override
         {
            return new Gradient(*this);
         }

      private:
         const GridFunction<H1>& m_u;
         mfem::GradientGridFunctionCoefficient m_mfemVectorCoefficient;
   };
   Gradient(const GridFunction<H1>&) -> Gradient<GridFunction<H1>>;

   template <>
   class Gradient<TrialFunction<H1>> : public TrialFunctionBase
   {
      public:
         Gradient(const TrialFunction<H1>& u);

         ValueType getValueType() const override;

         void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               VectorShape& values
               ) const override;

         const FiniteElementSpaceBase& getFiniteElementSpace() const override;

         Gradient* copy() const noexcept override
         {
            return new Gradient(*this);
         }
      private:
         const TrialFunction<H1>& m_u;
   };
   Gradient(const TrialFunction<H1>&) -> Gradient<TrialFunction<H1>>;

   template <>
   class Gradient<TestFunction<H1>> : public TestFunctionBase
   {
      public:
         Gradient(const TestFunction<H1>& v);

         ValueType getValueType() const override;

         void getValue(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               VectorShape& values
               ) const override;

         const FiniteElementSpaceBase& getFiniteElementSpace() const override;

         Gradient* copy() const noexcept override
         {
            return new Gradient(*this);
         }
      private:
         const TestFunction<H1>& m_v;
   };
   Gradient(const TestFunction<H1>&) -> Gradient<TestFunction<H1>>;
}

#endif
