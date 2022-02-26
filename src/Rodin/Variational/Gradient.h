/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRADIENT_H
#define RODIN_VARIATIONAL_GRADIENT_H

#include "ForwardDecls.h"

#include "H1.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
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
         Gradient(const GridFunction<H1>& u)
            : m_u(u),
              m_mfemVectorCoefficient(&m_u.getHandle())
         {}

         size_t getDimension() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getDimension();
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            m_mfemVectorCoefficient.Eval(value, trans, ip);
         }

         VectorCoefficientBase* copy() const noexcept override
         {
            return new Gradient(*this);
         }

      private:
         const GridFunction<H1>& m_u;
         mutable mfem::GradientGridFunctionCoefficient m_mfemVectorCoefficient;
   };
   Gradient(const GridFunction<H1>&) -> Gradient<GridFunction<H1>>;

   template <ShapeFunctionSpaceType Space>
   class Gradient<ShapeFunction<H1, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Gradient(ShapeFunction<H1, Space>& u)
            : m_u(u)
         {
            assert(m_u.getFiniteElementSpace().getVectorDimension() == 1);
         }

         size_t getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation& trans) const override
         {
            return trans.GetSpaceDim();
         }

         size_t getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const
         {
            return fe.GetDof();
         }

         size_t getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         Internal::Rank3Operator getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            size_t dofs = fe.GetDof();
            size_t sdim = trans.GetSpaceDim();
            Internal::Rank3Operator result(getRows(fe, trans), getDOFs(fe, trans), getColumns(fe, trans));
            mfem::DenseMatrix tmp(dofs, sdim);
            fe.CalcPhysDShape(trans, tmp);
            for (size_t i = 0; i < sdim; i++)
            {
               for (size_t j = 0; j < dofs; j++)
               {
                  result(i, j, 0) = tmp(j, i);
               }
            }
            return result;
         }

         H1& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const H1& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Gradient* copy() const noexcept override
         {
            return new Gradient(*this);
         }
      private:
         ShapeFunction<H1, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Gradient(ShapeFunction<H1, Space>&) -> Gradient<ShapeFunction<H1, Space>>;
}

#endif
