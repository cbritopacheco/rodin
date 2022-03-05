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
#include "Jacobian.h"
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
   class Grad<GridFunction<H1>> : public VectorCoefficientBase
   {
      public:
         /**
          * @brief Constructs the gradient of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Grad(const GridFunction<H1>& u)
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
            return new Grad(*this);
         }

      private:
         const GridFunction<H1>& m_u;
         mutable mfem::GradientGridFunctionCoefficient m_mfemVectorCoefficient;
   };
   Grad(const GridFunction<H1>&) -> Grad<GridFunction<H1>>;

   template <ShapeFunctionSpaceType Space>
   class Grad<ShapeFunction<H1, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Grad(ShapeFunction<H1, Space>& u)
            : m_u(u)
         {
            assert(m_u.getFiniteElementSpace().getVectorDimension() == 1);
         }

         int getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation& trans) const override
         {
            return trans.GetSpaceDim();
         }

         int getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const
         {
            return m_u.getDOFs(fe, trans);
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            int dofs = fe.GetDof();
            int sdim = trans.GetSpaceDim();
            mfem::DenseMatrix dshape;
            dshape.SetSize(dofs, sdim);
            fe.CalcPhysDShape(trans, dshape);
            return std::unique_ptr<Rank3Operator>(
                  new Internal::JacobianShapeR3O(std::move(dshape), sdim, 1));
         }

         H1& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const H1& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Grad* copy() const noexcept override
         {
            return new Grad(*this);
         }
      private:
         ShapeFunction<H1, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Grad(ShapeFunction<H1, Space>&) -> Grad<ShapeFunction<H1, Space>>;
}

#endif
