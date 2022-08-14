/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "H1.h"
#include "GridFunction.h"
#include "BasisOperator.h"
#include "ShapeFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the Jacobian matrix @f$ \mathbf{J}_u @f$ of the
    * function @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^s \rightarrow \mathbb{R}^d @f$, the Jacobian matrix
    * @f$ \mathbf{J}_u(x) @f$ at any point @f$ x = (x_1, \ldots, x_s) @f$ is
    * defined by the @f$ s \times d @f$ matrix
    * @f[
    * \mathbf{J}_u = \begin{bmatrix}
    * \dfrac{\partial u_1}{\partial x_1} & \ldots & \dfrac{\partial u_d}{\partial x_1}\\
    * \vdots & \ddots & \vdots\\
    * \dfrac{\partial u_1}{\partial x_s} & \ldots & \dfrac{\partial u_d}{\partial x_s}
    * \end{bmatrix} .
    * @f]
    * This class aids in the calculation of the Jacobian of a GridFunction<H1>.
    */
   template <class Trait>
   class Jacobian<GridFunction<H1, Trait>> : public MatrixFunctionBase
   {
      public:
         /**
          * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Jacobian(GridFunction<H1, Trait>& u)
            :  m_u(u)
         {}

         Jacobian(const Jacobian& other)
            :  MatrixFunctionBase(other),
               m_u(other.m_u)
         {}

         Jacobian(Jacobian&& other)
            : MatrixFunctionBase(std::move(other)),
              m_u(other.m_u)
         {}

         int getRows() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getDimension();
         }

         int getColumns() const override
         {
            return m_u.getFiniteElementSpace().getVectorDimension();
         }

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            m_u.getHandle().GetVectorGradient(
                  FunctionBase::getTraceElementTrans(trans, ip), value);
            value.Transpose();
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }

      private:
         GridFunction<H1, Trait>& m_u;
   };
   template <class Trait>
   Jacobian(GridFunction<H1, Trait>&) -> Jacobian<GridFunction<H1, Trait>>;

   template <ShapeFunctionSpaceType Space>
   class Jacobian<ShapeFunction<H1, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Jacobian(ShapeFunction<H1, Space>& u)
            : m_u(u)
         {}

         Jacobian(const Jacobian& other)
            : ShapeFunctionBase<Space>(other),
              m_u(other.m_u)
         {}

         Jacobian(Jacobian&& other)
            : ShapeFunctionBase<Space>(std::move(other)),
              m_u(other.m_u)
         {}

         FiniteElementSpace<H1>& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const FiniteElementSpace<H1>& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         const ShapeFunction<H1, Space>& getLeaf() const override
         {
            return m_u.getLeaf();
         }

         int getRows() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getSpaceDimension();
         }

         int getColumns() const override
         {
            return m_u.getFiniteElementSpace().getVectorDimension();
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_u.getDOFs(fe, trans);
         }

         std::unique_ptr<BasisOperator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            int dofs = fe.GetDof();
            int sdim = trans.GetSpaceDim();
            int vdim = m_u.getFiniteElementSpace().getVectorDimension();
            mfem::DenseMatrix dshape;
            dshape.SetSize(dofs, sdim);
            fe.CalcPhysDShape(trans, dshape);
            return std::unique_ptr<BasisOperator>(new JSSFBO(std::move(dshape), sdim, vdim));
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }
      private:
         ShapeFunction<H1, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Jacobian(ShapeFunction<H1, Space>&) -> Jacobian<ShapeFunction<H1, Space>>;
}

#endif
