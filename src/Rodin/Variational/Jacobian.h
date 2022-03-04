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
#include "Rank3Operator.h"
#include "ShapeFunction.h"
#include "VectorCoefficient.h"
#include "MatrixCoefficient.h"

namespace Rodin::Variational::Internal
{
   class JacobianShapeR3O : public Rank3Operator
   {
      public:
         JacobianShapeR3O(mfem::DenseMatrix dshape, int sdim, int vdim)
            : m_dshape(dshape),
              m_sdim(sdim),
              m_vdim(vdim)
         {}

         int GetRows() const override;

         int GetColumns() const override;

         int GetDOFs() const override;

         JacobianShapeR3O& operator=(double s) override;

         JacobianShapeR3O& operator*=(double s) override;

         double operator()(int row, int col, int dof) const override;

         std::unique_ptr<Rank3Operator> Trace() const override;

      private:
         mfem::DenseMatrix m_dshape;
         int m_sdim, m_vdim;
   };
}

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
    * \end{bmatrix}
    * @f]
    * 
    */
   template <>
   class Jacobian<GridFunction<H1>> : public MatrixCoefficientBase
   {
      public:
         /**
          * @brief Constructs the Jacobian matrix of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Jacobian(GridFunction<H1>& u)
            :  m_u(u)
         {}

         Jacobian(const Jacobian& other)
            : m_u(other.m_u)
         {}

         int getRows() const override
         {
            return m_u.getFiniteElementSpace().getVectorDimension();
         }

         int getColumns() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getDimension();
         }

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint&) override
         {
            m_u.getHandle().GetVectorGradient(trans, value);
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }

      private:
         GridFunction<H1>& m_u;
   };
   Jacobian(GridFunction<H1>&) -> Jacobian<GridFunction<H1>>;

   template <ShapeFunctionSpaceType Space>
   class Jacobian<ShapeFunction<H1, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Jacobian(const ShapeFunction<H1, Space>& u)
            : m_u(u)
         {}

         H1& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const H1& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
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
            return m_u.getFiniteElementSpace().getVectorDimension();
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
            int vdim = m_u.getFiniteElementSpace().getVectorDimension();
            mfem::DenseMatrix dshape;
            dshape.SetSize(dofs, sdim);
            return std::unique_ptr<Rank3Operator>(
                  new Internal::JacobianShapeR3O(std::move(dshape), sdim, vdim));
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }
      private:
         const ShapeFunction<H1, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Jacobian(ShapeFunction<H1, Space>&) -> Jacobian<ShapeFunction<H1, Space>>;
}

#endif
