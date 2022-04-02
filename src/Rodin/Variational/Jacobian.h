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
#include "VectorFunction.h"
#include "MatrixFunction.h"

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

         JacobianShapeR3O(const JacobianShapeR3O& other)
            : m_dshape(other.m_dshape),
              m_sdim(other.m_sdim),
              m_vdim(other.m_vdim)
         {}

         JacobianShapeR3O(JacobianShapeR3O&& other)
            : m_dshape(std::move(other.m_dshape)),
              m_sdim(other.m_sdim),
              m_vdim(other.m_vdim)
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
    * \end{bmatrix} .
    * @f]
    * This class aids in the calculation of the Jacobian of a GridFunction<H1>.
    */
   template <>
   class Jacobian<GridFunction<H1>> : public MatrixFunctionBase
   {
      public:
         /**
          * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Jacobian(GridFunction<H1>& u)
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
            const auto& traceDomain = getTraceDomain();
            switch (trans.ElementType)
            {
               case mfem::ElementTransformation::BDR_ELEMENT:
               {
                  int fn = trans.mesh->GetBdrFace(trans.ElementNo);
                  if (trans.mesh->FaceIsInterior(fn))
                  {
                     if (traceDomain.empty())
                     {
                        Alert::Exception()
                           << "Integration over an interior boundary element "
                           << "requires a trace domain."
                           << Alert::Raise;
                     }

                     // Extend the values on the trace domain up to the
                     // interior boundary.
                     mfem::FaceElementTransformations* ft =
                        trans.mesh->GetFaceElementTransformations(fn);
                     ft->SetAllIntPoints(&ip);
                     if (traceDomain.count(ft->GetElement1Transformation().Attribute))
                        m_u.getHandle().GetVectorGradient(ft->GetElement1Transformation(), value);
                     else if (traceDomain.count(ft->GetElement2Transformation().Attribute))
                        m_u.getHandle().GetVectorGradient(ft->GetElement2Transformation(), value);
                     else
                     {
                        // The boundary over which we are evaluating must be
                        // the interface between the trace domain and some
                        // other domain, i.e. it is not the boundary that was
                        // specified!
                        auto ex = Alert::Exception();
                        ex << "Boundary element " << trans.ElementNo
                           << " with attribute " << trans.Attribute
                           << " is not a boundary of the trace domain.";
                        ex.raise();
                     }
                  }
                  else
                  {
                     m_u.getHandle().GetVectorGradient(trans, value);
                  }
                  break;
               }
               default:
                  m_u.getHandle().GetVectorGradient(trans, value);
            }
            value.Transpose();
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
         Jacobian(ShapeFunction<H1, Space>& u)
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

         ShapeFunction<H1, Space>& getRoot() override
         {
            return m_u;
         }

         const ShapeFunction<H1, Space>& getRoot() const override
         {
            return m_u;
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
               const mfem::ElementTransformation& trans) const override
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
            fe.CalcPhysDShape(trans, dshape);
            return std::unique_ptr<Rank3Operator>(
                  new Internal::JacobianShapeR3O(std::move(dshape), sdim, vdim));
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
