#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"

#include "H1.h"
#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the divergence @f \text{div} u @f of a function @f$ u \in H^1 @f$.
    *
    * This class represents the sum of the first derivatives of a
    * vector valued H1 ShapeFunction.
    * @f[
    *    \text{div} u
    * @f]
    */
   template <ShapeFunctionSpaceType Space>
   class Div<ShapeFunction<H1, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         /**
          * @brief Constructs Div object
          * @param[in] u ShapeFunction to be differentiated
          */
         Div(ShapeFunction<H1, Space>& u)
            : m_u(u)
         {
            assert(m_u.getFiniteElementSpace().getVectorDimension() > 1);
         }

         ShapeFunction<H1, Space>& getLeaf() override
         {
            return m_u;
         }

         const ShapeFunction<H1, Space>& getLeaf() const override
         {
            return m_u;
         }

         int getRows(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
         }

         int getColumns(
               const mfem::FiniteElement&,
               const mfem::ElementTransformation&) const override
         {
            return 1;
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
            return Internal::JacobianShapeR3O(std::move(dshape), sdim, vdim).Trace();
         }

         FiniteElementSpace<H1>& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const FiniteElementSpace<H1>& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Div* copy() const noexcept override
         {
            return new Div(*this);
         }
      private:
         ShapeFunction<H1, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Div(ShapeFunction<H1, Space>&) -> Div<ShapeFunction<H1, Space>>;
}

#endif
