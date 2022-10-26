#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"

#include "H1.h"
#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "ScalarFunction.h"
#include "Trace.h"

namespace Rodin::Variational
{
   /**
    * @defgroup DivSpecializations Div Template Specializations
    * @brief Template specializations of the Div class.
    * @see Div
    */

   /**
    * @ingroup DivSpecializations
    */
   template <ShapeFunctionSpaceType Space, class Trait>
   class Div<ShapeFunction<H1<Trait>, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         /**
          * @brief Constructs Div object
          * @param[in] u ShapeFunction to be differentiated
          */
         Div(ShapeFunction<H1<Trait>, Space>& u)
            : m_u(u)
         {
            if (u.getRangeType() != RangeType::Scalar &&
                  u.getRangeType() != RangeType::Vector)
            {
               UnexpectedRangeTypeException(
                     {RangeType::Scalar, RangeType::Vector}, u.getRangeType());
            }
         }

         Div(const Div& other)
            : ShapeFunctionBase<Space>(other),
              m_u(other.m_u)
         {}

         Div(Div&& other)
            : ShapeFunctionBase<Space>(std::move(other)),
              m_u(other.m_u)
         {}

         const ShapeFunction<H1<Trait>, Space>& getLeaf() const override
         {
            return m_u.getLeaf();
         }

         int getRows() const override
         {
            return 1;
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_u.getDOFs(fe, trans);
         }

         void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip,
               ShapeComputator& comp) const override
         {
            const auto& dshape = comp.getPhysicalDShape(fe, trans, ip);
            const int opDofs = getDOFs(fe, trans);
            op.setSize(1, 1, opDofs);
            for (int i = 0; i < opDofs; i++)
               op(0, 0, i) = dshape.GetData()[i];
         }

         H1<Trait>& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const H1<Trait>& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Div* copy() const noexcept override
         {
            return new Div(*this);
         }
      private:
         ShapeFunction<H1<Trait>, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space, class Trait>
   Div(ShapeFunction<H1<Trait>, Space>&) -> Div<ShapeFunction<H1<Trait>, Space>>;
}

#endif
