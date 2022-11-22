/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_AVERAGE_H
#define RODIN_VARIATIONAL_AVERAGE_H

#include "ForwardDecls.h"
#include "BasisOperator.h"
#include "ShapeFunction.h"

namespace Rodin::Variational
{
   template <ShapeFunctionSpaceType Space>
   class Average<ShapeFunctionBase<Space>> : public ShapeFunctionBase<Space>
   {
      public:
         using Parent = ShapeFunctionBase<Space>;

         constexpr
         Average(ShapeFunctionBase<Space>& u)
            : m_u(u)
         {}

         constexpr
         Average(const Average& other)
            :  Parent(other),
               m_u(other.m_u)
         {}

         constexpr
         Average(Average&& other)
            :  Parent(std::move(other)),
               m_u(std::move(other.m_u))
         {}

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return m_u.get().getLeaf();
         }

         int getRows() const override
         {
            return m_u.get().getRows();
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_u.get().getDOFs(fe, trans);
         }

         int getColumns() const override
         {
            return m_u.get().getColumns();
         }

         // void getOperator(
         //       DenseBasisOperator& op,
         //       ShapeComputator& compute,
         //       const Bilinear::Assembly::Native::Interface& data) override
         // {
         //    getOperator(op, fe, trans.GetElement1Transformation(), ip1, compute);
         //    DenseBasisOperator tmp;
         //    getOperator(tmp, fe, trans.GetElement2Transformation(), ip2, compute);
         //    op += tmp;
         //    op *= 0.5;
         // }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return m_u.get().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return m_u.get().getFiniteElementSpace();
         }

         Average* copy() const noexcept override
         {
            return new Average(*this);
         }

      private:
         std::reference_wrapper<ShapeFunctionBase<Space>> m_u;
   };
}

#endif
