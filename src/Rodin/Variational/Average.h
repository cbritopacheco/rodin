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

         int getDOFs(const Geometry::Simplex& element) const override
         {
            assert(false);
            // assert(dynamic_cast<const Geometry::Interface*>(&element));
            // const auto& incident =
            //    static_cast<const Geometry::Interface&>(element).getElements();
            // assert(incident.size() == 2);
            // auto first = incident.begin();
            // auto second = std::next(first);
            // return m_u.get().getDOFs(*first) + m_u.get().getDOFs(*second);
         }

         int getColumns() const override
         {
            return m_u.get().getColumns();
         }

         void getOperator(
               DenseBasisOperator& op,
               ShapeComputator& compute,
               const Geometry::Simplex& element) override
         {
            assert(false);
            // switch (element.getRegion())
            // {
            //    case Geometry::Region::Interface:
            //    {
            //       assert(dynamic_cast<const Geometry::Interface*>(&element));
            //       const auto& incident = static_cast<const Geometry::Interface&>(element).getElements();
            //       assert(incident.size() == 2);
            //       auto first = incident.begin();
            //       auto second = std::next(first);

            //       first->getTransformation().SetIntPoint(&element.getTransformation().GetIntPoint());
            //       second->getTransformation().SetIntPoint(&element.getTransformation().GetIntPoint());

            //       const int ndofs1 = m_u.get().getDOFs(*first);
            //       const int ndofs2 = m_u.get().getDOFs(*second);

            //       DenseBasisOperator op1;
            //       m_u.get().getOperator(op1, compute, first);
            //       op1 *= 0.5;

            //       DenseBasisOperator op2;
            //       m_u.get().getOperator(op2, compute, second);
            //       op1 *= 0.5;

            //       DenseBasisOperator res;
            //       res.setSize(getRows(), getColumns(), ndofs1 + ndofs2);
            //       res = 0.0;

            //       for (int i = 0; i < ndofs1; i++)
            //          res(i) = std::move(op1(i));

            //       for (int i = 0; i < ndofs2; i++)
            //          res(i + ndofs1) = std::move(op2(i));

            //       break;
            //    }
            //    default:
            //    {
            //       assert(false);
            //       break;
            //    }
            // }
         }

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
