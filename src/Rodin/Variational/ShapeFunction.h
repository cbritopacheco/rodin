#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "BasisOperator.h"
#include "FiniteElementSpace.h"

#include "RangeShape.h"

namespace Rodin::Variational
{
   template <ShapeFunctionSpaceType Space>
   struct DualSpaceType;

   template <>
   struct DualSpaceType<TrialSpace>
   {
      static constexpr ShapeFunctionSpaceType Value = ShapeFunctionSpaceType::Test;
   };

   template <>
   struct DualSpaceType<TestSpace>
   {
      static constexpr ShapeFunctionSpaceType Value = ShapeFunctionSpaceType::Trial;
   };

   class ShapeComputator
   {
      public:
         using Key =
            std::tuple<
               const mfem::FiniteElement*,
               mfem::ElementTransformation*>;

         template <class Data>
         struct Value
         {
            const mfem::IntegrationPoint* pip;
            Data data;
         };

         ShapeComputator() = default;

         const mfem::Vector& getShape(
               const mfem::FiniteElement& el,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip)
         {
            auto it = m_shapeLookup.find({&el, &trans});
            if (it != m_shapeLookup.end())
            {
               if (it->second.pip != &ip)
               {
                  it->second.pip = &ip;
                  it->second.data.SetSize(el.GetDof());
                  el.CalcShape(ip, it->second.data);
               }
               return it->second.data;
            }
            else
            {
               auto itt =  m_shapeLookup.insert(
                     it, {{&el, &trans}, {&ip, mfem::Vector(el.GetDof())}});
               el.CalcShape(ip, itt->second.data);
               return itt->second.data;
            }
         }

         const mfem::Vector& getPhysicalShape(
               const mfem::FiniteElement& el,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip)
         {
            auto it = m_physShapeLookup.find({&el, &trans});
            if (it != m_physShapeLookup.end())
            {
               if (it->second.pip != &ip)
               {
                  it->second.pip = &ip;
                  it->second.data.SetSize(el.GetDof());
                  el.CalcPhysShape(trans, it->second.data);
               }
               return it->second.data;
            }
            else
            {
               auto itt =  m_physShapeLookup.insert(
                     it, {{&el, &trans}, {&ip, mfem::Vector(el.GetDof())}});
               el.CalcPhysShape(trans, itt->second.data);
               return itt->second.data;
            }
         }

         const mfem::DenseMatrix& getDShape(
               const mfem::FiniteElement& el,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip)
         {
            auto it = m_dShapeLookup.find({&el, &trans});
            if (it != m_dShapeLookup.end())
            {
               it->second.pip = &ip;
               it->second.data.SetSize(el.GetDof(), el.GetDim());
               el.CalcDShape(ip, it->second.data);
               return it->second.data;
            }
            else
            {
               auto itt =  m_dShapeLookup.insert(
                     it, {{&el, &trans}, {&ip, mfem::DenseMatrix(el.GetDof(), el.GetDim())}});
               el.CalcDShape(ip, itt->second.data);
               return itt->second.data;
            }
         }

         const mfem::DenseMatrix& getPhysicalDShape(
               const mfem::FiniteElement& el,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip)
         {
            auto it = m_physDShapeLookup.find({&el, &trans});
            if (it != m_physDShapeLookup.end())
            {
               it->second.pip = &ip;
               it->second.data.SetSize(el.GetDof(), trans.GetSpaceDim());
               el.CalcPhysDShape(trans, it->second.data);
               return it->second.data;
            }
            else
            {
               auto itt =  m_physDShapeLookup.insert(
                     it, {{&el, &trans}, {&ip, mfem::DenseMatrix(el.GetDof(), trans.GetSpaceDim())}});
               el.CalcPhysDShape(trans, itt->second.data);
               return itt->second.data;
            }
         }

      private:
         std::map<Key, Value<mfem::Vector>> m_shapeLookup;
         std::map<Key, Value<mfem::Vector>> m_physShapeLookup;

         std::map<Key, Value<mfem::DenseMatrix>> m_dShapeLookup;
         std::map<Key, Value<mfem::DenseMatrix>> m_physDShapeLookup;
   };

   /**
    * @brief Base class for operators
    */
   template <ShapeFunctionSpaceType Space>
   class ShapeFunctionBase : public FormLanguage::Base
   {
      public:
         constexpr
         ShapeFunctionBase()
            : FormLanguage::Base()
         {}

         constexpr
         ShapeFunctionBase(const ShapeFunctionBase& other)
            : FormLanguage::Base(other)
         {}

         constexpr
         ShapeFunctionBase(ShapeFunctionBase&& other)
            : FormLanguage::Base(std::move(other))
         {}

         constexpr
         ShapeFunctionSpaceType getSpaceType() const
         {
            return Space;
         }

         constexpr
         RangeType getRangeType() const
         {
            if (getRows() == 1 && getColumns() == 1)
               return RangeType::Scalar;
            else if (getRows() > 1 && getColumns() == 1)
               return RangeType::Vector;
            else
               return RangeType::Matrix;
         }

         constexpr
         RangeShape getRangeShape() const
         {
            return {getRows(), getColumns()};
         }

         virtual const ShapeFunctionBase<Space>& getLeaf() const = 0;

         virtual int getRows() const = 0;

         virtual int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const = 0;

         virtual int getColumns() const = 0;

         virtual void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip,
               ShapeComputator& compute) const = 0;

         Transpose<ShapeFunctionBase<Space>> T() const
         {
            return Transpose<ShapeFunctionBase<Space>>(*this);
         }

         virtual FiniteElementSpaceBase& getFiniteElementSpace() = 0;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;

         virtual ShapeFunctionBase<Space>* copy() const noexcept override = 0;
   };

   template <ShapeFunctionSpaceType Space, class Context>
   class ShapeFunction<H1<Context>, Space> : public ShapeFunctionBase<Space>
   {
      public:
         constexpr
         ShapeFunction(H1<Context>& fes)
            : m_fes(fes)
         {}

         constexpr
         ShapeFunction(const ShapeFunction& other)
            :  ShapeFunctionBase<Space>(other),
               m_fes(other.m_fes)
         {}

         constexpr
         ShapeFunction(ShapeFunction&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_fes(other.m_fes)
         {}

         H1<Context>& getFiniteElementSpace() override
         {
            return m_fes;
         }

         const H1<Context>& getFiniteElementSpace() const override
         {
            return m_fes;
         }

         int getRows() const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation&) const override
         {
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip,
               ShapeComputator& compute
               ) const override
         {
            const auto& shape = compute.getPhysicalShape(fe, trans, ip);
            const int n = shape.Size();
            const int vdim = getFiniteElementSpace().getVectorDimension();
            op.setSize(vdim, 1, vdim * n);
            op = 0.0;
            for (int i = 0; i < vdim; i++)
               for (int j = 0; j < n; j++)
                  op(i, 0, j + i * n) = shape(j);
         }

         virtual const ShapeFunction<H1<Context>, Space>& getLeaf() const override = 0;

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         H1<Context>& m_fes;
   };

   template <ShapeFunctionSpaceType Space, class Context>
   using H1ShapeFunction = ShapeFunction<H1<Context>, Space>;
}

#endif
