#ifndef RODIN_VARIATIONAL_SHAPEFUNCTION_H
#define RODIN_VARIATIONAL_SHAPEFUNCTION_H

#include <mfem.hpp>

#include "Rodin/Alert/Exception.h"
#include "Rodin/FormLanguage/Base.h"

#include "H1.h"
#include "ForwardDecls.h"
#include "BasisOperator.h"
#include "FiniteElementSpace.h"
#include "Integrator.h"
#include "Assembly.h"

#include "RangeShape.h"

namespace Rodin::Variational
{
   /**
    * @defgroup ShapeFunctionSpecializations ShapeFunction Template Specializations
    * @brief Template specializations of the ShapeFunction class.
    * @see ShapeFunction
    */

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

         constexpr
         Transpose<ShapeFunctionBase<Space>> T() const
         {
            return Transpose<ShapeFunctionBase<Space>>(*this);
         }

         virtual const ShapeFunctionBase<Space>& getLeaf() const = 0;

         virtual int getRows() const = 0;

         virtual int getColumns() const = 0;

         virtual int getDOFs(const Geometry::SimplexBase& element) const = 0;

         virtual void getOperator(
               DenseBasisOperator& op,
               ShapeComputator& compute,
               const Geometry::SimplexBase& element) const = 0;

         virtual FiniteElementSpaceBase& getFiniteElementSpace() = 0;

         virtual const FiniteElementSpaceBase& getFiniteElementSpace() const = 0;

         virtual ShapeFunctionBase<Space>* copy() const noexcept override = 0;
   };

   /**
    * @ingroup ShapeFunctionSpecializations
    * @brief L2 ShapeFunction
    */
   template <ShapeFunctionSpaceType Space, class ... Ts>
   class ShapeFunction<L2<Ts...>, Space> : public ShapeFunctionBase<Space>
   {
      public:
         using FES = L2<Ts...>;

         constexpr
         ShapeFunction(FES& fes)
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

         FES& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         const FES& getFiniteElementSpace() const override
         {
            return m_fes.get();
         }

         int getRows() const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(const Geometry::SimplexBase& element) const override
         {
            const auto& fe = getFiniteElementSpace().getFiniteElement(element);
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         void getOperator(
               DenseBasisOperator& op,
               ShapeComputator& compute,
               const Geometry::Point& point,
               const Geometry::Element& element) const override
         {
            const auto& shape =
               compute.getPhysicalShape(
                     getFiniteElementSpace().getFiniteElement(element),
                     element.getTransformation(),
                     element.getTransformation().GetIntPoint());
            const int n = shape.Size();
            const int vdim = getFiniteElementSpace().getVectorDimension();
            op.setSize(vdim, 1, vdim * n);
            op = 0.0;
            for (int i = 0; i < vdim; i++)
               for (int j = 0; j < n; j++)
                  op(i, 0, j + i * n) = shape(j);
         }

         virtual const ShapeFunction<FES, Space>& getLeaf() const override = 0;

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         std::reference_wrapper<FES> m_fes;
   };

   /**
    * @ingroup ShapeFunctionSpecializations
    * @brief H1 ShapeFunction
    */
   template <ShapeFunctionSpaceType Space, class ... Ts>
   class ShapeFunction<H1<Ts...>, Space> : public ShapeFunctionBase<Space>
   {
      public:
         using FES = H1<Ts...>;

         constexpr
         ShapeFunction(FES& fes)
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

         FES& getFiniteElementSpace() override
         {
            return m_fes.get();
         }

         const FES& getFiniteElementSpace() const override
         {
            return m_fes.get();
         }

         int getRows() const override
         {
            return getFiniteElementSpace().getVectorDimension();
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(const Geometry::SimplexBase& element) const override
         {
            const auto& fe = getFiniteElementSpace().getFiniteElement(element);
            return fe.GetDof() * getFiniteElementSpace().getVectorDimension();
         }

         void getOperator(
               DenseBasisOperator& op,
               ShapeComputator& compute,
               const Geometry::SimplexBase& element) const override
         {
            const auto& shape =
               compute.getPhysicalShape(
                     getFiniteElementSpace().getFiniteElement(element),
                     element.getTransformation(),
                     element.getTransformation().GetIntPoint());
            const int n = shape.Size();
            const int vdim = getFiniteElementSpace().getVectorDimension();
            op.setSize(vdim, 1, vdim * n);
            op = 0.0;
            for (int i = 0; i < vdim; i++)
               for (int j = 0; j < n; j++)
                  op(i, 0, j + i * n) = shape(j);
         }

         virtual const ShapeFunction<FES, Space>& getLeaf() const override = 0;

         virtual ShapeFunction* copy() const noexcept override = 0;

      private:
         std::reference_wrapper<FES> m_fes;
   };
}

#endif
