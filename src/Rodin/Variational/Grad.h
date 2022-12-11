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
#include "Utility.h"
#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
   /**
    * @defgroup GradSpecializations Grad Template Specializations
    * @brief Template specializations of the Grad class.
    * @see Grad
    */

   /**
    * @ingroup GradSpecializations
    */
   template <class Trait>
   class Grad<GridFunction<H1<Trait>>> : public VectorFunctionBase
   {
      public:
         /**
          * @brief Constructs the gradient of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Grad(const GridFunction<H1<Trait>>& u)
            : m_u(u)
         {}

         Grad(const Grad& other)
            :  VectorFunctionBase(other),
               m_u(other.m_u)
         {}

         Grad(Grad&& other)
            :  VectorFunctionBase(std::move(other)),
               m_u(other.m_u)
         {}

         int getDimension() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getSpaceDimension();
         }

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            FunctionValue::Vector grad;
            switch (p.getElement().getRegion())
            {
               case Geometry::Region::Domain:
               {
                  assert(dynamic_cast<const Geometry::Element*>(&p.getElement()));
                  const auto& element = p.getElement();
                  auto& trans = element.getTransformation();
                  m_u.getHandle().GetGradient(trans, grad);
                  break;
               }
               case Geometry::Region::Boundary:
               {
                  assert(dynamic_cast<const Geometry::Boundary*>(&p.getElement()));
                  const auto& boundary = static_cast<const Geometry::Boundary&>(p.getElement());
                  auto& ft = boundary.getFaceTransformations();
                  assert(ft.Elem1);
                  auto& trans1 = ft.GetElement1Transformation();
                  assert(ft.Elem2);
                  auto& trans2 = ft.GetElement2Transformation();
                  const auto& mesh = p.getElement().getMesh();
                  if (&mesh == &m_u.getFiniteElementSpace().getMesh())
                  {
                     if (getTraceDomain().count(trans1.Attribute))
                        m_u.getHandle().GetGradient(trans1, grad);
                     else if (getTraceDomain().count(trans2.Attribute))
                        m_u.getHandle().GetGradient(trans2, grad);
                     else
                        assert(false);
                  }
                  else if (mesh.isSubMesh())
                  {
                     // Element's mesh is the child
                     const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(mesh);
                     const auto& parent = submesh.getParent();

                     if (&parent == &m_u.getFiniteElementSpace().getMesh())
                     {
                        submesh.getElementMap().left.at(trans1.ElementNo);

                        if (getTraceDomain().count(trans1.Attribute))
                        {
                           int parentIdx = submesh.getElementMap().left.at(trans1.ElementNo);
                           trans1.ElementNo = parentIdx;
                           m_u.getHandle().GetGradient(trans1, grad);
                        }
                        else if (getTraceDomain().count(trans2.Attribute))
                        {
                           int parentIdx = submesh.getElementMap().left.at(trans2.ElementNo);
                           trans2.ElementNo = parentIdx;
                           m_u.getHandle().GetGradient(trans2, grad);
                        }
                        else
                           assert(false);
                     }
                     else
                     {
                        assert(false);
                     }
                  }
                  else if (m_u.getFiniteElementSpace().getMesh().isSubMesh())
                  {
                     // Element's mesh is the parent

                     const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(
                           m_u.getFiniteElementSpace().getMesh());
                     const auto& parent = submesh.getParent();
                     if (&parent == &mesh)
                     {
                        auto it1 = submesh.getElementMap().right.find(trans1.ElementNo);
                        auto it2 = submesh.getElementMap().right.find(trans2.ElementNo);

                        if (it1 != submesh.getElementMap().right.end() &&
                              (getTraceDomain().size() == 0 || getTraceDomain().count(trans1.Attribute)))
                        {
                           int parentIdx = submesh.getElementMap().left.at(trans1.ElementNo);
                           trans1.ElementNo = parentIdx;
                           m_u.getHandle().GetGradient(trans1, grad);
                        }
                        else if (it2 != submesh.getElementMap().right.end() &&
                              (getTraceDomain().size() == 0 || getTraceDomain().count(trans2.Attribute)))
                        {
                           int parentIdx = submesh.getElementMap().left.at(trans2.ElementNo);
                           trans2.ElementNo = parentIdx;
                           m_u.getHandle().GetGradient(trans2, grad);
                        }
                        else
                           assert(false);
                     }
                     else
                     {
                        assert(false);
                     }
                  }
                  else
                  {
                     assert(false);
                  }
                  break;
               }
               case Geometry::Region::Interface:
               {
                  assert(false);
                  break;
               }
            }
            return grad;
         }

         VectorFunctionBase* copy() const noexcept override
         {
            return new Grad(*this);
         }

      private:
         const GridFunction<H1<Trait>>& m_u;
   };
   template <class Trait>
   Grad(const GridFunction<H1<Trait>>&) -> Grad<GridFunction<H1<Trait>>>;

   /**
    * @ingroup GradSpecializations
    */
   template <ShapeFunctionSpaceType Space, class ... Ts>
   class Grad<ShapeFunction<H1<Ts...>, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         using FES = H1<Ts...>;

         constexpr
         Grad(ShapeFunction<H1<Ts...>, Space>& u)
            : m_u(u)
         {
            if (u.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, u.getRangeType()).raise();
         }

         constexpr
         Grad(const Grad& other)
            : ShapeFunctionBase<Space>(other),
              m_u(other.m_u)
         {}

         constexpr
         Grad(Grad&& other)
            : ShapeFunctionBase<Space>(std::move(other)),
              m_u(other.m_u)
         {}

         const ShapeFunction<H1<Ts...>, Space>& getLeaf() const override
         {
            return m_u.getLeaf();
         }

         int getRows() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getSpaceDimension();
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(const Geometry::Simplex& element) const override
         {
            return m_u.getDOFs(element);
         }

         void getOperator(
               DenseBasisOperator& op,
               ShapeComputator& compute,
               const Geometry::Point& p) const override
         {
            auto& trans = p.getElement().getTransformation();
            const auto& fe = getFiniteElementSpace().getFiniteElement(p.getElement());
            const auto& dshape = compute.getPhysicalDShape(fe, trans, trans.GetIntPoint());
            const int n = dshape.NumRows();
            const int sdim = trans.GetSpaceDim();
            op.setSize(sdim, 1, n);
            for (int j = 0; j < n; j++)
               for (int k = 0; k < sdim; k++)
                  op(k, 0, j) = dshape(j, k);
         }

         H1<Ts...>& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const H1<Ts...>& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Grad* copy() const noexcept override
         {
            return new Grad(*this);
         }
      private:
         ShapeFunction<H1<Ts...>, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space, class ... Ts>
   Grad(ShapeFunction<H1<Ts...>, Space>&) -> Grad<ShapeFunction<H1<Ts...>, Space>>;

   /**
    * @ingroup GradSpecializations
    *
    * @brief Represents the broken gradient.
    */
   template <ShapeFunctionSpaceType Space, class ... Ts>
   class Grad<ShapeFunction<L2<Ts...>, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         using FES = L2<Ts...>;

         constexpr
         Grad(ShapeFunction<FES, Space>& u)
            : m_u(u)
         {
            if (u.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, u.getRangeType()).raise();
         }

         constexpr
         Grad(const Grad& other)
            : ShapeFunctionBase<Space>(other),
              m_u(other.m_u)
         {}

         constexpr
         Grad(Grad&& other)
            : ShapeFunctionBase<Space>(std::move(other)),
              m_u(other.m_u)
         {}

         const ShapeFunction<FES, Space>& getLeaf() const override
         {
            return m_u.getLeaf();
         }

         int getRows() const override
         {
            return m_u.getFiniteElementSpace().getMesh().getSpaceDimension();
         }

         int getColumns() const override
         {
            return 1;
         }

         int getDOFs(const Geometry::Simplex& element) const override
         {
            return m_u.getDOFs(element);
         }

         void getOperator(
               DenseBasisOperator& op,
               ShapeComputator& compute,
               const Geometry::Simplex& element) const override
         {
            assert(element.getRegion() == Geometry::Region::Domain);
            auto& trans = element.getTransformation();
            const auto& fe = getFiniteElementSpace().getFiniteElement(element);
            const auto& dshape = compute.getPhysicalDShape(fe, trans, trans.GetIntPoint());
            const int n = dshape.NumRows();
            const int sdim = trans.GetSpaceDim();
            op.setSize(sdim, 1, n);
            for (int j = 0; j < n; j++)
               for (int k = 0; k < sdim; k++)
                  op(k, 0, j) = dshape(j, k);
         }

         FES& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const FES& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Grad* copy() const noexcept override
         {
            return new Grad(*this);
         }
      private:
         ShapeFunction<FES, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space, class ... Ts>
   Grad(ShapeFunction<L2<Ts...>, Space>&) -> Grad<ShapeFunction<L2<Ts...>, Space>>;
}

#endif
