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
    * @defgroup JacobianSpecializations Jacobian Template Specializations
    * @brief Template specializations of the Jacobian class.
    * @see Jacobian
    */

   /**
    * @ingroup JacobianSpecializations
    * @brief Jacobian of an H1 GridFunction object.
    */
   template <class Trait>
   class Jacobian<GridFunction<H1<Trait>>> : public MatrixFunctionBase
   {
      public:
         /**
          * @brief Constructs the Jacobian matrix of an @f$ H^1 (\Omega)^d @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Jacobian(GridFunction<H1<Trait>>& u)
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

         FunctionValue getValue(const Geometry::Point& p) const override
         {
            FunctionValue::Matrix grad;
            const auto& element = p.getSimplex();
            const auto& mesh = m_u.getFiniteElementSpace().getMesh();
            if (element.getDimension() == mesh.getDimension())
            {
               assert(dynamic_cast<const Geometry::Element*>(&p.getSimplex()));
               const auto& element = p.getSimplex();
               auto& trans = element.getTransformation();
               m_u.getHandle().GetVectorGradient(trans, grad);
            }
            else if (element.getDimension() == mesh.getDimension() - 1)
            {
               assert(dynamic_cast<const Geometry::Face*>(&p.getSimplex()));
               const auto& face = static_cast<const Geometry::Face&>(p.getSimplex());
               if (face.isInterface())
               {
                  assert(false);
               }
               else if (face.isBoundary())
               {
                  assert(dynamic_cast<const Geometry::Boundary*>(&p.getSimplex()));
                  const auto& boundary = static_cast<const Geometry::Boundary&>(p.getSimplex());
                  auto& ft = boundary.getFaceTransformations();
                  assert(ft.Elem1);
                  auto& trans1 = ft.GetElement1Transformation();
                  assert(ft.Elem2);
                  auto& trans2 = ft.GetElement2Transformation();
                  const auto& mesh = p.getSimplex().getMesh();
                  if (mesh.isSubMesh())
                  {
                     const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(mesh);
                     const auto& parent = submesh.getParent();

                     if (getTraceDomain().count(trans1.Attribute))
                     {
                        int parentIdx = submesh.getElementMap().left.at(trans1.ElementNo);
                        trans1.ElementNo = parentIdx;
                        m_u.getHandle().GetVectorGradient(trans1, grad);
                     }
                     else if (getTraceDomain().count(trans2.Attribute))
                     {
                        int parentIdx = submesh.getElementMap().left.at(trans2.ElementNo);
                        const auto& parentElement = parent.getElement(parentIdx);
                        assert(!parentElement.end());
                        m_u.getHandle().GetVectorGradient(parentElement->getTransformation(), grad);
                     }
                     else
                        assert(false);
                  }
                  else
                  {
                     if (getTraceDomain().count(trans1.Attribute))
                        m_u.getHandle().GetVectorGradient(trans1, grad);
                     else if (getTraceDomain().count(trans2.Attribute))
                        m_u.getHandle().GetVectorGradient(trans2, grad);
                     else
                        assert(false);
                  }
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
            return grad;
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }

      private:
         GridFunction<H1<Trait>>& m_u;
   };
   template <class Trait>
   Jacobian(GridFunction<H1<Trait>>&) -> Jacobian<GridFunction<H1<Trait>>>;

   /**
    * @ingroup JacobianSpecializations
    * @brief Jacobian of an H1 ShapeFunction object.
    */
   template <ShapeFunctionSpaceType Space, class Trait>
   class Jacobian<ShapeFunction<H1<Trait>, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Jacobian(ShapeFunction<H1<Trait>, Space>& u)
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

         H1<Trait>& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const H1<Trait>& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         const ShapeFunction<H1<Trait>, Space>& getLeaf() const override
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

         int getDOFs(const Geometry::Simplex& element) const override
         {
            return m_u.getDOFs(element);
         }

         void getOperator(
               DenseBasisOperator& op, ShapeComputator& compute,
               const Geometry::Point& p) const override
         {
            const auto& element = p.getSimplex();
            auto& trans = element.getTransformation();
            const auto& fe = getFiniteElementSpace().getFiniteElement(element);
            const auto& dshape = compute.getPhysicalDShape(fe, trans, trans.GetIntPoint());
            const int n = dshape.NumRows();
            const int sdim = trans.GetSpaceDim();
            const int vdim = m_u.getFiniteElementSpace().getVectorDimension();
            op.setSize(sdim, vdim, vdim * n);
            op = 0.0;
            for (int i = 0; i < vdim; i++)
               for (int j = 0; j < n; j++)
                  for (int k = 0; k < sdim; k++)
                     op(k, i, j + i * n) = dshape(j, k);
         }

         Jacobian* copy() const noexcept override
         {
            return new Jacobian(*this);
         }
      private:
         ShapeFunction<H1<Trait>, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space, class Trait>
   Jacobian(ShapeFunction<H1<Trait>, Space>&) -> Jacobian<ShapeFunction<H1<Trait>, Space>>;
}

#endif
