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
            const auto& simplex = p.getSimplex();
            const auto& simplexMesh = simplex.getMesh();
            const auto& fesMesh = m_u.getFiniteElementSpace().getMesh();
            if (simplex.getDimension() == fesMesh.getDimension())
            {
               assert(dynamic_cast<const Geometry::Element*>(&p.getSimplex()));
               const auto& element = p.getSimplex();
               auto& trans = element.getTransformation();
               m_u.getHandle().GetVectorGradient(trans, grad);
            }
            else if (simplex.getDimension() == fesMesh.getDimension() - 1)
            {
               assert(dynamic_cast<const Geometry::Face*>(&p.getSimplex()));
               const auto& face = static_cast<const Geometry::Face&>(p.getSimplex());
               mfem::FaceElementTransformations* ft =
                  const_cast<Geometry::MeshBase&>(simplexMesh).getHandle()
                  .GetFaceElementTransformations(face.getIndex());
               if (simplexMesh.isSubMesh())
               {
                  const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(simplexMesh);
                  assert(submesh.getParent() == fesMesh);
                  if (ft->Elem1 && getTraceDomain() == ft->Elem1->Attribute)
                  {
                     Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
                     ft->Elem1->ElementNo = parentIdx;
                     ft->Elem1No = parentIdx;
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
                  }
                  else if (ft->Elem2 && getTraceDomain() == ft->Elem2->Attribute)
                  {
                     Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem2No);
                     ft->Elem2->ElementNo = parentIdx;
                     ft->Elem2No = parentIdx;
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem2, grad);
                  }
                  else if (face.isBoundary())
                  {
                     assert(ft->Elem1);
                     Geometry::Index parentIdx = submesh.getElementMap().left.at(ft->Elem1No);
                     ft->Elem1->ElementNo = parentIdx;
                     ft->Elem1No = parentIdx;
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
                  }
                  else
                  {
                     assert(false);
                  }
               }
               else if (fesMesh.isSubMesh())
               {
                  const auto& submesh = static_cast<const Geometry::SubMesh<Context::Serial>&>(fesMesh);
                  assert(submesh.getParent() == simplexMesh);
                  const auto& s2pe = submesh.getElementMap();
                  if (ft->Elem1 && s2pe.right.count(ft->Elem1No) && getTraceDomain() == ft->Elem1->Attribute)
                  {
                     Geometry::Index idx = s2pe.right.at(ft->Elem1No);
                     ft->Elem1->ElementNo = idx;
                     ft->Elem1No = idx;
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
                  }
                  else if (ft->Elem2 && s2pe.right.count(ft->Elem2No) && getTraceDomain() == ft->Elem2->Attribute)
                  {
                     Geometry::Index idx = s2pe.right.at(ft->Elem2No);
                     ft->Elem2->ElementNo = idx;
                     ft->Elem2No = idx;
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem2, grad);
                  }
                  else if (face.isBoundary())
                  {
                     assert(ft->Elem1);
                     Geometry::Index idx = s2pe.right.at(ft->Elem1No);
                     ft->Elem1->ElementNo = idx;
                     ft->Elem1No = idx;
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
                  }
                  else
                  {
                     assert(false);
                  }
               }
               else
               {
                  if (ft->Elem1 && getTraceDomain() == ft->Elem1->Attribute)
                  {
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
                  }
                  else if (ft->Elem2 && getTraceDomain() == ft->Elem2->Attribute)
                  {
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     m_u.getHandle().GetVectorGradient(*ft->Elem2, grad);
                  }
                  else if (face.isBoundary())
                  {
                     ft->SetAllIntPoints(&p.getIntegrationPoint());
                     assert(ft->Elem1);
                     m_u.getHandle().GetVectorGradient(*ft->Elem1, grad);
                  }
                  else
                  {
                     assert(false);
                  }
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
