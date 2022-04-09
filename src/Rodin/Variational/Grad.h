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
#include "Jacobian.h"
#include "GridFunction.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "VectorFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represents the gradient @f$ \nabla u @f$ of the scalar function
    * @f$ u @f$.
    *
    * For @f$ u : \mathbb{R}^n \rightarrow \mathbb{R} @f$, the gradient
    * @f$ \nabla u : \mathbb{R}^n \rightarrow \mathbb{R} @f$ at the point
    * @f$ x = (x_1, \ldots, x_n) @f$ is defined by
    * @f[
    *    \nabla u (x) =
    *    \left[
    *       \dfrac{\partial u}{\partial x_1}(x), \ldots,
    *       \dfrac{\partial u}{\partial x_n}(x)
    *    \right]^T
    * @f]
    */
   template <class Trait>
   class Grad<GridFunction<H1, Trait>> : public VectorFunctionBase
   {

      mfem::ElementTransformation *RefinedToCoarse(
        mfem::Mesh &coarse_mesh, const mfem::ElementTransformation &T,
        const mfem::IntegrationPoint &ip, mfem::IntegrationPoint &coarse_ip) const
      {
        mfem::Mesh &fine_mesh = *T.mesh;
        // Get the element transformation of the coarse element containing the
        // fine element.
        int fine_element = T.ElementNo;
        const mfem::CoarseFineTransformations &cf = fine_mesh.GetRefinementTransforms();
        int coarse_element = cf.embeddings[fine_element].parent;
        mfem::ElementTransformation *coarse_T =
           coarse_mesh.GetElementTransformation(coarse_element);
        // Transform the integration point from fine element coordinates to coarse
        // element coordinates.
        mfem::Geometry::Type geom = T.GetGeometryType();
        mfem::IntegrationPointTransformation fine_to_coarse;
        mfem::IsoparametricTransformation &emb_tr = fine_to_coarse.Transf;
        emb_tr.SetIdentityTransformation(geom);
        emb_tr.SetPointMat(cf.point_matrices[geom](cf.embeddings[fine_element].matrix));
        fine_to_coarse.Transform(ip, coarse_ip);
        coarse_T->SetIntPoint(&coarse_ip);
        return coarse_T;
      }

      void GetGradient(
            mfem::Vector& grad, mfem::ElementTransformation& trans,
            const mfem::IntegrationPoint& ip) const
      {
         mfem::Mesh* gf_mesh = m_u.getHandle().FESpace()->GetMesh();
         if (trans.mesh == gf_mesh)
         {
            m_u.getHandle().GetGradient(trans, grad);
         }
         else
         {
            mfem::IntegrationPoint coarse_ip;
            mfem::ElementTransformation *coarse_T = RefinedToCoarse(*gf_mesh, trans, ip, coarse_ip);
            m_u.getHandle().GetGradient(*coarse_T, grad);
         }
      }

      public:
         /**
          * @brief Constructs the gradient of an @f$ H^1 @f$ function
          * @f$ u @f$.
          * @param[in] u Grid function to be differentiated
          */
         Grad(const GridFunction<H1, Trait>& u)
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

         void getValue(
               mfem::Vector& value,
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
                        GetGradient(value, ft->GetElement1Transformation(), ip);
                     else if (traceDomain.count(ft->GetElement2Transformation().Attribute))
                        GetGradient(value, ft->GetElement2Transformation(), ip);
                     else
                     {
                        // The boundary over which we are evaluating must be
                        // the interface between the trace domain and some
                        // other domain, i.e. it is not the boundary that was
                        // specified!
                        Alert::Exception()
                           << "Boundary element " << trans.ElementNo
                           << " with attribute " << trans.Attribute
                           << " is not a boundary of the trace domain."
                           << Alert::Raise;
                     }
                  }
                  else
                  {
                     GetGradient(value, trans, ip);
                  }
                  break;
               }
               default:
                  GetGradient(value, trans, ip);
            }
         }

         VectorFunctionBase* copy() const noexcept override
         {
            return new Grad(*this);
         }

      private:
         const GridFunction<H1, Trait>& m_u;
   };
   template <class Trait>
   Grad(const GridFunction<H1, Trait>&) -> Grad<GridFunction<H1, Trait>>;

   template <ShapeFunctionSpaceType Space>
   class Grad<ShapeFunction<H1, Space>> : public ShapeFunctionBase<Space>
   {
      public:
         Grad(ShapeFunction<H1, Space>& u)
            : m_u(u)
         {
            assert(m_u.getFiniteElementSpace().getVectorDimension() == 1);
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
               const mfem::ElementTransformation& trans) const override
         {
            return trans.GetSpaceDim();
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
            mfem::DenseMatrix dshape;
            dshape.SetSize(dofs, sdim);
            fe.CalcPhysDShape(trans, dshape);
            return std::unique_ptr<Rank3Operator>(
                  new Internal::JacobianShapeR3O(std::move(dshape), sdim, 1));
         }

         FiniteElementSpace<H1>& getFiniteElementSpace() override
         {
            return m_u.getFiniteElementSpace();
         }

         const FiniteElementSpace<H1>& getFiniteElementSpace() const override
         {
            return m_u.getFiniteElementSpace();
         }

         Grad* copy() const noexcept override
         {
            return new Grad(*this);
         }
      private:
         ShapeFunction<H1, Space>& m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Grad(ShapeFunction<H1, Space>&) -> Grad<ShapeFunction<H1, Space>>;
}

#endif
