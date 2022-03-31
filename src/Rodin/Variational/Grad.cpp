#include "Grad.h"

namespace Rodin::Variational
{
   mfem::ElementTransformation* Grad<GridFunction<H1>>::RefinedToCoarse(
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
}
