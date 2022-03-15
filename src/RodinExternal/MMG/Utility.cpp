#include <algorithm>

#include <mmg2d/mmg2d.h>
#include <mmg3d/mmg3d.h>

#include "Rodin/Alert.h"

#include "Utility.h"

namespace Rodin::External::MMG
{
  void MMG5_Mesh_Copy(const MMG5_pMesh src, MMG5_pMesh dst)
  {
     assert(src);
     assert(dst);

     // Copy all non-pointer fields
     *dst = *src;

     // Fields which have been deprecated (or seem to have anyways)
     dst->ipar = nullptr;

     // int* : mesh->adja
     if (src->adja)
     {
       MMG5_SAFE_CALLOC(dst->adja, 4 * src->nemax + 5, int, /* No op */);
       std::copy(
           src->adja, src->adja + 4 * src->nemax + 5,
           dst->adja);
     }
     else
       dst->adja = nullptr;

     // int* : mesh->adjt
     if (src->adjt)
     {
       MMG5_SAFE_CALLOC(dst->adjt, 3 * src->nt + 4, int, /* No op */);
       std::copy(
           src->adjt, src->adjt + 3 * src->nt + 4,
           dst->adjt);
     }
     else
       dst->adjt = nullptr;

     // int* : mesh->adjapr
     if (src->adjapr)
     {
       MMG5_SAFE_CALLOC(dst->adjapr, 5 * src->nprism + 6, int, /* No op */);
       std::copy(
           src->adjapr, src->adjapr + 5 * src->nprism + 6,
           dst->adjapr);
     }
     else
       dst->adjapr = nullptr;

     // int* : mesh->adjq
     if (src->adjq)
     {
       MMG5_SAFE_CALLOC(dst->adjq, 4 * src->nquad + 5, int, /* No op */);
       std::copy(
           src->adjq, src->adjq + 4 * src->nquad + 5,
           dst->adjq);
     }
     else
       dst->adjq = nullptr;

     // MMG5_pPoint : mesh->point
     if (src->point)
     {
       /* Why MMG fields are 1-indexed is beyond me, but it does not pose so much
        * difficulty. We will copy from the 0-index, if it contains garbage then
        * it should not matter because MMG code should only access memory from the
        * 1-index.
        */
       MMG5_SAFE_CALLOC(dst->point, src->npmax + 1,
           MMG5_Point, /* No op */);
       std::copy(
           src->point, src->point + src->npmax + 1,
           dst->point);
     }
     else
       dst->point = nullptr;

     // MMG5_pxPoint : mesh->xpoint
     if (src->xpoint)
     {
       MMG5_SAFE_CALLOC(dst->xpoint, src->xpmax + 1, MMG5_xPoint, /* No op */);
       std::copy(
           src->xpoint, src->xpoint + src->xpmax + 1,
           dst->xpoint);
     }
     else
       dst->xpoint = nullptr;

     // MMG5_pTetra : mesh->tetra
     if (src->tetra)
     {
       MMG5_SAFE_CALLOC(dst->tetra, src->nemax + 1, MMG5_Tetra, /* No op */);
       std::copy(
           src->tetra, src->tetra + src->nemax + 1,
           dst->tetra);
     }
     else
       dst->tetra = nullptr;

     // MMG5_pxTetra : mesh->xtetra
     if (src->tetra)
     {
       MMG5_SAFE_CALLOC(dst->xtetra, src->xtmax + 1, MMG5_xTetra, /* No op */);
       std::copy(
           src->xtetra, src->xtetra + src->xtmax + 1,
           dst->xtetra);
     }
     else
       dst->xtetra = nullptr;

     // MMG5_pPrism : mesh->prism
     if (src->prism)
     {
       MMG5_SAFE_CALLOC(dst->prism, src->nprism + 1, MMG5_Prism, /* No op */);
       std::copy(
           src->prism, src->prism + src->nprism + 1,
           dst->prism);
     }
     else
       dst->prism = nullptr;

     // MMG5_pxprism : mesh->xprism
     if (src->xprism)
     {
       MMG5_SAFE_CALLOC(dst->xprism, src->nprism + 1, MMG5_xPrism, /* No op */);
       std::copy(
           src->xprism, src->xprism + src->nprism + 1,
           dst->xprism);
     }
     else
       dst->xprism = nullptr;

     // MMG5_pTria : mesh->tria
     if (src->tria)
     {
       MMG5_SAFE_CALLOC(dst->tria, src->ntmax + 1, MMG5_Tria, /* No op */);
       std::copy(
           src->tria, src->tria + src->ntmax + 1,
           dst->tria);
     }
     else
       dst->tria = nullptr;

     // MMG5_pQuad : mesh->quadra
     if (src->quadra)
     {
       MMG5_SAFE_CALLOC(dst->quadra, src->nquad + 1, MMG5_Quad, /* No op */);
       std::copy(
           src->quadra, src->quadra + src->nquad + 1,
           dst->quadra);
     }
     else
       dst->quadra = nullptr;

     // MMG5_pEdge : mesh->edge
     if (src->edge)
     {
       MMG5_SAFE_CALLOC(dst->edge, src->namax + 1, MMG5_Edge, /* No op */);
       std::copy(
           src->edge, src->edge + src->namax + 1,
           dst->edge);
     }
     else
       dst->edge = nullptr;

     // MMG5_HGeom : mesh->htab
     dst->htab = src->htab;

     // MMG5_hgeom : mesh->htab.geom
     if (src->htab.geom)
     {
       MMG5_SAFE_CALLOC(dst->htab.geom, src->htab.max + 1, MMG5_hgeom, /* No op */);
       std::copy(
           src->htab.geom, src->htab.geom + src->htab.max + 1, dst->htab.geom);
     }
     else
       dst->htab.geom = nullptr;

     // MMG5_Info : mesh->info
     dst->info = src->info;

     // MMG5_Par : mesh->par
     if (src->info.par)
     {
       MMG5_SAFE_CALLOC(dst->info.par, src->info.npar, MMG5_Par, /* No op */);
       std::copy(src->info.par, src->info.par + src->info.npar, dst->info.par);
     }
     else
       dst->info.par = nullptr;

     // int* : mesh->info.br
     if (src->info.br)
     {
       MMG5_SAFE_CALLOC(dst->info.br, src->info.nbr, int, /* No op */);
       std::copy(src->info.br, src->info.br + src->info.nbr, dst->info.br);
     }
     else
       dst->info.br = nullptr;

     // MMG5_Mat : mesh->info.mat
     if (src->info.mat)
     {
       MMG5_SAFE_CALLOC(dst->info.mat, src->info.nmat, MMG5_Mat, /* No op */);
       std::copy(src->info.mat, src->info.mat + src->info.nmat, dst->info.mat);
     }
     else
       dst->info.mat = nullptr;

     // MMG5_Info : mesh->info.invmat.lookup
     if (src->info.invmat.lookup)
     {
       MMG5_SAFE_CALLOC(dst->info.invmat.lookup, src->info.invmat.size, int, /* No op */);
       std::copy(
           src->info.invmat.lookup,
           src->info.invmat.lookup + src->info.invmat.size,
           dst->info.invmat.lookup
           );
     }
     else
       dst->info.invmat.lookup = nullptr;

     // char* : mesh->namein
     MMG5_SAFE_CALLOC(dst->namein, std::strlen(src->namein),
         char, /* No op */);
     std::strcpy(dst->namein, src->namein);

     // char* : mesh->nameout
     MMG5_SAFE_CALLOC(dst->nameout, std::strlen(src->nameout),
         char, /* No op */);
     std::strcpy(dst->nameout, src->nameout);
  }

  void MMG5_Mesh_To_MFEM_Mesh_Cast(const MMG5_pMesh src, mfem::Mesh& dst)
  {
      bool shiftTriAttr = false;
      for (int i = 1; i <= src->nt; i++)
      {
        if (src->tria[i].ref == 0)
          shiftTriAttr = true;
        if (src->tria[i].ref < 0)
          Alert::Exception(
              "Negative triangle attributes are not supported").raise();
      }

      bool shiftEdgeAttr = false;
      for (int i = 1; i <= src->na; i++)
      {
        if (src->edge[i].ref == 0)
          shiftEdgeAttr = true;
        if (src->edge[i].ref < 0)
          Alert::Exception(
              "Negative edge element attributes are not supported.").raise();
      }

      if (shiftTriAttr)
        Alert::Warning(
            "Triangles with attribute equal to 0 are not supported. "
            "All triangle attributes will be incremented by 1.").raise();

      if (shiftEdgeAttr)
        Alert::Warning(
            "Edges with attribute equal to 0 are not supported. "
            "All edge attributes will be incremented by 1.").raise();

      /* So for some reason mmg types are 1 indexed. So when accessing the
       * arrays make sure to start at 1 and not 0. I don't know why this is the
       * case and I'm not sure if it's for every array in the library.
       */
      for (int i = 1; i <= src->np; i++)
      {
        dst.AddVertex(
            src->point[i].c[0],
            src->point[i].c[1],
            src->point[i].c[2]
            );
      }

      for (int i = 1; i <= src->nt; i++)
      {
        dst.AddTriangle(
            src->tria[i].v[0] - 1,
            src->tria[i].v[1] - 1,
            src->tria[i].v[2] - 1,
            src->tria[i].ref + shiftTriAttr
            );
      }

      for (int i = 1; i <= src->na; i++)
      {
         dst.AddBdrSegment(
               src->edge[i].a - 1,
               src->edge[i].b - 1,
               src->edge[i].ref + shiftEdgeAttr
               );
      }
      dst.FinalizeMesh(0, true);
  }
}
