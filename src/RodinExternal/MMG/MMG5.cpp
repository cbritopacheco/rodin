/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/IO/MeshLoader.h"

#include "MMG5.h"

namespace Rodin::External::MMG
{
   MMG5_pMesh MMG5::createMesh(int version, int dim, std::optional<int> spaceDim)
   {
      MMG5_pMesh res;
      MMG5_SAFE_CALLOC(res, 1, MMG5_Mesh,
            Alert::Exception("Failed to allocate memory for the mesh").raise());

      if (!spaceDim)
         spaceDim = dim;
      bool isSurface = (*spaceDim - 1) == dim;
      assert(isSurface || (dim == *spaceDim));
      if (isSurface)
      {
         MMGS_Set_commonFunc();
         MMGS_Init_parameters(res);
      }
      else if (dim == 2)
      {
         MMG2D_Set_commonFunc();
         MMG2D_Init_parameters(res);
      }
      else if (dim == 3)
      {
         MMG3D_Set_commonFunc();
         MMG3D_Init_parameters(res);
      }
      else
      {
         Alert::Exception("Unhandled case").raise();
         return nullptr;
      }
      res->np = 0;
      res->ver = version;
      res->dim = dim;
      res->npmax = std::max({MMG2D_NPMAX, MMG3D_NPMAX, MMGS_NPMAX});
      res->info.imprim = getMMGVerbosityLevel();
      return res;
   }

   bool MMG5::isSurfaceMesh(MMG5_pMesh mesh)
   {
      switch (mesh->dim)
      {
         case 2:
            return mesh->nt == 0;
         case 3:
            return mesh->ne == 0;
         default:
            return false;
      }
   }

   void MMG5::copyMesh(const MMG5_pMesh src, MMG5_pMesh dst)
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
        MMG5_SAFE_CALLOC(dst->adja, 4 * src->ne + 5, int, /* No op */);
        std::copy(src->adja, src->adja + 4 * src->ne + 5, dst->adja);
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
        MMG5_SAFE_CALLOC(dst->point, src->np + 1, MMG5_Point, /* No op */);
        std::copy(src->point, src->point + src->np + 1, dst->point);
      }
      else
        dst->point = nullptr;

      // MMG5_pxPoint : mesh->xpoint
      if (src->xpoint)
      {
        MMG5_SAFE_CALLOC(dst->xpoint, src->xp + 1, MMG5_xPoint, /* No op */);
        std::copy(src->xpoint, src->xpoint + src->xp + 1, dst->xpoint);
      }
      else
        dst->xpoint = nullptr;

      // MMG5_pTetra : mesh->tetra
      if (src->tetra)
      {
        MMG5_SAFE_CALLOC(dst->tetra, src->ne + 1, MMG5_Tetra, /* No op */);
        std::copy(src->tetra, src->tetra + src->ne + 1, dst->tetra);
      }
      else
        dst->tetra = nullptr;

      // MMG5_pxTetra : mesh->xtetra
      if (src->tetra)
      {
        MMG5_SAFE_CALLOC(dst->xtetra, src->xt + 1, MMG5_xTetra, /* No op */);
        std::copy(src->xtetra, src->xtetra + src->xt + 1, dst->xtetra);
      }
      else
        dst->xtetra = nullptr;

      // MMG5_pPrism : mesh->prism
      if (src->prism)
      {
        MMG5_SAFE_CALLOC(dst->prism, src->nprism + 1, MMG5_Prism, /* No op */);
        std::copy(src->prism, src->prism + src->nprism + 1, dst->prism);
      }
      else
        dst->prism = nullptr;

      // MMG5_pxprism : mesh->xprism
      if (src->xprism)
      {
        MMG5_SAFE_CALLOC(dst->xprism, src->nprism + 1, MMG5_xPrism, /* No op */);
        std::copy(src->xprism, src->xprism + src->nprism + 1, dst->xprism);
      }
      else
        dst->xprism = nullptr;

      // MMG5_pTria : mesh->tria
      if (src->tria)
      {
        MMG5_SAFE_CALLOC(dst->tria, src->nt + 1, MMG5_Tria, /* No op */);
        std::copy(src->tria, src->tria + src->nt + 1, dst->tria);
      }
      else
        dst->tria = nullptr;

      // MMG5_pQuad : mesh->quadra
      if (src->quadra)
      {
        MMG5_SAFE_CALLOC(dst->quadra, src->nquad + 1, MMG5_Quad, /* No op */);
        std::copy(src->quadra, src->quadra + src->nquad + 1, dst->quadra);
      }
      else
        dst->quadra = nullptr;

      // MMG5_pEdge : mesh->edge
      if (src->edge)
      {
        MMG5_SAFE_CALLOC(dst->edge, src->na + 1, MMG5_Edge, /* No op */);
        std::copy(src->edge, src->edge + src->na + 1, dst->edge);
      }
      else
        dst->edge = nullptr;

      // MMG5_HGeom : mesh->htab
      dst->htab = src->htab;

      // MMG5_hgeom : mesh->htab.geom
      if (src->htab.geom)
      {
        MMG5_SAFE_CALLOC(dst->htab.geom, src->htab.max + 1, MMG5_hgeom, /* No op */);
        std::copy(src->htab.geom, src->htab.geom + src->htab.max + 1, dst->htab.geom);
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
      MMG5_SAFE_CALLOC(dst->namein, std::strlen(src->namein), char, /* No op */);
      std::strcpy(dst->namein, src->namein);

      // char* : mesh->nameout
      MMG5_SAFE_CALLOC(dst->nameout, std::strlen(src->nameout), char, /* No op */);
      std::strcpy(dst->nameout, src->nameout);
   }

   MMG5_pMesh MMG5::loadMesh(const boost::filesystem::path& filename)
   {
      Mesh mesh;
      mesh.load(filename, IO::MeshFormat::MEDIT);
      return rodinToMesh(mesh);
   }

   void MMG5::saveMesh(MMG5_pMesh mesh, const boost::filesystem::path& filename)
   {
      // TODO: Change this so we use the MMG2D, MMG3D, MMGS functions
      meshToRodin(mesh).save(filename);
   }

   MMG5_pSol MMG5::createSolution(MMG5_pMesh mesh, int vdim)
   {
      if (vdim < 1 || vdim > 3)
         return nullptr;
      MMG5_pSol res = nullptr;
      MMG5_SAFE_CALLOC(res, 1, MMG5_Sol,
            Alert::Exception("Failed to allocate memory for MMG5_pSol").raise());
      res->ver  = mesh->ver;
      res->size = vdim;
      if (vdim == 1)
         res->type = MMG5_Scalar;
      else
         res->type = MMG5_Vector;
      res->dim = mesh->dim;
      res->np  = mesh->np;
      res->npi = 0;
      res->npmax = mesh->npmax;
      if (res->np)
      {
         // So (res->size + 1) * (res->np + 1) seems to work for most
         // applications
         MMG5_SAFE_CALLOC(res->m, (res->size + 1) * (res->np + 1), double,
               Alert::Exception("Failed to allocate memory for MMG5_pSol->m").raise());
      }
      else
      {
         res->m = nullptr;
      }
      return res;
   }

   MMG5_pMesh MMG5::rodinToMesh(const Rodin::Mesh<Traits::Serial>& src)
   {
      bool isSurface = src.isSurface();
      assert(isSurface || (src.getSpaceDimension() == src.getDimension()));

      MMG5_pMesh res = createMesh(
            s_meshVersionFormatted, src.getDimension(), src.getSpaceDimension());

      res->npi = 0;
      res->nai = 0;
      res->nti = 0;
      res->xp  = 0;

      res->np = src.getHandle().GetNV();

      if (isSurface) // Use MMGS
      {
         res->nt = src.getHandle().GetNE();
         res->na = src.getHandle().GetNBE();

         MMGS_Set_commonFunc();
         if (!MMGS_zaldy(res))
            Alert::Exception("Memory allocation for MMG5_pMesh failed.").raise();

         // Copy points
         for (int i = 1; i <= res->np; i++)
         {
           const double* coords = src.getHandle().GetVertex(i - 1);
           std::copy(coords, coords + 3, res->point[i].c);
         }

         // Copy edges
         for (int i = 1; i <= res->na; i++)
         {
            mfem::Array<int> vertices;
            src.getHandle().GetBdrElementVertices(i - 1, vertices);
            res->edge[i].a = vertices[0] + 1;
            res->edge[i].b = vertices[1] + 1;
            res->edge[i].ref = src.getHandle().GetBdrAttribute(i - 1);
         }

         // Copy triangles
         for (int i = 1; i <= res->nt; i++)
         {
           MMG5_pTria pt = &res->tria[i];
           mfem::Array<int> vertices;
           src.getHandle().GetElementVertices(i - 1, vertices);
           for (int j = 0; j < 3; j++)
             pt->v[j] = vertices[j] + 1;

           pt->ref = src.getHandle().GetAttribute(i - 1);
           for (int j = 0; j < 3; j++)
              res->point[pt->v[j]].tag &= ~MG_NUL;
         }
      }
      else if (src.getDimension() == 2) // Use MMG2D
      {
         res->nt = src.getHandle().GetNE();
         res->na = src.getHandle().GetNBE();

         MMG2D_Set_commonFunc();
         if (!MMG2D_zaldy(res))
            Alert::Exception("Memory allocation for MMG5_pMesh failed.").raise();

         // Copy points
         for (int i = 1; i <= res->np; i++)
         {
           const double* coords = src.getHandle().GetVertex(i - 1);
           std::copy(coords, coords + 3, res->point[i].c);
         }

         // Copy edges
         for (int i = 1; i <= res->na; i++)
         {
            mfem::Array<int> vertices;
            src.getHandle().GetBdrElementVertices(i - 1, vertices);
            res->edge[i].a = vertices[0] + 1;
            res->edge[i].b = vertices[1] + 1;
            res->edge[i].ref = src.getHandle().GetBdrAttribute(i - 1);
            res->edge[i].tag |= MG_REF + MG_BDY;
         }

         // Copy triangles
         int reorientedCount = 0;
         for (int i = 1; i <= res->nt; i++)
         {
            MMG5_pTria pt = &res->tria[i];
            mfem::Array<int> vertices;
            src.getHandle().GetElementVertices(i - 1, vertices);
            for (int j = 0; j < 3; j++)
              pt->v[j] = vertices[j] + 1;

            pt->ref = src.getHandle().GetAttribute(i - 1);
            for (int j = 0; j < 3; j++)
            {
               res->point[pt->v[j]].tag &= ~MG_NUL;
               pt->edg[j] = 0;
            }

            // Check orientation
            double orientation = MMG2D_quickarea(
                  res->point[pt->v[0]].c,
                  res->point[pt->v[1]].c,
                  res->point[pt->v[2]].c);

            if(orientation < 0)
            {
              int tmp = pt->v[2];
              pt->v[2] = pt->v[1];
              pt->v[1] = tmp;
              reorientedCount++;
              Alert::Warning()
                << "Element reoriented: " << i << Alert::Raise;
            }
          }
          if (reorientedCount > 0)
          {
            Alert::Warning()
              << "Number of elements reoriented: " << std::to_string(reorientedCount)
              << Alert::Raise;
          }
      }
      else if (src.getDimension() == 3) // Use MMG3D
      {
         res->ne = src.getHandle().GetNE();
         res->nt = src.getHandle().GetNBE();

         MMG3D_Set_commonFunc();
         if (!MMG3D_zaldy(res))
            Alert::Exception("Memory allocation for MMG5_pMesh failed.").raise();

         // Copy points
         for (int i = 1; i <= res->np; i++)
         {
           const double* coords = src.getHandle().GetVertex(i - 1);
           std::copy(coords, coords + 3, res->point[i].c);
         }

         // Copy triangles
         for (int i = 1; i <= res->nt; i++)
         {
           MMG5_pTria pt = &res->tria[i];
           mfem::Array<int> vertices;
           src.getHandle().GetBdrElementVertices(i - 1, vertices);
           for (int j = 0; j < 3; j++)
             pt->v[j] = vertices[j] + 1;
           pt->ref = src.getHandle().GetBdrAttribute(i - 1);
         }

         // Copy tetrahedra
         int reorientedCount = 0;
         for (int i = 1; i <= res->ne; i++)
         {
            MMG5_pTetra pt = &res->tetra[i];
            mfem::Array<int> vertices;
            src.getHandle().GetElementVertices(i - 1, vertices);
            for (int j = 0; j < 4; j++)
              pt->v[j] = vertices[j] + 1;

            pt->ref = src.getHandle().GetAttribute(i - 1);

            for (int j = 0; j < 4; j++)
              res->point[pt->v[j]].tag &= ~MG_NUL;

            if (MMG5_orvol(res->point, pt->v) < 0.0)
            {
              reorientedCount++;
              int aux = pt->v[2];
              pt->v[2] = pt->v[3];
              pt->v[3] = aux;
              Alert::Warning() << "Element reoriented: " << i << Alert::Raise;
            }
            if (reorientedCount > 0)
            {
              Alert::Warning()
                << "Number of elements reoriented: " << std::to_string(reorientedCount)
                << Alert::Raise;
            }
         }
      }
      else
      {
         Alert::Exception("Unhandled case.").raise();
         return nullptr;
      }
      return res;
   }

   Rodin::Mesh<Traits::Serial> MMG5::meshToRodin(const MMG5_pMesh src)
   {
      bool isSurface = isSurfaceMesh(src);
      mfem::Mesh dst;

      switch (src->dim) // Switch over space dimension
      {
         case 2:
         {
            assert(!isSurface); // Surface embedded in 2D space not handled yet
            dst = mfem::Mesh(
                  src->dim, // Dimension
                  src->np,  // Vertex count
                  src->nt,  // Element count,
                  src->na   // Boundary element count
                  );
            for (int i = 1; i <= src->np; i++)
            {
               dst.AddVertex(
                     src->point[i].c[0],
                     src->point[i].c[1],
                     src->point[i].c[2]);
            }

            for (int i = 1; i <= src->nt; i++)
            {
               dst.AddTriangle(
                     src->tria[i].v[0] - 1,
                     src->tria[i].v[1] - 1,
                     src->tria[i].v[2] - 1,
                     src->tria[i].ref);
            }

            for (int i = 1; i <= src->na; i++)
            {
               dst.AddBdrSegment(
                     src->edge[i].a - 1,
                     src->edge[i].b - 1,
                     src->edge[i].ref);
            }
            break;
         }
         case 3:
         {
            if (isSurface)
            {
               dst = mfem::Mesh(
                     2,       // Dimension
                     src->np, // Vertex count
                     src->nt, // Element count
                     src->na, // Boundary element count
                     src->dim // Space dim
                     );
               for (int i = 1; i <= src->np; i++)
               {
                  dst.AddVertex(
                       src->point[i].c[0],
                       src->point[i].c[1],
                       src->point[i].c[2]);
               }

               for (int i = 1; i <= src->nt; i++)
               {
                  dst.AddTriangle(
                        src->tria[i].v[0] - 1,
                        src->tria[i].v[1] - 1,
                        src->tria[i].v[2] - 1,
                        src->tria[i].ref);
               }

               for (int i = 1; i <= src->na; i++)
               {
                  dst.AddBdrSegment(
                        src->edge[i].a - 1,
                        src->edge[i].b - 1,
                        src->edge[i].ref);
               }
            }
            else
            {
               dst = mfem::Mesh(
                     2,       // Dimension
                     src->np, // Vertex count
                     src->ne, // Element count
                     src->nt  // Boundary element count
                     );
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
                  dst.AddBdrTriangle(
                        src->tria[i].v[0] - 1,
                        src->tria[i].v[1] - 1,
                        src->tria[i].v[2] - 1,
                        src->tria[i].ref);
               }

               for (int i = 1; i <= src->ne; i++)
               {
                  dst.AddTet(
                        src->tetra[i].v[0] - 1,
                        src->tetra[i].v[1] - 1,
                        src->tetra[i].v[2] - 1,
                        src->tetra[i].v[3] - 1,
                        src->tetra[i].ref);
               }
            }
            break;
         }
         default:
         {
            Alert::Exception("Unhandled case").raise();
         }
      }

      dst.FinalizeMesh(0, true);
      return Rodin::Mesh<Traits::Serial>(std::move(dst));
   }

   void MMG5::copySolution(const MMG5_pSol src, MMG5_pSol dst)
   {
      // Copy all non-pointer fields
      *dst = *src;

      if (dst->np)
      {
         // So (dst->size + 1) * (dst->np + 1) seems to work for most
         // applications
         MMG5_SAFE_CALLOC(dst->m, (dst->size + 1) * (dst->np + 1), double,
               Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
         std::copy(src->m, src->m + dst->size * (dst->np + 1), dst->m);
      }
      else
      {
         dst->m = nullptr;
      }

      if (src->namein)
      {
         auto nameInLength = std::strlen(src->namein);
         MMG5_SAFE_CALLOC(dst->namein, nameInLength, char,
               Alert::Exception("Failed to allocate memory for the MMG5_pSol->namein").raise());
         std::memcpy(dst->namein, src->namein, nameInLength + 1);
      }
      else
      {
         dst->namein = nullptr;
      }

      if (src->nameout)
      {
         auto nameOutLength = std::strlen(src->nameout);
         MMG5_SAFE_CALLOC(dst->nameout, nameOutLength, char,
               Alert::Exception("Failed to allocate memory for the MMG5_pSol->nameout").raise());
         std::memcpy(dst->nameout, src->nameout, nameOutLength + 1);
      }
      else
      {
         dst->nameout = nullptr;
      }
   }

   void MMG5::swapSolution(MMG5_pSol a, MMG5_pSol b)
   {
      std::swap(a->dim, b->dim);
      std::swap(a->entities, b->entities);
      std::swap(a->m, b->m);
      std::swap(a->namein, b->namein);
      std::swap(a->nameout, b->nameout);
      std::swap(a->np, b->np);
      std::swap(a->npi, b->npi);
      std::swap(a->npmax, b->npmax);
      std::swap(a->size, b->size);
      std::swap(a->type, b->type);
      std::swap(a->umax, b->umax);
      std::swap(a->umin, b->umin);
      std::swap(a->ver, b->ver);
   }

   void MMG5::destroySolution(MMG5_pSol sol)
   {
      if (sol)
      {
         if (sol->m)
            MMG5_SAFE_FREE(sol->m);
         if (sol->namein)
            MMG5_SAFE_FREE(sol->namein);
         if (sol->nameout)
            MMG5_SAFE_FREE(sol->nameout);
         MMG5_SAFE_FREE(sol);
      }
      sol = nullptr;
   }

   void MMG5::destroyMesh(MMG5_pMesh mesh)
   {
      if (isSurfaceMesh(mesh))
      {
         MMGS_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_end);
      }
      else if (mesh->dim == 2)
      {
         MMG2D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_end);
      }
      else if (mesh->dim == 3)
      {
         MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_end);
      }
      else
      {
         assert(false);
      }
   }

   MMG5& MMG5::setParameters(MMG5_pMesh mesh)
   {
      bool isSurface = isSurfaceMesh(mesh);
      switch(mesh->dim)
      {
         case 2:
         {
            assert(!isSurface);
            if (m_hmin)
            {
               MMG2D_Set_dparameter(
                  mesh, nullptr, MMG2D_DPARAM_hmin, *m_hmin);
            }
            if (m_hmax)
            {
               MMG2D_Set_dparameter(
                  mesh, nullptr, MMG2D_DPARAM_hmax, *m_hmax);
            }
            if (m_hgrad)
            {
               MMG2D_Set_dparameter(
                  mesh, nullptr, MMG2D_DPARAM_hgrad, *m_hgrad);
            }
            if (m_hausd)
            {
               MMG2D_Set_dparameter(
                  mesh, nullptr, MMG2D_DPARAM_hausd, *m_hausd);
            }
            break;
         }
         case 3:
         {
            if (isSurface)
            {
               if (m_hmin)
               {
                  MMGS_Set_dparameter(
                     mesh, nullptr, MMGS_DPARAM_hmin, *m_hmin);
               }
               if (m_hmax)
               {
                  MMGS_Set_dparameter(
                     mesh, nullptr, MMGS_DPARAM_hmax, *m_hmax);
               }
               if (m_hgrad)
               {
                  MMGS_Set_dparameter(
                     mesh, nullptr, MMGS_DPARAM_hgrad, *m_hgrad);
               }
               if (m_hausd)
               {
                  MMGS_Set_dparameter(
                     mesh, nullptr, MMGS_DPARAM_hausd, *m_hausd);
               }
            }
            else
            {
               if (m_hmin)
               {
                  MMG3D_Set_dparameter(
                     mesh, nullptr, MMG3D_DPARAM_hmin, *m_hmin);
               }
               if (m_hmax)
               {
                  MMG3D_Set_dparameter(
                     mesh, nullptr, MMG3D_DPARAM_hmax, *m_hmax);
               }
               if (m_hgrad)
               {
                  MMG3D_Set_dparameter(
                     mesh, nullptr, MMG3D_DPARAM_hgrad, *m_hgrad);
               }
               if (m_hausd)
               {
                  MMG3D_Set_dparameter(
                     mesh, nullptr, MMG3D_DPARAM_hausd, *m_hausd);
               }
            }
            break;
            default:
               Alert::Exception("Unhandled case").raise();
         }
      }
      return *this;
   }
}

