#include <algorithm>

#include "Rodin/Alert.h"

#include "Utility.h"
#include "Configure.h"

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

  void MMG5_Sol_Copy(const MMG5_pSol src, MMG5_pSol dst)
  {
    // Copy all non-pointer fields
    *dst = *src;

    if (dst->np)
    {
       MMG5_SAFE_CALLOC(dst->m, dst->size * (dst->npmax + 1), double,
             Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
       std::copy(src->m, src->m + dst->size * (dst->npmax + 1), dst->m);
    }

    if (src->namein)
    {
       auto nameInLength = std::strlen(src->namein);
       MMG5_SAFE_CALLOC(dst->namein, nameInLength, char,
             Alert::Exception("Failed to allocate memory for the MMG5_pSol->namein").raise());
       std::memcpy(dst->namein, src->namein, nameInLength + 1);
    }

    if (src->nameout)
    {
       auto nameOutLength = std::strlen(src->nameout);
       MMG5_SAFE_CALLOC(dst->nameout, nameOutLength, char,
             Alert::Exception("Failed to allocate memory for the MMG5_pSol->nameout").raise());
       std::memcpy(dst->nameout, src->nameout, nameOutLength + 1);
    }
  }

  void MMG5_Sol_To_Rodin_IncompleteGridFunction(
        const MMG5_pSol src, Variational::IncompleteGridFunction& dst)
  {
    assert(src->type == MMG5_Scalar);
    double* data = new double[src->np];
    // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
    // + np + 1.
    std::copy(src->m + 1, src->m + src->np + 1, data);
    dst.getHandle().SetDataAndSize(data, src->np);
    dst.getHandle().MakeDataOwner();
  }

  void Load_MMG5_Sol(const boost::filesystem::path& filename, int meshDim, MMG5_pSol sol)
  {
    /*
     * To load the solution file, we use basically the same methodology in
     * MMG2D_loadSol in mmg2d/inout_2d.c
     *
     * We cannot use the actual function since it requires an MMG5_pMesh
     * object which is not used in the call and subcalls of this method,
     * except for keeping track of some maximum memory constraints and
     * propagating the level of verbosity. As such, we have reimplemented the
     * loading utilizing bits of the same code.
     */
    FILE *inm;
    long posnp;
    int iswp, ier, dim, ver, bin, np, nsols;
    int *type;

    // Verbosity is high when in Debug
    ier =  MMG5_loadSolHeader(
          filename.c_str(), meshDim,
          &inm, &ver, &bin, &iswp, &np, &dim, &nsols, &type, &posnp,
          VERBOSITY_LEVEL);

    switch (ier)
    {
       case -1:
          fclose(inm);
          MMG5_SAFE_FREE(type);
          Alert::Exception(
                "Failed to load solution. Invalid data.").raise();
          break;
       case 0:
          fclose(inm);
          if (type)
             MMG5_SAFE_FREE(type);
          Alert::Exception(
                "Failed to load solution. File not found: " + filename.string()).raise();
          break;
       case 1:
          // Success
          break;
       default:
          fclose(inm);
          MMG5_SAFE_FREE(type);
          Alert::Exception(
                "Failed to load solution. Invalid error code returned.").raise();
    }

    if (nsols != 1)
    {
       fclose(inm);
       MMG5_SAFE_FREE(type);
       Alert::Exception(
             "Failed to load solution. Multiple solutions not supported.").raise();
    }

    sol->type = type[0];
    switch (type[0])
    {
       case MMG5_Scalar:
          sol->size = 1;
          break;
       case MMG5_Vector:
          sol->size = 2;
          Alert::Exception(
                "Failed to load solution. Expected sol->size == 1, got 2.").raise();
          break;
       case MMG5_Tensor:
          sol->size = 3;
          Alert::Exception(
                "Failed to load solution. Expected sol->size == 1, got 3.").raise();
          break;
       default:
          Alert::Exception(
                "Failed to load solution. Unknown solution type.").raise();
    }

    assert(!sol->m);
    if (np)
    {
       sol->np  = np;
       sol->npi = np;
       sol->npmax = std::max(static_cast<int>(1.5 * sol->np), MMG2D_NPMAX);
       MMG5_SAFE_CALLOC(
             sol->m, (sol->size * (sol->npmax + 1)), double,
             Alert::Exception("Failed to allocate memory for the MMG5_pSol->m").raise());
    }

    MMG5_SAFE_FREE(type);

    // Read the solutions
    rewind(inm);
    fseek(inm, posnp, SEEK_SET);

    assert(sol->ver > 1);
    assert(!bin);
    for (int k = 1; k <= sol->np; k++)
    {
      double dbuf;
      for (int i = 0; i < sol->size; i++)
      {
         auto read = [&inm, &dbuf] () { MMG_FSCANF(inm, "%lf", &dbuf); return 1; };
         if (read() < 0)
            Alert::Exception("Failed to load solution. Error while reading.").raise();
          sol->m[sol->size * k + i] = dbuf;
      }
    }

    fclose(inm);
  }

  void MMG5_Sol_Swap(MMG5_pSol a, MMG5_pSol b)
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

  void MMG5_Sol_Free(MMG5_pSol sol)
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
}
