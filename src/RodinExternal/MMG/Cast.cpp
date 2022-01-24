/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

#include <mmg2d/mmg2d.h>

/*
 * We have to undef the I macro (from complex.h) since it clashes with mfem
 * code (i.e. table.hpp) where the I variable is defined and causes all sorts
 * of nasty errors.
 */
#ifdef I
#undef I
#endif

#include "Rodin/Alert.h"
#include "Rodin/Variational.h"

#include "ScalarSolution2D.h"

#include "Cast.h"


namespace Rodin
{
  using namespace External;

  // ---- <MMG::Mesh2D, Rodin::Mesh> -----------------------------------------
  template <>
  template <>
  Rodin::Mesh
  Cast<MMG::Mesh2D>::to<Rodin::Mesh>()
  const
  {
    auto& mesh = from();
    mfem::Mesh res(
         2,
         mesh.count<MMG::Mesh2D::Entity::Vertex>(),
         mesh.count<MMG::Mesh2D::Entity::Triangle>());

    auto& mmgMesh = mesh.getHandle();

    if (mmgMesh->nt == 0)
       Alert::Exception("Converting from a non-triangular MMG::Mesh2D to a"
             " Rodin::Mesh is not supported ").raise();

    bool shiftAttr = false;
    for (int i = 1; i <= mmgMesh->nt; i++)
    {
      if (mmgMesh->tria[i].ref == 0)
        shiftAttr = true;
      if (mmgMesh->tria[i].ref < 0)
        Alert::Exception(
            "Negative element attributes are not supported").raise();
    }

    bool shiftBdrAttr = false;
    for (int i = 1; i <= mmgMesh->na; i++)
    {
      if (mmgMesh->edge[i].ref == 0)
        shiftBdrAttr = true;
      if (mmgMesh->edge[i].ref < 0)
        Alert::Exception(
            "Negative boundary element attributes are not supported.").raise();
    }

    if (shiftAttr)
      Alert::Warning(
          "Elements with attribute equal to 0 are not supported. "
          "All element attributes will be incremented by 1.").raise();

    if (shiftBdrAttr)
      Alert::Warning(
          "Boundary elements with attribute equal to 0 are not supported. "
          "All boundary element attributes will be incremented by 1.").raise();

    /* So for some reason mmg types are 1 indexed. So when accessing the
     * arrays make sure to start at 1 and not 0. I don't know why this is the
     * case and I'm not sure if it's for every array in the library.
     */
    for (int i = 1; i <= mesh.getHandle()->np; i++)
    {
      res.AddVertex(
          mmgMesh->point[i].c[0],
          mmgMesh->point[i].c[1]
          );
    }

    for (int i = 1; i <= mmgMesh->nt; i++)
    {
      res.AddTriangle(
          mmgMesh->tria[i].v[0] - 1,
          mmgMesh->tria[i].v[1] - 1,
          mmgMesh->tria[i].v[2] - 1,
          mmgMesh->tria[i].ref + shiftAttr
          );
    }
    for (int i = 1; i <= mmgMesh->na; i++)
    {
       res.AddBdrSegment(
             mmgMesh->edge[i].a - 1,
             mmgMesh->edge[i].b - 1,
             mmgMesh->edge[i].ref + shiftBdrAttr
             );
    }
    res.FinalizeMesh(0, true);

    return Rodin::Mesh(res);
  }

  // ---- From: Rodin::Mesh - To: MMG::Mesh2D> -------------------------------
  template <>
  template <>
  MMG::Mesh2D
  Cast<Rodin::Mesh>::to<MMG::Mesh2D>()
  const
  {
    auto& mesh = from();
    auto& mfemMesh = mesh.getHandle();

    if (mesh.getDimension() != 2)
      Alert::Exception("Mesh must be two dimensional.").raise();

    if (mfemMesh.GetNE() == 0)
      Alert::Exception("Converting from an empty mesh is not supported.").raise();

    if (mfemMesh.NURBSext)
       Alert::Exception(
             "Converting from a NURBS mfem::Mesh to MMG::Mesh2D is"
             " not supported.").raise();

    mfem::Array<mfem::Geometry::Type> geoms;
    mfemMesh.GetGeometries(2, geoms);
    if (std::any_of(geoms.begin(), geoms.end(),
             [](mfem::Geometry::Type geometry)
             {
               return geometry != mfem::Geometry::Type::TRIANGLE;
             }))
    {
      Alert::Exception(
            "Converting from a non-triangular mfem::Mesh to MMG::Mesh2D is not"
            " not supported.").raise();
    }

    /*
     * To build the MMG mesh we follow the same procedure as that of the
     * function MMG2D_loadMesh in inout_2d.c
     */
    MMG::Mesh2D res;
    auto& mmgMesh = res.getHandle();

    mmgMesh->np = mmgMesh->nt = mmgMesh->na = mmgMesh->xp = 0;
    mmgMesh->np = mfemMesh.GetNV();
    mmgMesh->nt = mfemMesh.GetNE();
    mmgMesh->na = mfemMesh.GetNBE();

    mmgMesh->npi  = mmgMesh->np;
    mmgMesh->nai  = mmgMesh->na;
    mmgMesh->nti  = mmgMesh->nt;

    MMG2D_Set_commonFunc();
    if (!MMG2D_zaldy(mmgMesh))
       Alert::Exception("Memory allocation for the MMG mesh failed.").raise();

    // Copy points
    for (int i = 1; i <= mmgMesh->np; i++)
    {
      const double* coords = mfemMesh.GetVertex(i - 1);
      std::copy(coords, coords + 3, mmgMesh->point[i].c);
    }

    // Copy edges
    for (int i = 1; i <= mmgMesh->na; i++)
    {
       mfem::Array<int> vertices;
       mfemMesh.GetBdrElementVertices(i - 1, vertices);
       mmgMesh->edge[i].a = vertices[0] + 1;
       mmgMesh->edge[i].b = vertices[1] + 1;
       mmgMesh->edge[i].ref = mfemMesh.GetBdrAttribute(i - 1);
       mmgMesh->edge[i].tag |= MG_REF + MG_BDY;
    }

    // Copy triangles
    int reorientedCount = 0;
    for (int i = 1; i <= mmgMesh->nt; i++)
    {
      MMG5_pTria pt = &mmgMesh->tria[i];
      mfem::Array<int> vertices;
      mfemMesh.GetElementVertices(i - 1, vertices);
      for (int j = 0; j < 3; j++)
        pt->v[j] = vertices[j] + 1;

      pt->ref = mfemMesh.GetAttribute(i - 1);
      for (int j = 0; j < 3; j++)
      {
         mmgMesh->point[pt->v[j]].tag &= ~MG_NUL;
         pt->edg[j] = 0;
      }

      // Check orientation
      double orientation = MMG2D_quickarea(
            mmgMesh->point[pt->v[0]].c,
            mmgMesh->point[pt->v[1]].c,
            mmgMesh->point[pt->v[2]].c);

      if(orientation < 0)
      {
        int tmp = pt->v[2];
        pt->v[2] = pt->v[1];
        pt->v[1] = tmp;
        reorientedCount++;
        Alert::Warning()
          << "Bad orientation in element " << i << ". "
          << "Number of elements reoriented: " << std::to_string(reorientedCount)
          << Alert::Raise;
      }
    }
    return res;
  }

  // ---- From: Rodin::GridFunction - To: MMG::ScalarSolution2D<false>> ------
  template <>
  template <>
  External::MMG::IncompleteScalarSolution2D
  Cast<Variational::IncompleteGridFunction>
  ::to<External::MMG::IncompleteScalarSolution2D>()
  const
  {
    auto& gf = from();
    int size = gf.getHandle().Size();
    const double* data = gf.getHandle().GetData();

    if (!size)
      return External::MMG::IncompleteScalarSolution2D();
    else
    {
      External::MMG::IncompleteScalarSolution2D res(size);
      // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m + size + 1.
      std::copy(data, data + size, res.getHandle()->m + 1);
      return res;
    }
  }

   template<>
   template <>
   Variational::IncompleteGridFunction
   Cast<External::MMG::ScalarSolution2D>
   ::to<Variational::IncompleteGridFunction>() const
   {
     auto& sol = from();
     MMG5_pSol mmgSol = sol.getHandle();
     assert(mmgSol->type == MMG5_Scalar);
     Variational::IncompleteGridFunction res;
     double* data = new double[mmgSol->np];
     // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
     // + np + 1.
     std::copy(mmgSol->m + 1, mmgSol->m + mmgSol->np + 1, data);
     res.getHandle().SetDataAndSize(data, mmgSol->np);
     res.getHandle().MakeDataOwner();
     return res;
   }

   template <>
   template <>
   Variational::IncompleteGridFunction
   Cast<External::MMG::IncompleteScalarSolution2D>
   ::to<Variational::IncompleteGridFunction>() const
   {
     auto& sol = from();
     MMG5_pSol mmgSol = sol.getHandle();
     assert(mmgSol->type == MMG5_Scalar);
     Variational::IncompleteGridFunction res;
     double* data = new double[mmgSol->np];
     // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at m
     // + np + 1.
     std::copy(mmgSol->m + 1, mmgSol->m + mmgSol->np + 1, data);
     res.getHandle().SetDataAndSize(data, mmgSol->np);
     res.getHandle().MakeDataOwner();
     return res;
   }

   template <>
   template <>
   External::MMG::IncompleteScalarSolution2D
   Cast<Variational::GridFunction<Variational::H1>>
   ::to<External::MMG::IncompleteScalarSolution2D>()
   const
   {
     auto& gf = from();
     auto [data, size] = gf.getData();
     if (!size)
       return External::MMG::IncompleteScalarSolution2D();
     else
     {
       External::MMG::IncompleteScalarSolution2D res(size);
       // MMG5_pSol->m is 1 indexed. We must start at m + 1 and finish at
       // m + size + 1.
       std::copy(data, data + size, res.getHandle()->m + 1);
       return res;
     }
   }
}
