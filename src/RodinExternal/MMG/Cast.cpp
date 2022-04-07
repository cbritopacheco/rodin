/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Alert.h"

#include "Rodin/Mesh.h"
#include "Rodin/Variational.h"

#include "Common.h"
#include "Utility.h"

#include "Mesh2D.h"
#include "ScalarSolution2D.h"
#include "VectorSolution2D.h"

#include "Mesh3D.h"
#include "ScalarSolution3D.h"

#include "MeshS.h"
#include "ScalarSolutionS.h"
#include "VectorSolutionS.h"

#include "Cast.h"

namespace Rodin
{
  using namespace External;

   // ---- mmg2d -------------------------------------------------------------
  template <>
  template <>
  Rodin::Mesh<>
  Cast<MMG::Mesh2D>::to<Rodin::Mesh<>>()
  const
  {
    auto& mesh = from();
    auto& src = mesh.getHandle();

    mfem::Mesh dst(
         2,
         mesh.count(MMG::Mesh2D::Entity::Vertex),
         mesh.count(MMG::Mesh2D::Entity::Triangle));

    if (src->nt == 0)
    {
       Alert::Exception()
         << "MMG::Mesh2D is empty (triangle count equals zero)"
         << Alert::Raise;
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

    bool shiftTriAttr = false;
    for (int i = 1; i <= src->nt; i++)
    {
      if (src->tria[i].ref == 0)
        shiftTriAttr = true;
      if (src->tria[i].ref < 0)
        Alert::Exception(
            "Negative triangle attributes are not supported").raise();
    }

    if (shiftEdgeAttr)
      Alert::Warning(
          "Edges with attribute equal to 0 are not supported. "
          "All edge attributes will be incremented by 1.").raise();

    if (shiftTriAttr)
      Alert::Warning(
          "Triangles with attribute equal to 0 are not supported. "
          "All triangle attributes will be incremented by 1.").raise();

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

    return Rodin::Mesh<>(std::move(dst));
  }

  template <>
  template <>
  MMG::Mesh2D
  Cast<Rodin::Mesh<>>::to<MMG::Mesh2D>()
  const
  {
    auto& mesh = from();
    auto& mfemMesh = mesh.getHandle();

    if (mesh.getDimension() != 2)
      Alert::Exception("Mesh must be two dimensional.").raise();

    if (mesh.getSpaceDimension() != 2)
      Alert::Exception("Mesh must be embedded in two dimensional space.").raise();

    if (mfemMesh.GetNE() == 0)
      Alert::Exception("Converting from an empty mesh is not supported.").raise();

    if (mfemMesh.NURBSext)
       Alert::Exception(
             "Converting from a NURBS mfem::Mesh to an MMG::Mesh2D is not supported.").raise();

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
       Alert::Exception("Memory allocation for MMG::Mesh2D mesh failed.").raise();

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
          << "Element reoriented: " << i << Alert::Raise;
      }
    }
    if (reorientedCount > 0)
    {
      Alert::Warning()
        << "Number of elements reoriented: " << std::to_string(reorientedCount)
        << Alert::Raise;
    }
    return res;
  }

  template <>
  template <>
  Variational::IncompleteGridFunction
  Cast<External::MMG::ScalarSolution2D>
  ::to<Variational::IncompleteGridFunction>() const
  {
    Variational::IncompleteGridFunction res;
    External::MMG::MMG5_Sol_To_Rodin_IncompleteGridFunction(from().getHandle(), res);
    return res;
  }

  template <>
  template <>
  External::MMG::IncompleteScalarSolution2D
  Cast<Variational::GridFunction<Variational::H1>>
  ::to<External::MMG::IncompleteScalarSolution2D>()
  const
  {
    assert(from().getFiniteElementSpace().getVectorDimension() == 1);
    External::MMG::IncompleteScalarSolution2D res;
    External::MMG::Rodin_GridFunction_To_MMG5_Sol(from(), res.getHandle());
    return res;
  }

  template <>
  template <>
  External::MMG::IncompleteVectorSolution2D
  Cast<Variational::GridFunction<Variational::H1>>
  ::to<External::MMG::IncompleteVectorSolution2D>() const
  {
    assert(from().getFiniteElementSpace().getVectorDimension() == 2);
    External::MMG::IncompleteVectorSolution2D res;
    External::MMG::Rodin_GridFunction_To_MMG5_Sol(from(), res.getHandle());
    return res;
  }

  // ---- mmg3d -------------------------------------------------------------
  template <>
  template <>
  Rodin::Mesh<>
  Cast<MMG::Mesh3D>::to<Rodin::Mesh<>>() const
  {
     auto& mesh = from();
     mfem::Mesh dst(
          3,
          mesh.count(MMG::Mesh3D::Entity::Vertex),
          mesh.count(MMG::Mesh3D::Entity::Tetrahedra));

     auto& src = mesh.getHandle();
     if (src->ne == 0)
     {
        Alert::Exception()
          << "MMG::Mesh3D is empty (tetrahedron count equals zero)"
          << Alert::Raise;
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

     bool shiftTriAttr = false;
     for (int i = 1; i <= src->nt; i++)
     {
       if (src->tria[i].ref == 0)
         shiftTriAttr = true;
       if (src->tria[i].ref < 0)
         Alert::Exception(
             "Negative triangle attributes are not supported").raise();
     }

     bool shiftTetAttr = false;
     for (int i = 1; i <= src->ne; i++)
     {
       if (src->tetra[i].ref == 0)
         shiftTetAttr = true;
       if (src->tetra[i].ref < 0)
         Alert::Exception(
             "Negative tetrahedron attributes are not supported").raise();
     }

     if (shiftEdgeAttr)
       Alert::Warning(
           "Edges with attribute equal to 0 are not supported. "
           "All edge attributes will be incremented by 1.").raise();

     if (shiftTriAttr)
       Alert::Warning(
           "Triangles with attribute equal to 0 are not supported. "
           "All triangle attributes will be incremented by 1.").raise();

     if (shiftTetAttr)
       Alert::Warning(
           "Tetrahedron with attribute equal to 0 are not supported. "
           "All tetrahedron attributes will be incremented by 1.").raise();

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
       dst.AddBdrTriangle(
           src->tria[i].v[0] - 1,
           src->tria[i].v[1] - 1,
           src->tria[i].v[2] - 1,
           src->tria[i].ref + shiftTriAttr
           );
     }

     for (int i = 1; i <= src->ne; i++)
     {
        dst.AddTet(
              src->tetra[i].v[0] - 1,
              src->tetra[i].v[1] - 1,
              src->tetra[i].v[2] - 1,
              src->tetra[i].v[3] - 1,
              src->tetra[i].ref + shiftTetAttr);
     }

     dst.FinalizeMesh(0, true);

     return Rodin::Mesh<>(std::move(dst));
  }

  template <>
  template <>
  External::MMG::Mesh3D
  Cast<Rodin::Mesh<>>::to<External::MMG::Mesh3D>() const
  {
    auto& mesh = from();
    auto& mfemMesh = mesh.getHandle();

    if (mesh.getDimension() != 3)
      Alert::Exception("Mesh must be three dimensional.").raise();

    if (mesh.getSpaceDimension() != 3)
      Alert::Exception("Mesh must be embedded in three dimensional space.").raise();

    if (mfemMesh.GetNE() == 0)
      Alert::Exception("Converting from an empty mesh is not supported.").raise();

    if (mfemMesh.NURBSext)
       Alert::Exception(
             "Converting from a NURBS mfem::Mesh to an MMG::Mesh3D is not supported.").raise();

    mfem::Array<mfem::Geometry::Type> geoms;
    mfemMesh.GetGeometries(3, geoms);
    if (std::any_of(geoms.begin(), geoms.end(),
             [](mfem::Geometry::Type geometry)
             {
               return geometry != mfem::Geometry::Type::TETRAHEDRON;
             }))
    {
      Alert::Exception()
        << "Converting from a non-tetrahedral mfem::Mesh to an MMG::Mesh3D is not"
        << " not supported."
        << Alert::Raise;
    }

    /*
     * To build the MMG mesh we follow the same procedure as that of the
     * function MMG2D_loadMesh in inout_2d.c
     */
    MMG::Mesh3D res;
    auto& mmgMesh = res.getHandle();

    mmgMesh->np = mmgMesh->nt = mmgMesh->na = mmgMesh->xp = 0;
    mmgMesh->np = mfemMesh.GetNV();
    mmgMesh->ne = mfemMesh.GetNE();
    mmgMesh->nt = mfemMesh.GetNBE();

    mmgMesh->npi  = mmgMesh->np;
    mmgMesh->nai  = mmgMesh->na;
    mmgMesh->nti  = mmgMesh->nt;

    MMG3D_Set_commonFunc();
    if (!MMG3D_zaldy(mmgMesh))
       Alert::Exception("Memory allocation for MMG::Mesh3D mesh failed.").raise();

    // Copy points
    for (int i = 1; i <= mmgMesh->np; i++)
    {
      const double* coords = mfemMesh.GetVertex(i - 1);
      std::copy(coords, coords + 3, mmgMesh->point[i].c);
    }

    // Copy triangles
    for (int i = 1; i <= mmgMesh->nt; i++)
    {
      MMG5_pTria pt = &mmgMesh->tria[i];
      mfem::Array<int> vertices;
      mfemMesh.GetBdrElementVertices(i - 1, vertices);
      for (int j = 0; j < 3; j++)
        pt->v[j] = vertices[j] + 1;
      pt->ref = mfemMesh.GetBdrAttribute(i - 1);
    }

    // Copy tetrahedra
    int reorientedCount = 0;
    for (int i = 1; i <= mmgMesh->ne; i++)
    {
      MMG5_pTetra pt = &mmgMesh->tetra[i];
      mfem::Array<int> vertices;
      mfemMesh.GetElementVertices(i - 1, vertices);
      for (int j = 0; j < 4; j++)
        pt->v[j] = vertices[j] + 1;

      pt->ref = mfemMesh.GetAttribute(i - 1);

      for (int j = 0; j < 4; j++)
        mmgMesh->point[pt->v[j]].tag &= ~MG_NUL;

      if (MMG5_orvol(mmgMesh->point, pt->v) < 0.0)
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
    return res;
  }

  template <>
  template <>
  External::MMG::IncompleteScalarSolution3D
  Cast<Variational::GridFunction<Variational::H1>>
  ::to<External::MMG::IncompleteScalarSolution3D>()
  const
  {
    assert(from().getFiniteElementSpace().getVectorDimension() == 1);
    External::MMG::IncompleteScalarSolution3D res;
    External::MMG::Rodin_GridFunction_To_MMG5_Sol(from(), res.getHandle());
    return res;
  }

  // ---- mmgs --------------------------------------------------------------
  template <>
  template <>
  Rodin::Mesh<>
  Cast<MMG::MeshS>::to<Rodin::Mesh<>>() const
  {
     auto& mesh = from();
     mfem::Mesh dst(
          2,
          mesh.count(MMG::MeshS::Entity::Vertex),
          mesh.count(MMG::MeshS::Entity::Triangle),
          0,
          3);

     auto& src = mesh.getHandle();
     if (src->nt == 0)
     {
        Alert::Exception()
          << "MMG::MeshS is empty (triangle count equals zero)"
          << Alert::Raise;
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

     bool shiftTriAttr = false;
     for (int i = 1; i <= src->nt; i++)
     {
       if (src->tria[i].ref == 0)
         shiftTriAttr = true;
       if (src->tria[i].ref < 0)
         Alert::Exception(
             "Negative triangle attributes are not supported").raise();
     }

     if (shiftEdgeAttr)
       Alert::Warning(
           "Edges with attribute equal to 0 are not supported. "
           "All edge attributes will be incremented by 1.").raise();

     if (shiftTriAttr)
       Alert::Warning(
           "Triangles with attribute equal to 0 are not supported. "
           "All triangle attributes will be incremented by 1.").raise();

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

     return Rodin::Mesh<>(std::move(dst));
  }

  template <>
  template <>
  MMG::MeshS
  Cast<Rodin::Mesh<>>::to<MMG::MeshS>()
  const
  {
    auto& mesh = from();
    auto& mfemMesh = mesh.getHandle();

    if (!(mesh.getDimension() == 2 && mesh.getHandle().SpaceDimension() == 3))
    {
      Alert::Exception()
        << "Mesh must have two dimensional elements, "
        << "embedded in three dimensional space"
        << Alert::Raise;
    }

    if (mfemMesh.GetNE() == 0)
      Alert::Exception(
          "Converting from an empty mesh is not supported.").raise();

    if (mfemMesh.NURBSext)
       Alert::Exception(
             "Converting from a NURBS mfem::Mesh to an MMG::MeshS is"
             " not supported.").raise();

    mfem::Array<mfem::Geometry::Type> geoms;
    mfemMesh.GetGeometries(2, geoms);
    if (std::any_of(geoms.begin(), geoms.end(),
             [](mfem::Geometry::Type geometry)
             {
               return geometry != mfem::Geometry::Type::TRIANGLE;
             }))
    {
      Alert::Exception()
        << "Converting from a non-triangular mfem::Mesh to MMG::Mesh2D is not"
        << " not supported."
        << Alert::Raise;
    }

    /*
     * To build the MMG mesh we follow the same procedure as that of the
     * function MMGS_loadMesh in inout_s.c
     */
    MMG::MeshS res;
    auto& mmgMesh = res.getHandle();

    mmgMesh->np = mmgMesh->nt = mmgMesh->na = mmgMesh->xp = 0;
    mmgMesh->np = mfemMesh.GetNV();
    mmgMesh->nt = mfemMesh.GetNE();
    mmgMesh->na = mfemMesh.GetNBE();

    mmgMesh->npi  = mmgMesh->np;
    mmgMesh->nai  = mmgMesh->na;
    mmgMesh->nti  = mmgMesh->nt;

    MMGS_Set_commonFunc();
    if (!MMGS_zaldy(mmgMesh))
       Alert::Exception("Memory allocation for the MMG::SufaceMesh failed.").raise();

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
    }

    // Copy triangles
    for (int i = 1; i <= mmgMesh->nt; i++)
    {
      MMG5_pTria pt = &mmgMesh->tria[i];
      mfem::Array<int> vertices;
      mfemMesh.GetElementVertices(i - 1, vertices);
      for (int j = 0; j < 3; j++)
        pt->v[j] = vertices[j] + 1;

      pt->ref = mfemMesh.GetAttribute(i - 1);
      for (int j = 0; j < 3; j++)
         mmgMesh->point[pt->v[j]].tag &= ~MG_NUL;
    }
    return res;
  }

  template <>
  template <>
  External::MMG::IncompleteScalarSolutionS
  Cast<Variational::GridFunction<Variational::H1>>
  ::to<External::MMG::IncompleteScalarSolutionS>()
  const
  {
    assert(from().getFiniteElementSpace().getVectorDimension() == 1);
    External::MMG::IncompleteScalarSolutionS res;
    External::MMG::Rodin_GridFunction_To_MMG5_Sol(from(), res.getHandle());
    return res;
  }

  template <>
  template <>
  External::MMG::IncompleteVectorSolutionS
  Cast<Variational::GridFunction<Variational::H1>>
  ::to<External::MMG::IncompleteVectorSolutionS>()
  const
  {
    assert(from().getFiniteElementSpace().getVectorDimension() == 3);
    External::MMG::IncompleteVectorSolutionS res;
    External::MMG::Rodin_GridFunction_To_MMG5_Sol(from(), res.getHandle());
    return res;
  }

  template<>
  template <>
  Variational::IncompleteGridFunction
  Cast<External::MMG::ScalarSolutionS>
  ::to<Variational::IncompleteGridFunction>()
  const
  {
    Variational::IncompleteGridFunction res;
    External::MMG::MMG5_Sol_To_Rodin_IncompleteGridFunction(from().getHandle(), res);
    return res;
  }
}
