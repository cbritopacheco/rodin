#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Utility/Overloaded.h"

#include "LevelSetDiscretizer.h"

namespace Rodin::MMG
{
  LevelSetDiscretizer& LevelSetDiscretizer::surface(bool meshTheSurface)
  {
    m_meshTheSurface = meshTheSurface;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::setLevelSet(Real ls)
  {
    m_ls = ls;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::setRMC(Real rmc)
  {
    m_rmc = rmc;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::setBaseReferences(
    const FlatSet<Geometry::Attribute>& refs)
  {
    m_lsBaseReferences = refs;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::setBoundaryReference(
    const Geometry::Attribute& ref)
  {
    m_isoref = ref;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::split(
    const Geometry::Attribute& ref, const Split& s)
  {
    m_split[ref] = s;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::noSplit(
    const Geometry::Attribute& ref)
  {
    m_split[ref] = NoSplit;
    return *this;
  }

  LevelSetDiscretizer& LevelSetDiscretizer::setSplit(const SplitMap& split)
  {
    assert(split.size() > 0);
    m_split = split;
    return *this;
  }

  ReturnCode LevelSetDiscretizer::discretizeMMG2D(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      MMG2D_Set_dparameter(mesh, sol, MMG2D_DPARAM_rmc, *m_rmc);
    if (getSplitMap().size() > 0)
    {
      assert(getSplitMap().size() == m_uniqueSplit.size());
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_numberOfMat, m_uniqueSplit.size());
      for (const auto& v : m_uniqueSplit)
      {
        const auto& ref = v.first;
        const auto& split = v.second;

        std::visit(
            Utility::Overloaded
            {
              [&](const NoSplitT&)
              {
                if (!MMG2D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_NoSplit, ref, ref))
                {
                  Alert::MemberFunctionException(*this, __func__)
                    << "Could not set the multi-material reference lookup."
                    << Alert::Raise;
                }
              },
              [&](const Split& s)
              {
                if (!MMG2D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_Split,
                      s.interior, s.exterior))
                {
                  Alert::MemberFunctionException(*this, __func__)
                    << "Could not set the multi-material reference lookup."
                    << Alert::Raise;
                }
              }
          }, split);
      }
    }

    if (m_isoref)
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_isoref, *m_isoref);
    if (m_meshTheSurface)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Meshing the surface for a 2D mesh is not supported."
        << Alert::Raise;
    }
    else
    {
      MMG2D_Set_iparameter(mesh, sol, MMG2D_IPARAM_iso, 1);
    }
    MMG2D_Set_dparameter(mesh, sol, MMG2D_DPARAM_ls, m_ls);
    return MMG2D_mmg2dls(mesh, sol, nullptr);
  }

  ReturnCode LevelSetDiscretizer::discretizeMMG3D(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      MMG3D_Set_dparameter(mesh, sol, MMG3D_DPARAM_rmc, *m_rmc);

    if (m_lsBaseReferences.size())
    {
      MMG3D_Set_iparameter(
          mesh, sol, MMG3D_IPARAM_numberOfLSBaseReferences, m_lsBaseReferences.size());
      for (const auto& br : m_lsBaseReferences)
        MMG3D_Set_lsBaseReference(mesh, sol, br);
    }

    if (m_isoref)
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_isoref, *m_isoref);
    if (getSplitMap().size() > 0)
    {
      mesh->memMax *= 2.0; // Real allowed memory because of bug
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_numberOfMat, m_uniqueSplit.size());
      for (const auto& v : m_uniqueSplit)
      {
        const auto& ref = v.first;
        const auto& split = v.second;

        std::visit(
          Utility::Overloaded
          {
            [&](const NoSplitT&)
            {
              if (!MMG3D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_NoSplit, ref, ref))
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Could not set the multi-material reference lookup."
                  << Alert::Raise;
              }
            },
            [&](const Split& s)
            {
              if (!MMG3D_Set_multiMat(
                    mesh, sol, ref, MMG5_MMAT_Split, s.interior, s.exterior))
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Could not set the multi-material reference lookup."
                  << Alert::Raise;
              }
            }
          }, split);
      }
    }

    if (m_meshTheSurface)
    {
      MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_isosurf, 1);
    }
    else
    {
      // MMG3D_Set_iparameter(mesh, sol, MMG3D_IPARAM_iso, 1);
      // We set it manually because otherwise it messes up the mesh references.
      mesh->info.iso = 1;
    }
    MMG3D_Set_dparameter(mesh, sol, MMG3D_DPARAM_ls, m_ls);
    return MMG3D_mmg3dls(mesh, sol, nullptr);
  }

  int LevelSetDiscretizer::discretizeMMGS(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      Alert::Warning() << "Warning RMC option is not supported for surfaces." << Alert::Raise;
    if (m_isoref)
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_isoref, *m_isoref);
    if (getSplitMap().size() > 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Material splitting is not supported for surfaces." << Alert::Raise;
    }

    if (m_meshTheSurface)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Meshing the surface for a surface mesh is not supported."
        << Alert::Raise;
    }
    else
    {
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_iso, 1);
    }
    MMGS_Set_dparameter(mesh, sol, MMGS_DPARAM_ls, m_ls);
    return MMGS_mmgsls(mesh, sol, nullptr);
  }

  void LevelSetDiscretizer::generateUniqueSplit(const FlatSet<Geometry::Attribute>& attrs)
  {
    m_uniqueSplit.clear();

    // Compute existing attributes
    FlatSet<Geometry::Attribute> existingAttributes;
    existingAttributes.insert(attrs.begin(), attrs.end());

    // Add the attributes from the split map to the existing attributes
    for (const auto& [attr, split] : getSplitMap())
    {
      existingAttributes.insert(attr);
      std::visit(
        Utility::Overloaded
        {
          [&](const Split& s)
          {
            existingAttributes.insert(s.interior);
            existingAttributes.insert(s.exterior);
          },
          [&](const NoSplitT&) {}
        }, split);
    }

    // Generate unique splits for each reference
    Geometry::Attribute gen = 1;
    for (const auto it : getSplitMap())
    {
      const auto& attribute = it.first;
      const auto& split = it.second;
      std::visit(
        Utility::Overloaded
        {
          [&](const Split& s)
          {
            // Generate unique interior attribute
            Geometry::Attribute interior;
            do { interior = gen++; } while (existingAttributes.count(interior));
            existingAttributes.insert(interior);
            m_g2om.insert({ interior, s.interior });

            // Generate unique exterior attribute
            Geometry::Attribute exterior;
            do { exterior = gen++; } while (existingAttributes.count(exterior));
            existingAttributes.insert(exterior);
            m_g2om.insert({ exterior, s.exterior });

            // Add the unique split to the map
            m_uniqueSplit.insert({ attribute, Split{ interior, exterior } });
          },
          [&](const NoSplitT&)
          {
            m_uniqueSplit.insert({ attribute, NoSplit });
          }
        }, split);
    }
  }

  void LevelSetDiscretizer::deleteBoundaryRef(MMG5_pMesh mesh, Geometry::Attribute ref)
  {
    if (m_meshTheSurface || MMG5::isSurfaceMesh(mesh) || mesh->dim == 2)
    {
      size_t oldna = mesh->na;
      std::vector<size_t> ids;
      ids.reserve(oldna);
      for (int i = 1; i <= mesh->na; i++)
      {
        if (ref != mesh->edge[i].ref)
          ids.push_back(i);
      }
      size_t newna = ids.size();
      // MMG treats a non-null edge array as holding at least one entry, so an
      // empty result must leave the pointer null. Likewise, namax must stay
      // positive as MMG uses it to size internal allocations.
      MMG5_pEdge edges = nullptr;
      if (newna > 0)
      {
        MMG5_SAFE_CALLOC(edges, newna + 1, MMG5_Edge,
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to reallocate edge memory." << Alert::Raise);
        for (size_t i = 1; i <= newna; i++)
          edges[i] = mesh->edge[ids[i - 1]];
        mesh->namax = newna;
      }
      MMG5_SAFE_FREE(mesh->edge);
      mesh->na = newna;
      mesh->edge = edges;
    }
    else if (mesh->dim == 3)
    {
      size_t oldnt = mesh->nt;
      std::vector<size_t> ids;
      ids.reserve(oldnt);
      for (int i = 1; i <= mesh->nt; i++)
      {
        if (ref != mesh->tria[i].ref)
          ids.push_back(i);
      }
      size_t newnt = ids.size();
      // MMG treats a non-null triangle array as holding at least one entry,
      // so an empty result must leave the pointer null. Likewise, ntmax must
      // stay positive as MMG5_bdrySet sizes xtetra from it.
      MMG5_pTria triangles = nullptr;
      if (newnt > 0)
      {
        MMG5_SAFE_CALLOC(triangles, newnt + 1, MMG5_Tria,
            Alert::MemberFunctionException(*this, __func__)
              << "Failed to reallocate triangles." << Alert::Raise);
        for (size_t i = 1; i <= newnt; i++)
          triangles[i] = mesh->tria[ids[i - 1]];
        mesh->ntmax = newnt;
      }
      MMG5_SAFE_FREE(mesh->tria);
      mesh->nt = newnt;
      mesh->tria = triangles;
    }
    else
    {
      assert(false); // Unhandled case
    }
  }

  MMG::Mesh LevelSetDiscretizer::discretize(const MMG::RealGridFunction& ls)
  {
    const auto& fes = ls.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    IndexSet requiredVertices;
    IndexSet requiredEdges;
    IndexSet requiredTriangles;
    IndexSet requiredTetrahedra;
    if (const auto* inputMesh = dynamic_cast<const MMG::Mesh*>(&mesh))
    {
      requiredVertices = inputMesh->getRequiredVertices();
      requiredEdges = inputMesh->getRequiredEdges();
      requiredTriangles = inputMesh->getRequiredTriangles();
      requiredTetrahedra = inputMesh->getRequiredTetrahedra();
    }

    MMG5_pMesh mmgMesh = nullptr;
    mmgMesh = rodinToMesh(ls.getFiniteElementSpace().getMesh());

    // Erase boundary elements which have the isoref
    if (m_isoref)
      deleteBoundaryRef(mmgMesh, *m_isoref);

    MMG5_pSol sol = createSolution(mmgMesh, ls.getFiniteElementSpace().getVectorDimension());
    copySolution(ls, sol);
    MMG5::setParameters(mmgMesh);
    const bool isSurface = ls.getFiniteElementSpace().getMesh().isSurface();
    const size_t meshDim = mesh.getDimension();

    FlatSet<Geometry::Attribute> attrs;
    if (m_meshTheSurface)
      attrs = mesh.getAttributeIndex().getAttributes(meshDim - 1);
    else
      attrs = mesh.getAttributeIndex().getAttributes(meshDim);

    generateUniqueSplit(attrs);

    ReturnCode retcode = MMG5_STRONGFAILURE;
    switch (mmgMesh->dim)
    {
      case 2:
      {
        assert(!isSurface);
        retcode = discretizeMMG2D(mmgMesh, sol);
        break;
      }
      case 3:
      {
        if (isSurface)
        {
          retcode = discretizeMMGS(mmgMesh, sol);
        }
        else
        {
          retcode = discretizeMMG3D(mmgMesh, sol);
        }
        break;
      }
    }

    if (retcode != MMG5_SUCCESS)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Failed to discretize the implicit domain."
        << Alert::Raise;
    }

    // Delete zero ref
    // deleteBoundaryRef(mmgMesh, 0);

    auto rodinMesh = meshToRodin(mmgMesh);
    rodinMesh.getRequiredVertices().clear();
    const size_t vertexCount = rodinMesh.getVertexCount();
    for (const auto& idx : requiredVertices)
    {
      if (idx < vertexCount)
        rodinMesh.setRequiredVertex(idx);
    }
    rodinMesh.getRequiredEdges().clear();
    const size_t edgeCount =
      rodinMesh.getPolytopeCount(Geometry::Polytope::Type::Segment);
    for (const auto& idx : requiredEdges)
    {
      if (idx < edgeCount)
        rodinMesh.setRequiredEdge(idx);
    }
    rodinMesh.getRequiredTriangles().clear();
    const size_t triangleCount =
      rodinMesh.getPolytopeCount(Geometry::Polytope::Type::Triangle);
    for (const auto& idx : requiredTriangles)
    {
      if (idx < triangleCount)
        rodinMesh.setRequiredTriangle(idx);
    }
    rodinMesh.getRequiredTetrahedra().clear();
    const size_t tetrahedronCount =
      rodinMesh.getPolytopeCount(Geometry::Polytope::Type::Tetrahedron);
    for (const auto& idx : requiredTetrahedra)
    {
      if (idx < tetrahedronCount)
        rodinMesh.setRequiredTetrahedron(idx);
    }
    destroySolution(sol);
    destroyMesh(mmgMesh);

    // Recover original attributes
    if (!isSurface)
    {
      const size_t meshDim = rodinMesh.getDimension();
      size_t d;
      if (m_meshTheSurface)
        d = meshDim - 1;
      else
        d = meshDim;
      for (auto it = rodinMesh.getPolytope(d); !it.end(); ++it)
      {
        const auto& polytope = *it;
        const Index idx = polytope.getIndex();
        const Geometry::Attribute attr = polytope.getAttribute().value_or(0);
        auto attrIt = m_g2om.find(attr);
        if (attrIt != m_g2om.end())
          rodinMesh.setAttribute({ d, idx }, attrIt->second);
        else
          assert(std::holds_alternative<NoSplitT>(m_split.at(attr)));
      }
    }

    return rodinMesh;
  }
}
