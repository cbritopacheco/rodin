#include "ImplicitDomainMesher.h"

namespace Rodin::External::MMG
{
  ImplicitDomainMesher& ImplicitDomainMesher::surface(bool meshTheSurface)
  {
    m_meshTheSurface = meshTheSurface;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setLevelSet(double ls)
  {
    m_ls = ls;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setRMC(double rmc)
  {
    m_rmc = rmc;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setBaseReferences(
    const FlatSet<MaterialAttribute>& refs)
  {
    m_lsBaseReferences = refs;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setBoundaryReference(
    const MaterialAttribute& ref)
  {
    m_isoref = ref;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::split(
    const MaterialAttribute& ref, const Split& s)
  {
    m_split[ref] = s;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::noSplit(
    const MaterialAttribute& ref)
  {
    m_split[ref] = NoSplit;
    return *this;
  }

  ImplicitDomainMesher& ImplicitDomainMesher::setSplit(const SplitMap& split)
  {
    assert(split.size() > 0);
    m_split = split;
    return *this;
  }

  ReturnCode ImplicitDomainMesher::discretizeMMG2D(MMG5_pMesh mesh, MMG5_pSol sol)
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
                  Alert::Exception() << "Could not set the multi-material reference lookup."
                                     << Alert::Raise;
                }
              },
              [&](const Split& s)
              {
                if (!MMG2D_Set_multiMat(mesh, sol, ref, MMG5_MMAT_Split,
                      s.interior, s.exterior))
                {
                  Alert::Exception() << "Could not set the multi-material reference lookup."
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
      Alert::Exception()
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

  ReturnCode ImplicitDomainMesher::discretizeMMG3D(MMG5_pMesh mesh, MMG5_pSol sol)
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
      mesh->memMax *= 2.0; // Double allowed memory because of bug
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
                Alert::Exception() << "Could not set the multi-material reference lookup."
                                   << Alert::Raise;
              }
            },
            [&](const Split& s)
            {
              if (!MMG3D_Set_multiMat(
                    mesh, sol, ref, MMG5_MMAT_Split, s.interior, s.exterior))
              {
                Alert::Exception() << "Could not set the multi-material reference lookup."
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

  int ImplicitDomainMesher::discretizeMMGS(MMG5_pMesh mesh, MMG5_pSol sol)
  {
    if (m_rmc)
      Alert::Warning("Warning RMC option is not supported for surfaces").raise();
    if (m_isoref)
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_isoref, *m_isoref);
    if (getSplitMap().size() > 0)
      Alert::Exception("Material splitting is not supported for surfaces.").raise();

    if (m_meshTheSurface)
    {
      Alert::Exception() << "Meshing the surface for a surface mesh is not supported."
                         << Alert::Raise;
    }
    else
    {
      MMGS_Set_iparameter(mesh, sol, MMGS_IPARAM_iso, 1);
    }
    MMGS_Set_dparameter(mesh, sol, MMGS_DPARAM_ls, m_ls);
    return MMGS_mmgsls(mesh, sol, nullptr);
  }

  void ImplicitDomainMesher::generateUniqueSplit(const FlatSet<Geometry::Attribute>& attr)
  {
    m_uniqueSplit.clear();

    // Compute existing attributes
    FlatSet<Geometry::Attribute> existingAttributes;
    existingAttributes.insert(attr.begin(), attr.end());

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
    MaterialAttribute gen = 1;
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
            MaterialAttribute interior;
            do { interior = gen++; } while (existingAttributes.count(interior));
            existingAttributes.insert(interior);
            m_g2om.insert({ interior, s.interior });

            // Generate unique exterior attribute
            MaterialAttribute exterior;
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

  void ImplicitDomainMesher::deleteBoundaryRef(MMG5_pMesh mesh, MaterialAttribute ref)
  {
    if (m_meshTheSurface || mesh->dim == 2)
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
      MMG5_pEdge edges = nullptr;
      MMG5_SAFE_MALLOC(edges, newna + 1, MMG5_Edge,
          Alert::Exception("Failed to reallocate edges.").raise());
      for (int i = 1; i <= newna; i++)
        edges[i] = mesh->edge[ids[i - 1]];
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
      MMG5_pTria triangles = nullptr;
      MMG5_SAFE_MALLOC(triangles, newnt + 1, MMG5_Tria,
          Alert::Exception("Failed to reallocate triangles.").raise());
      for (int i = 1; i <= newnt; i++)
        triangles[i] = mesh->tria[ids[i - 1]];
      MMG5_SAFE_FREE(mesh->tria);
      mesh->nt = newnt;
      mesh->tria = triangles;
    }
    else
    {
      assert(false); // Unhandled case
    }
  }

  MMG::Mesh ImplicitDomainMesher::discretize(const MMG::ScalarGridFunction& ls)
  {
    const auto& fes = ls.getFiniteElementSpace();
    const auto& mesh = fes.getMesh();
    MMG5_pMesh mmgMesh = nullptr;
    try
    {
      mmgMesh = rodinToMesh(
          dynamic_cast<const MMG::Mesh&>(ls.getFiniteElementSpace().getMesh()));
    } catch (std::bad_cast&)
    {
      Alert::Exception() << "Mesh must be of type MMG::Mesh." << Alert::Raise;
    }
    // Erase boundary elements which have the isoref
    // if (m_isoref)
    //  deleteBoundaryRef(mesh, *m_isoref);

    MMG5_pSol sol = createSolution(mmgMesh, ls.getFiniteElementSpace().getVectorDimension());
    copySolution(ls, sol);
    MMG5::setParameters(mmgMesh);
    const bool isSurface = ls.getFiniteElementSpace().getMesh().isSurface();
    if (m_meshTheSurface)
    {
      assert(false);
      // generateUniqueSplit(ls.getFiniteElementSpace().getMesh().getBoundaryAttributes());
    }
    else
    {
      const size_t meshDim = mesh.getDimension();
      const auto& attributeIndex = mesh.getAttributeIndex();
      FlatSet<Geometry::Attribute> attrs;
      for (auto it = attributeIndex.begin(meshDim); it != attributeIndex.end(meshDim); ++it)
        attrs.insert(it->second);
      generateUniqueSplit(attrs);
    }

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
      Alert::Exception() << "Failed to discretize the implicit domain."
                         << Alert::Raise;
    }

    auto rodinMesh = meshToRodin(mmgMesh);
    destroySolution(sol);
    destroyMesh(mmgMesh);

    // Recover original attributes
    const size_t meshDim = rodinMesh.getDimension();
    for (auto it = rodinMesh.getElement(); !it.end(); ++it)
    {
      const auto& cell = *it;
      const Index idx = cell.getIndex();
      const Geometry::Attribute attr = cell.getAttribute();
      auto attrIt = m_g2om.find(attr);
      if (attrIt != m_g2om.end())
        rodinMesh.setAttribute({ meshDim, idx }, attrIt->second);
      else
        assert(std::holds_alternative<NoSplitT>(m_split.at(attr)));
    }
    return rodinMesh;
  }
}
