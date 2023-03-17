/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "SubMesh.h"

namespace Rodin::Geometry
{
  SubMesh<Context::Serial>::Builder::Builder()
  {}

  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::setReference(
      Mesh<Context::Serial>::Builder&& build, SubMesh<Context::Serial>& mesh)
  {
    const size_t dim = mesh.getDimension();
    m_mbuild.emplace(std::move(build));
    m_ref.emplace(std::ref(mesh));
    m_s2ps.resize(dim + 1);
    m_sidx.resize(dim + 1, 0);
    return *this;
  }

  SubMesh<Context::Serial>::Builder&
  SubMesh<Context::Serial>::Builder::include(std::set<Index> indices)
  {
    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    assert(m_mbuild.has_value());
    auto& build = m_mbuild.value();

    assert(m_s2ps.size() == ref.getDimension() + 1);
    const auto& parent = m_ref->get().getParent();

    assert(ref.getDimension() <= parent.getDimension());
    const size_t dim = ref.getDimension();

    /* First we add the simplices that coincide with the submesh dimension (the
     * elements in the submesh) along with the vertices that are incident. Note
     * that the element dimension in the submesh is not necessarily equal to
     * the element dimension in the parent mesh. For example:
     * - If the submesh dimension is 2 and the parent dimension is 2, the
     *   elements of the submesh are elements (e.g. triangles) of the parent mesh.
     * - If the submesh dimension is 2 and the parent dimension is 3, the
     *   elements of the submesh are faces (e.g. triangles) of the elements
     *   (e.g. tetrahedra) of the parent mesh.
     * - If the submesh dimension is 1 and the parent dimension is 3, the
     *   elements of the submesh are the facets (e.g. edges) of the faces (e.g.
     *   triangles) of the elements (e.g. tetrahedra) of the parent mesh.
     */
    for (const auto& idx : indices)
    {
      auto simplex = parent.getSimplex(dim, idx);

      // Add simplex vertices to the resulting mesh
      const Array<Index>& pvs =
        parent.getConnectivity().getIncidence({dim, 0}, idx);

      Array<Index> vs(pvs.size());
      for (Index i = 0; i < static_cast<Index>(vs.size()); i++)
      {
        Index pvidx = pvs[i];
        if (m_s2ps[0].right.count(pvidx) == 0) // Only add vertex if it is not in the map
        {
          build.vertex(parent.getVertex(pvidx)->getCoordinates());
          vs[i] = m_sidx[0]++;
          m_s2ps[0].insert({ vs[i], pvidx });
        }
        else // Else get the id of the vertex in the submesh
        {
          vs[i] = m_s2ps[0].right.at(pvidx);
        }
      }

      // Add element with the new vertex ordering
      build.element(simplex->getGeometry(), vs);
      m_s2ps[dim].insert({ m_sidx[dim]++, idx });
    }

    /* Next we add all the information of the lower dimensional simplices
     * incident to the element. For example:
     *  - If the submesh dimension is 3, we should add the information (e.g.
     *  attribute) of the triangles and edges incident to the tetrahedra (i.e.
     *  the elements) from the parent mesh.
     *  - If the submesh dimension is 2 and the parent dimension is 2, we
     *  should add information about the edges incident to the triangles (i.e.
     *  the elements) from the parent mesh.
     *  - If the submesh dimension is 2 and the parent dimension is 3, we
     *  should add the information about the edges incident to the triangles
     *  (i.e. the faces) from the parent mesh.
     */
    // for (const auto& idx : indices)
    // {
    //   for (size_t d = 1; d <= dim - 1; d++)
    //   {
    //     const Array<Index>& pconnectivity =
    //       parent.getConnectivity().getIncidence({ dim, d }, idx);
    //     for (const auto& i : pconnectivity)
    //     {
    //       if (m_s2ps[d].right.count(i)) // Only add simplex if it is not in the map
    //         continue;
    //       else if (parent.getAttributeIndex().isTracked(d, i))
    //       {
    //         const SimplexIterator simplex = parent.getSimplex(d, i);
    //         build.simplex(simplex->getGeometry(), simplex->getVertices());
    //       }
    //     }
    //   }
    // }

    if (dim == parent.getDimension()) // We are not in the surface case
    {
      for (auto it = parent.getBoundary(); !it.end(); ++it)
      {
        int el1 = -1, el2 = -1;
        parent.getHandle().GetFaceElements(it->getIndex(), &el1, &el2);
        if (indices.count(el1) || indices.count(el2))
        {
          assert(el1 >= 0 || el2 >= 0);
          mfem::Array<int> pvs;
          parent.getHandle().GetFaceVertices(it->getIndex(), pvs);
          Array<Index> vs(pvs.Size());
          for (int i = 0; i < pvs.Size(); i++)
            vs[i] = m_s2ps[0].right.at(pvs[i]);
          build.face(it->getGeometry(), vs);
        }
      }
    }

    return *this;
  }

  void SubMesh<Context::Serial>::Builder::finalize()
  {
    assert(m_ref.has_value());
    auto& ref = m_ref->get();

    assert(m_mbuild.has_value());
    auto& build = m_mbuild.value();

    build.finalize();

    ref.m_s2ps = std::move(m_s2ps);
  }
}
