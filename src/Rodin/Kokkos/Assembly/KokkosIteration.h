/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_KOKKOSITERATION_H
#define RODIN_KOKKOS_ASSEMBLY_KOKKOSITERATION_H

/**
 * @file
 * @brief Indexed region traversal over a resident Kokkos mesh.
 */

#include "Rodin/Geometry/Region.h"
#include "Rodin/Kokkos/ForwardDecls.h"
#include "Rodin/Kokkos/Geometry/Mesh.h"

namespace Rodin::Assembly
{
  /**
   * @brief Kokkos traversal of the existing mesh assembly regions.
   * @tparam Backend Kokkos execution backend.
   *
   * ## Architecture
   *
   * This descriptor follows the indexed contract of @ref SequentialIteration
   * and @ref OpenMPIteration "OpenMPIteration". It supplies a candidate count,
   * an entity dimension, a region predicate and polytope access. Cells selects
   * all maximal-dimensional entities. Faces, Boundary and Interface scan all
   * codimension-one entities;
   * Boundary accepts exactly one incident cell, and Interface accepts more
   * than one, including nonmanifold adjacency. getCount() is the candidate
   * count before filtering. Attribute selection remains an assembly concern.
   *
   * Construction on the host retains scalar metadata and a private borrowed
   * mesh handle. Capture this descriptor by value in a kernel on the mesh's
   * execution space. filter() and getPolytope() access resident memory and
   * must execute in that memory space; metadata access is also valid on the
   * host. No incidence is computed and no mesh topology is downloaded.
   * Nonempty Boundary and Interface traversal require explicitly materialized
   * face-to-cell incidence. Face regions require positive mesh dimension.
   *
   * The descriptor owns no mesh allocations. Finish kernels before mutating,
   * moving or destroying the mesh, and reconstruct the descriptor afterward.
   * Kernel scheduling and synchronization belong to the assembly backend.
   */
  template <class Backend>
  class KokkosIteration<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>> final
  {
    public:
      /// Resident mesh supplying this traversal.
      using MeshType = Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>;

      /**
       * @brief Borrows a mesh region for indexed device traversal.
       * @param mesh Mesh that must outlive the descriptor and its kernels.
       * @param region Existing Rodin region selecting candidates and filtering.
       */
      KokkosIteration(const MeshType& mesh, const Geometry::Region& region)
        : m_mesh(),
          m_region(region),
          m_dimension(0),
          m_count(0)
      {
        const size_t dimension = mesh.getDimension();
        switch (region)
        {
          case Geometry::Region::Cells:
          {
            m_dimension = dimension;
            m_count = mesh.getCellCount();
            break;
          }
          case Geometry::Region::Faces:
          case Geometry::Region::Boundary:
          case Geometry::Region::Interface:
          {
            if (dimension == 0)
            {
              Alert::MemberFunctionException(*this, __func__)
                << "A zero-dimensional mesh has no faces." << Alert::Raise;
            }
            m_dimension = dimension - 1;
            m_count = mesh.getFaceCount();
            break;
          }
          default:
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Invalid mesh assembly region." << Alert::Raise;
          }
        }
        if (m_count != 0 &&
          (region == Geometry::Region::Boundary ||
            region == Geometry::Region::Interface) &&
          !mesh.getConnectivity().hasIncidence(m_dimension, dimension))
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Face-to-cell incidence has not been computed and is required "
               "to traverse this region."
            << Alert::Raise;
        }
        m_mesh = mesh.getDeviceView();
      }

      /**
       * @brief Gets the dimension of candidate entities.
       * @returns Topological dimension of the selected region.
       */
      KOKKOS_INLINE_FUNCTION
      size_t getDimension() const
      {
        return m_dimension;
      }

      /**
       * @brief Gets the number of candidates before region filtering.
       * @returns All cell candidates or all face candidates.
       */
      KOKKOS_INLINE_FUNCTION
      size_t getCount() const
      {
        return m_count;
      }

      /**
       * @brief Tests region membership in the mesh execution space.
       * @param i Candidate index, less than getCount().
       * @returns Whether the candidate belongs to the selected region.
       */
      KOKKOS_INLINE_FUNCTION
      bool filter(Index i) const
      {
        KOKKOS_ASSERT(i < m_count);
        switch (m_region)
        {
          case Geometry::Region::Cells:
          case Geometry::Region::Faces:
            return true;
          case Geometry::Region::Boundary:
            return m_mesh.getConnectivity()
                     .getIncidence(m_dimension, m_dimension + 1, i)
                     .size() == 1;
          case Geometry::Region::Interface:
            return m_mesh.getConnectivity()
                     .getIncidence(m_dimension, m_dimension + 1, i)
                     .size() > 1;
        }
        KOKKOS_ASSERT(false);
        return false;
      }

      /**
       * @brief Gets a candidate entity in the mesh execution space.
       * @param i Candidate index, less than getCount(), independently of filter().
       * @returns Borrowed device polytope descriptor.
       */
      KOKKOS_INLINE_FUNCTION
      Rodin::Kokkos::Geometry::PolytopeView getPolytope(Index i) const
      {
        KOKKOS_ASSERT(i < m_count);
        return m_mesh.getPolytope(m_dimension, i);
      }

    private:
      Rodin::Kokkos::Geometry::MeshView m_mesh;
      Geometry::Region m_region;
      size_t m_dimension;
      size_t m_count;
  };

  /**
   * @brief Deduces resident traversal from its mesh and region.
   * @tparam Backend Kokkos execution backend.
   * @param mesh Resident mesh to traverse.
   * @param region Existing Rodin assembly region.
   */
  template <class Backend>
  KokkosIteration(const Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>& mesh,
    const Geometry::Region& region)
    -> KokkosIteration<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>;
}

#endif
