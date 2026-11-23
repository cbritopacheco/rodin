/*
 *          Copyright Carlos BRITO PACHECO 2026.
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_ASSEMBLY_POLYTOPES_H
#define RODIN_KOKKOS_ASSEMBLY_POLYTOPES_H

/**
 * @file
 * @brief Prepared entity selection, quadrature and finite-element entry offsets.
 */

#include <memory>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Kokkos/Geometry/Region.h"
#include "Rodin/Kokkos/Variational/Quadrature.h"
#include "Rodin/Revision.h"

namespace Rodin::Assembly
{
  /**
   * @brief Polytopes integrated by a device integrator and the position of
   * their entries.
   *
   * ## Architecture
   *
   * The polytopes are the cells of the mesh or its faces, as selected by the
   * integration region. The quadrature tables and the attribute selection
   * are uploaded on the first integration and again only when the mesh, the
   * region, the quadrature orders or the selection change, so that repeated
   * integrations do not allocate device memory. The entry offsets are
   * located again only when the mesh revision, the region, the selection or
   * the vector dimensions change, so that repeated integrations neither launch
   * the counting kernels nor wait for them.
   *
   * @tparam Mesh Resident Kokkos mesh type.
   */
  template <class Mesh>
  class Polytopes
  {
    public:
      /// @brief Kokkos execution space.
      using ExecutionSpace = typename Mesh::ExecutionSpace;

      /// @brief Kokkos memory space.
      using MemorySpace = typename Mesh::MemorySpace;

      /// @brief Device quadrature tables.
      using QuadratureType = ::Rodin::Kokkos::Variational::Quadrature<Mesh>;

      /// @brief Device attribute selection.
      using AttributesType = ::Kokkos::View<::Rodin::Geometry::Attribute*, MemorySpace>;

      Polytopes() = default;

      /**
       * @brief Starts empty: a copy must not share the offsets it may locate again.
       * @param other Polytopes not copied.
       */
      Polytopes([[maybe_unused]] const Polytopes& other) {}

      /**
       * @brief Move constructor, keeping the located offsets.
       * @param other Polytopes to move.
       */
      Polytopes(Polytopes&& other) = default;

      /**
       * @brief Empties this object, see the copy constructor.
       * @param other Polytopes to assign, only compared to this object.
       * @return This object.
       */
      Polytopes& operator=(const Polytopes& other)
      {
        if (this != &other)
          *this = Polytopes();
        return *this;
      }

      /**
       * @brief Move assignment, keeping the located offsets.
       * @param other Polytopes to move.
       * @return This object.
       */
      Polytopes& operator=(Polytopes&& other) = default;

      /**
       * @brief Locates the entries of every selected polytope.
       * @param mesh Resident mesh.
       * @param region Integration region.
       * @param attributes Selected attributes; empty selects every polytope
       * of the region.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param order Callable giving the quadrature order on a host polytope.
       * @return Number of entries.
       * @throws Alert::Exception If the region is the boundary or the
       * interface and the face to cell incidence is missing.
       */
      template <class TrialDOFs, class TestDOFs, class Order>
      size_t prepare(const Mesh& mesh, ::Rodin::Geometry::Region region,
        const FlatSet<::Rodin::Geometry::Attribute>& attributes,
        const TrialDOFs& trialDOFs, Index trialVDim, const TestDOFs& testDOFs,
        Index testVDim, const Order& order)
      {
        const size_t localDimension = mesh.getDimension();
        const size_t meshDimension = region == ::Rodin::Geometry::Region::Cells
          ? localDimension
          : ::Rodin::Kokkos::Geometry::getDomainDimension(mesh);
        bool checkIncidence = localDimension == meshDimension;
        if constexpr (requires { mesh.getGlobalDimension(); })
          checkIncidence = checkIncidence && mesh.getCellCount() != 0;
        if ((region == ::Rodin::Geometry::Region::Boundary ||
              region == ::Rodin::Geometry::Region::Interface) &&
          checkIncidence)
        {
          if (meshDimension == 0)
          {
            Alert::Exception() << "A zero-dimensional mesh has no facets for integration."
                               << Alert::Raise;
          }
          if (!mesh.getConnectivity().isPrepared(meshDimension - 1, meshDimension))
          {
            Alert::Exception()
              << "Prepare face-to-cell incidence with mesh.getConnectivity().compute("
              << meshDimension - 1 << ", " << meshDimension
              << ") before integrating over the boundary or interface." << Alert::Raise;
          }
        }
        const auto& execution = mesh.getContext().getExecutionSpace();
        const size_t dimension =
          ::Rodin::Kokkos::Geometry::getRegionDimension(meshDimension, region);
        const Index count =
          ::Rodin::Kokkos::Geometry::getRegionPolytopeCount(mesh, region);
        if (m_mesh != &mesh || m_dimension != dimension || m_polytopeCount != count ||
          !m_quadrature || !m_quadrature->hasOrders(mesh, dimension, order))
        {
          if (count == 0)
            m_quadrature.emplace();
          else
            m_quadrature.emplace(mesh, dimension, order);
          m_offsets = ::Kokkos::View<Index*, MemorySpace>(
            ::Kokkos::view_alloc(execution, ::Kokkos::WithoutInitializing,
              std::string("Rodin entry offsets")),
            count);
          m_mesh = &mesh;
          m_dimension = dimension;
          m_polytopeCount = count;
          m_located.reset();
        }
        if (!m_selection || *m_selection != attributes)
        {
          m_attributes =
            ::Rodin::Kokkos::Geometry::uploadAttributes(execution, attributes);
          m_selection = attributes;
          m_located.reset();
        }
        // The entries of a polytope depend on its attributes, on the
        // entities of the mesh and on the spaces, which are fixed for the
        // integrator owning this object; the offsets of the last integration
        // hold as long as the mesh keeps its revision.
        const Located located{mesh.getRevision(), region, trialVDim, testVDim};
        if (m_located && *m_located == located)
          return m_count;
        m_region = region;
        m_count = locate(execution, mesh.getView(), region, dimension, m_attributes,
          m_offsets, count, trialDOFs, trialVDim, testDOFs, testVDim, m_spacing);
        m_located = located;
        return m_count;
      }

      /// @brief Gets the quadrature tables. @return Tables of the last integration.
      const QuadratureType& getQuadrature() const
      {
        assert(m_quadrature);
        return *m_quadrature;
      }

      /// @brief Gets the selection. @return Device attributes.
      const AttributesType& getAttributes() const
      {
        return m_attributes;
      }

      /// @brief Gets the position of the first entry of every polytope. @return Offsets.
      const ::Kokkos::View<Index*, MemorySpace>& getOffsets() const
      {
        return m_offsets;
      }

      /**
       * @brief Gets the spacing of the entries of the last integration.
       * @return Distance between consecutive entries of one polytope.
       */
      size_t getSpacing() const
      {
        return m_spacing;
      }

      /// @brief Gets the integration region. @return Region of the last integration.
      ::Rodin::Geometry::Region getRegion() const
      {
        return m_region;
      }

      /// @brief Gets the dimension of the polytopes. @return Polytope dimension.
      Index getDimension() const
      {
        return m_dimension;
      }

      /// @brief Gets the number of polytopes of the dimension. @return Polytope count.
      Index getPolytopeCount() const
      {
        return m_polytopeCount;
      }

      /// @brief Gets the entry count of the last integration. @return Count.
      size_t getCount() const
      {
        return m_count;
      }

      /**
       * @brief Number of entries of the local matrix of a polytope, zero
       * when the polytope is not integrated.
       */
      template <class TrialDOFs, class TestDOFs>
      struct EntryCount
      {
        /// @brief Device mesh.
        ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase> view;

        /// @brief Integration region.
        ::Rodin::Geometry::Region region;

        /// @brief Dimension of the integrated polytopes.
        Index dimension;

        /// @brief Selected attributes; empty selects every polytope of the region.
        AttributesType attributes;

        /// @brief Device degrees of freedom of the trial space.
        TrialDOFs trialDOFs;

        /// @brief Vector dimension of the trial space.
        Index trialVDim;

        /// @brief Device degrees of freedom of the test space.
        TestDOFs testDOFs;

        /// @brief Vector dimension of the test space.
        Index testVDim;

        /**
         * @brief Counts the entries of a polytope.
         * @param p Polytope index.
         * @return Number of entries, zero when the polytope is not
         * integrated.
         */
        KOKKOS_INLINE_FUNCTION Index operator()(Index p) const
        {
          const auto polytope = view.getPolytope(dimension, p);
          if (!Rodin::Kokkos::Geometry::isIntegrated(
                view, region, polytope, attributes))
            return 0;
          return trialDOFs.bind(polytope).size() * trialVDim *
            testDOFs.bind(polytope).size() * testVDim;
        }
      };

      /// @brief Counts the integrated polytopes and bounds their entry counts.
      template <class TrialDOFs, class TestDOFs>
      struct EntryRange
      {
        /// @brief Entry count of a polytope.
        EntryCount<TrialDOFs, TestDOFs> count;

        /**
         * @brief Adds a polytope to the reduction.
         * @param p Polytope index.
         * @param m Number of integrated polytopes.
         * @param lo Smallest entry count of an integrated polytope.
         * @param hi Largest entry count of an integrated polytope.
         */
        KOKKOS_INLINE_FUNCTION void operator()(
          Index p, Index& m, Index& lo, Index& hi) const
        {
          const Index n = count(p);
          if (n == 0)
            return;
          m += 1;
          lo = n < lo ? n : lo;
          hi = n > hi ? n : hi;
        }
      };

      /**
       * @brief Writes the entry offsets of the polytopes, one entry per
       * integrated polytope when their entry counts are uniform.
       */
      template <class TrialDOFs, class TestDOFs>
      struct EntryOffsets
      {
        /// @brief Entry count of a polytope.
        EntryCount<TrialDOFs, TestDOFs> count;

        /// @brief Entry offset of every polytope.
        ::Kokkos::View<Index*, MemorySpace> offsets;

        /// @brief Whether the integrated polytopes own the same entry count.
        bool uniform;

        /**
         * @brief Adds a polytope to the exclusive scan of the offsets.
         * @param p Polytope index.
         * @param update Running offset, the offset of the polytope on input.
         * @param final Whether the pass writes the offsets.
         */
        KOKKOS_INLINE_FUNCTION void operator()(Index p, Index& update, bool final) const
        {
          if (final)
            offsets(p) = update;
          const Index n = count(p);
          update += uniform && n > 0 ? 1 : n;
        }
      };

      /**
       * @brief Computes the entry offsets of the polytopes.
       *
       * A polytope owns the product of its trial and test local counts. When
       * every integrated polytope owns the same count and there are more
       * than one, the offset of a polytope is its rank among them and the
       * spacing is their number, as described by getEntryPosition().
       * Otherwise the offset of a polytope is the count of the previous
       * ones and the spacing is one.
       *
       * @param execution Execution space instance.
       * @param view Device mesh.
       * @param region Integration region.
       * @param dimension Dimension of the polytopes.
       * @param attributes Selected attributes; empty selects every polytope of
       * the region.
       * @param offsets Entry offset of every polytope, written.
       * @param count Number of polytopes of the dimension.
       * @param trialDOFs Device degrees of freedom of the trial space.
       * @param trialVDim Vector dimension of the trial space.
       * @param testDOFs Device degrees of freedom of the test space.
       * @param testVDim Vector dimension of the test space.
       * @param spacing Set to the spacing of the entries.
       * @return Number of entries.
       */
      template <class TrialDOFs, class TestDOFs>
      static size_t locate(const ExecutionSpace& execution,
        const ::Rodin::Kokkos::View<::Rodin::Geometry::MeshBase>& view,
        ::Rodin::Geometry::Region region, Index dimension,
        const AttributesType& attributes,
        const ::Kokkos::View<Index*, MemorySpace>& offsets, Index count,
        const TrialDOFs& trialDOFs, Index trialVDim, const TestDOFs& testDOFs,
        Index testVDim, size_t& spacing)
      {
        if (count == 0)
        {
          spacing = 1;
          return 0;
        }
        const EntryCount<TrialDOFs, TestDOFs> local{
          view, region, dimension, attributes, trialDOFs, trialVDim, testDOFs, testVDim};
        Index selected = 0;
        Index smallest = 0;
        Index largest = 0;
        ::Kokkos::parallel_reduce("Rodin entry counts",
          QuadratureType::policy(execution, count),
          EntryRange<TrialDOFs, TestDOFs>{local}, ::Kokkos::Sum<Index>(selected),
          ::Kokkos::Min<Index>(smallest), ::Kokkos::Max<Index>(largest));
        const bool uniform = selected > 1 && smallest == largest;
        Index total = 0;
        ::Kokkos::parallel_scan("Rodin entry offsets",
          QuadratureType::policy(execution, count),
          EntryOffsets<TrialDOFs, TestDOFs>{local, offsets, uniform}, total);
        spacing = uniform ? selected : 1;
        return uniform ? size_t(selected) * smallest : total;
      }

    private:
      /// @brief What the offsets of the last integration were computed for.
      struct Located
      {
        Revision revision;
        ::Rodin::Geometry::Region region;
        Index trialVDim;
        Index testVDim;

        /**
         * @brief Compares two locations member by member.
         * @param other Location to compare with.
         * @return Whether the locations are equal.
         */
        bool operator==(const Located& other) const = default;
      };

      const Mesh* m_mesh = nullptr;
      Optional<Located> m_located;
      ::Rodin::Geometry::Region m_region = ::Rodin::Geometry::Region::Cells;
      Index m_dimension = 0;
      Index m_polytopeCount = 0;
      Optional<QuadratureType> m_quadrature;
      Optional<FlatSet<::Rodin::Geometry::Attribute>> m_selection;
      AttributesType m_attributes;
      ::Kokkos::View<Index*, MemorySpace> m_offsets;
      size_t m_count = 0;
      size_t m_spacing = 1;
  };
}

#endif
