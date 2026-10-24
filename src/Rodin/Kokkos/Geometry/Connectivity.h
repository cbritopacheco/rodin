/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_GEOMETRY_CONNECTIVITY_H
#define RODIN_KOKKOS_GEOMETRY_CONNECTIVITY_H

/**
 * @file
 * @brief Resident Kokkos incidence discovery and computation.
 */

#include <array>
#include <map>
#include <mutex>
#include <Kokkos_Sort.hpp>
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Kokkos/Context/Kokkos.h"
#include "Rodin/Kokkos/Utility.h"
#include "MeshView.h"

namespace Rodin::Geometry
{
  /**
   * @brief Resident Kokkos graded incidence complex with device computation.
   * @tparam Backend Kokkos execution backend.
   *
   * ## Architecture
   *
   * Entity vertex lists are dimension-major resident arrays. Each incidence
   * relation owns independent CSR arrays. compute() uses Rodin's recursive
   * dimension-pair dispatch, device transpose and ordered intersection.
   * Discovery emits reference subentities, sorts symmetric vertex keys, and
   * scans first occurrences to preserve existing IDs and CPU discovery order.
   * Vertex orientation is retained independently of the symmetric keys.
   *
   * Only counts, allocation sizes and reference-cell constants cross the host
   * boundary during computation; actual mesh topology remains resident.
   * Legacy reference-returning queries materialize read-only host entries on
   * demand. Mutation requires exclusive access and invalidates borrowed views
   * and host references. Copies own independent resident arrays.
   */
  template <class Backend>
  class Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>> final
    : public ConnectivityBase
  {
    public:
      /// Distribution and execution association.
      using ContextType = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;

      /// Native execution space.
      using ExecutionSpace = typename ContextType::ExecutionSpaceType;

      /// Resident memory space.
      using MemorySpace = typename ContextType::MemorySpaceType;

      /// Borrowed device incidence complex.
      using DeviceView = Rodin::Kokkos::Geometry::ConnectivityView;

      /// The existing discovery/restriction vocabulary.
      using Mode = Connectivity<Rodin::Context::Local>::Mode;

      /// @param context Execution association; the runtime must be initialized.
      explicit Connectivity(ContextType context = ContextType())
        : m_context(std::move(context))
      {
        m_entityOffsets = upload(IndexVector{0, 0}, "entity offsets");
        m_vertexOffsets = upload(IndexVector{0}, "vertex offsets");
        publishRelations();
      }

      /**
       * @param context Execution association.
       * @param source Construction topology; it is not retained.
       */
      Connectivity(ContextType context, const ConnectivityBase& source)
        : m_context(std::move(context)),
          m_dimension(source.getDimension())
      {
        if (m_dimension > 3)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Kokkos connectivity supports dimensions zero through three."
            << Alert::Raise;
        }
        IndexVector offsets(m_dimension + 2, 0), vertexOffsets(1, 0), vertices;
        std::vector<Polytope::Type> geometries;
        for (size_t d = 0; d <= m_dimension; ++d)
        {
          offsets[d + 1] = offsets[d] + source.getCount(d);
          for (Index i = 0; i < source.getCount(d); ++i)
          {
            const auto geometry = source.getGeometry(d, i);
            const auto& key = source.getPolytope(d, i);
            const Polytope::Traits traits(geometry);
            if (traits.getDimension() != d || traits.getVertexCount() != key.size())
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Inconsistent input entity." << Alert::Raise;
            }
            geometries.push_back(geometry);
            ++m_geometryCounts[static_cast<size_t>(geometry)];
            for (const Index vertex : key)
            {
              if (vertex >= source.getCount(0))
              {
                Alert::MemberFunctionException(*this, __func__)
                  << "Invalid input vertex index." << Alert::Raise;
              }
              vertices.push_back(vertex);
            }
            vertexOffsets.push_back(vertices.size());
          }
        }
        std::copy(offsets.begin(), offsets.end(), m_offsets.begin());
        m_entityOffsets = upload(offsets, "entity offsets");
        m_geometries = upload(geometries, "entity geometries");
        m_vertexOffsets = upload(vertexOffsets, "vertex offsets");
        m_vertices = upload(vertices, "entity vertices");
        const auto* local =
          dynamic_cast<const Connectivity<Rodin::Context::Local>*>(&source);
        const auto* resident = dynamic_cast<const Connectivity*>(&source);
        for (size_t d = 0; d <= m_dimension; ++d)
        {
          for (size_t dp = 0; dp <= m_dimension; ++dp)
          {
            auto& relation = m_incidence[d][dp];
            if (source.getCount(0) == 0)
              continue;
            const auto& input = source.getIncidence(d, dp);
            relation.dirty = local ? local->isDirty(d, dp)
              : resident           ? resident->isDirty(d, dp)
                                   : input.empty();
            if (input.empty() && relation.dirty)
              continue;
            if (input.size() != source.getCount(d))
            {
              Alert::MemberFunctionException(*this, __func__)
                << "Invalid input incidence row count." << Alert::Raise;
            }
            IndexVector rows(1, 0), indices;
            for (const auto& row : input)
            {
              for (const Index target : row)
              {
                if (target >= source.getCount(dp))
                {
                  Alert::MemberFunctionException(*this, __func__)
                    << "Invalid input incidence target." << Alert::Raise;
                }
                indices.push_back(target);
              }
              rows.push_back(indices.size());
            }
            relation.offsets = upload(rows, "incidence offsets");
            relation.indices = upload(indices, "incidence indices");
            relation.present = true;
          }
        }
        publishRelations();
      }

      /// @param other Connectivity to deep-copy.
      Connectivity(const Connectivity& other)
        : m_context(other.m_context),
          m_dimension(other.m_dimension),
          m_offsets(other.m_offsets),
          m_geometryCounts(other.m_geometryCounts)
      {
        const auto& execution = m_context.getExecutionSpace();
        m_entityOffsets =
          Rodin::Kokkos::Internal::clone(execution, other.m_entityOffsets);
        m_geometries = Rodin::Kokkos::Internal::clone(execution, other.m_geometries);
        m_vertexOffsets =
          Rodin::Kokkos::Internal::clone(execution, other.m_vertexOffsets);
        m_vertices = Rodin::Kokkos::Internal::clone(execution, other.m_vertices);
        for (size_t d = 0; d <= m_dimension; ++d)
        {
          for (size_t dp = 0; dp <= m_dimension; ++dp)
          {
            auto& relation = m_incidence[d][dp];
            const auto& source = other.m_incidence[d][dp];
            relation.present = source.present;
            relation.dirty = source.dirty;
            relation.offsets = Rodin::Kokkos::Internal::clone(execution, source.offsets);
            relation.indices = Rodin::Kokkos::Internal::clone(execution, source.indices);
          }
        }
        publishRelations();
      }

      /// @param other Connectivity whose allocations are transferred.
      Connectivity(Connectivity&& other)
        : Connectivity(other.m_context)
      {
        *this = std::move(other);
      }

      /// @param other Connectivity whose allocations are transferred. @returns This object.
      Connectivity& operator=(Connectivity&& other)
      {
        if (this != &other)
        {
          m_context = other.m_context;
          m_dimension = other.m_dimension;
          m_offsets = other.m_offsets;
          m_geometryCounts = other.m_geometryCounts;
          m_entityOffsets = std::move(other.m_entityOffsets);
          m_geometries = std::move(other.m_geometries);
          m_vertexOffsets = std::move(other.m_vertexOffsets);
          m_vertices = std::move(other.m_vertices);
          m_descriptors = std::move(other.m_descriptors);
          m_incidence = std::move(other.m_incidence);
          m_reference = std::move(other.m_reference);
          m_keys.clear();
          m_relations.clear();
        }
        return *this;
      }

      /// @returns Borrowed resident tables, invalidated by mutation.
      DeviceView getDeviceView() const
      {
        return {m_dimension, m_entityOffsets.data(), m_geometries.data(),
          m_vertexOffsets.data(), m_vertices.data(), m_descriptors.data()};
      }

      /// @returns Associated execution context.
      const ContextType& getContext() const
      {
        return m_context;
      }

      /// @returns Topological dimension.
      size_t getDimension() const override
      {
        return m_dimension;
      }

      /// @returns Maximal dimension supported by this resident complex.
      size_t getMaximalDimension() const
      {
        return m_dimension;
      }

      /// @param d Entity dimension. @returns Entity count, or zero above the dimension.
      size_t getCount(size_t d) const override
      {
        return d <= m_dimension ? m_offsets[d + 1] - m_offsets[d] : 0;
      }

      /// @param geometry Reference family. @returns Entity count.
      size_t getCount(Polytope::Type geometry) const override
      {
        return m_geometryCounts[static_cast<size_t>(geometry)];
      }

      /// @param d Entity dimension. @returns Flattened dimension offset.
      Index getOffset(size_t d) const
      {
        assert(d <= m_dimension + 1);
        return m_offsets[d];
      }

      /// @param d Entity dimension. @param i Local index. @returns Reference family.
      Polytope::Type getGeometry(size_t d, Index i) const override
      {
        assert(i < getCount(d));
        return Rodin::Kokkos::Internal::read(
          m_context.getExecutionSpace(), m_geometries, m_offsets[d] + i);
      }

      /// @param d Entity dimension. @param i Local index. @returns Host key, materialized on demand.
      const Polytope::Key& getPolytope(size_t d, Index i) const override
      {
        assert(i < getCount(d));
        std::lock_guard<std::mutex> lock(m_mutex);
        const Index entity = m_offsets[d] + i;
        auto it = m_keys.find(entity);
        if (it == m_keys.end())
        {
          const auto& execution = m_context.getExecutionSpace();
          const Index first =
            Rodin::Kokkos::Internal::read(execution, m_vertexOffsets, entity);
          const Index end =
            Rodin::Kokkos::Internal::read(execution, m_vertexOffsets, entity + 1);
          const auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{},
            ::Kokkos::subview(m_vertices, std::make_pair(first, end)));
          Polytope::Key key(end - first);
          for (Index a = 0; a < end - first; ++a)
            key[a] = host(a);
          it = m_keys.emplace(entity, std::move(key)).first;
        }
        return it->second;
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Host relation, materialized on demand.
      const Incidence& getIncidence(size_t d, size_t dp) const override
      {
        checkDimensions(d, dp);
        std::lock_guard<std::mutex> lock(m_mutex);
        const auto pair = std::make_pair(d, dp);
        auto it = m_relations.find(pair);
        if (it == m_relations.end())
        {
          Incidence relation;
          const auto& resident = m_incidence[d][dp];
          if (resident.present)
          {
            const auto rows = ::Kokkos::create_mirror_view_and_copy(
              ::Kokkos::HostSpace{}, resident.offsets);
            const auto indices = ::Kokkos::create_mirror_view_and_copy(
              ::Kokkos::HostSpace{}, resident.indices);
            relation.resize(resident.offsets.extent(0) - 1);
            for (Index i = 0; i < relation.size(); ++i)
            {
              for (Index j = rows(i); j < rows(i + 1); ++j)
                relation[i].push_back(indices(j));
            }
          }
          it = m_relations.emplace(pair, std::move(relation)).first;
        }
        return it->second;
      }

      /// @param pair Dimension pair. @param i Source index. @returns Incident indices.
      const IndexVector& getIncidence(
        std::pair<size_t, size_t> pair, Index i) const override
      {
        return getIncidence(pair.first, pair.second).at(i);
      }

      /**
       * @param d Source dimension.
       * @param dp Target dimension.
       * @param mode Discover missing entities, or restrict to existing entities.
       * @returns This connectivity after resident computation completes.
       */
      Connectivity& compute(size_t d, size_t dp, Mode mode = Mode::Discover)
      {
        if (!getCount(0))
          return *this;
        checkDimensions(d, dp);
        const size_t D = getDimension();
        if (d == D && dp == 0)
          return *this;
        if (isDirty(D, D))
          transpose(0, D).intersection(D, D, 0);
        if (d != D && d != 0 && (isDirty(D, d) || isDirty(d, 0)))
          build(d, mode);
        if (dp != D && dp != 0 && (isDirty(D, dp) || isDirty(dp, 0)))
          build(dp, mode);
        if (isDirty(d, dp))
        {
          if (d < dp)
            compute(dp, d, mode).transpose(d, dp);
          else
          {
            const size_t pivot = d == 0 && dp == 0 ? D : 0;
            compute(d, pivot).compute(pivot, dp, mode).intersection(d, dp, pivot);
          }
        }
        return *this;
      }

      /// @param d Source dimension. @param dp Target dimension. @returns This connectivity.
      Connectivity& discover(size_t d, size_t dp)
      {
        return compute(d, dp, Mode::Discover);
      }

      /// @param d Source dimension. @param dp Target dimension. @returns This connectivity.
      Connectivity& restrict(size_t d, size_t dp)
      {
        return compute(d, dp, Mode::Restrict);
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Whether resident incidence arrays are materialized.
      bool hasIncidence(size_t d, size_t dp) const
      {
        checkDimensions(d, dp);
        return m_incidence[d][dp].present;
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Whether recomputation is required.
      bool isDirty(size_t d, size_t dp) const
      {
        checkDimensions(d, dp);
        return m_incidence[d][dp].dirty;
      }

      /// @param d Source dimension. @param dp Target dimension. @param dirty New dirty state. @returns This connectivity.
      Connectivity& setDirty(size_t d, size_t dp, bool dirty = true)
      {
        checkDimensions(d, dp);
        m_incidence[d][dp].dirty = dirty;
        return *this;
      }

      /// @param d Source dimension. @param dp Target dimension. @returns This connectivity after clearing the relation.
      Connectivity& clear(size_t d, size_t dp)
      {
        checkDimensions(d, dp);
        m_incidence[d][dp] = Relation{};
        m_relations.erase({d, dp});
        publishRelations();
        return *this;
      }

      /// @param d Lower dimension. @param dp Upper dimension. @returns This connectivity after device transpose.
      Connectivity& transpose(size_t d, size_t dp);

      /// @param d Source dimension. @param dp Target dimension. @param pivot Intermediate dimension. @returns This connectivity after ordered device composition.
      Connectivity& intersection(size_t d, size_t dp, size_t pivot);

      /// @param d Missing intermediate dimension. @param mode Discovery policy. @returns This connectivity after building cell-to-entity incidence.
      Connectivity& build(size_t d, Mode mode = Mode::Discover);

    private:
      using Descriptor = Rodin::Kokkos::Geometry::IncidenceDescriptor;
      using IncidenceView = Rodin::Kokkos::Geometry::IncidenceView;
      template <class T>
      using DeviceArray = ::Kokkos::View<T*, MemorySpace>;

      /// @brief Independently owned dimension-pair relation.
      struct Relation
      {
        DeviceArray<Index> offsets, indices; ///< Resident CSR arrays.
        bool present = false; ///< Whether arrays are materialized.
        bool dirty = true; ///< Whether compute must rebuild the relation.
      };

      /// @brief Candidate incidence with its deterministic traversal position.
      struct Entry
      {
        Index row, target, order; ///< Source, target and original occurrence.
      };

      /// @brief Groups equal incidences while retaining the earliest occurrence.
      struct EntryLess
      {
        /// @param a First entry. @param b Second entry. @returns Lexicographic ordering.
        KOKKOS_INLINE_FUNCTION bool operator()(const Entry& a, const Entry& b) const
        {
          return a.row != b.row    ? a.row < b.row
            : a.target != b.target ? a.target < b.target
                                   : a.order < b.order;
        }
      };

      /// @brief Ordered vertices of an intermediate reference subentity.
      struct Entity
      {
        ::Kokkos::Array<Index, 4>
          vertices{}; ///< At most four vertices for edges and faces.
        Index size = 0; ///< Active vertex count.
        Polytope::Type geometry = Polytope::Type::Point; ///< Reference family.
      };

      /// @brief Symmetric identity key and original candidate position.
      struct Candidate
      {
        ::Kokkos::Array<Index, 4> vertices{}; ///< Sorted vertex multiset.
        Index size = 0, order = 0; ///< Active count and occurrence.
      };

      /// @brief Orders symmetric keys, breaking ties by occurrence.
      struct CandidateLess
      {
        /// @param a First key. @param b Second key. @returns Lexicographic ordering.
        KOKKOS_INLINE_FUNCTION bool operator()(
          const Candidate& a, const Candidate& b) const
        {
          if (a.size != b.size)
            return a.size < b.size;
          for (Index k = 0; k < a.size; ++k)
          {
            if (a.vertices[k] != b.vertices[k])
              return a.vertices[k] < b.vertices[k];
          }
          return a.order < b.order;
        }
      };

      /// @brief Reference-cell subentities in Rodin local ordering.
      struct Reference
      {
        ::Kokkos::Array<Entity, 12> entities{}; ///< Maximum is twelve hexahedron edges.
        Index size = 0; ///< Number of local subentities.
      };

      /// @brief Exclusive row prefix sum with a terminal entry.
      struct Prefix
      {
        DeviceArray<Index> counts, offsets; ///< Counts and prefix output.

        /// @param i Prefix index. @param sum Running total. @param final Whether output may be written.
        KOKKOS_INLINE_FUNCTION void operator()(Index i, Index& sum, bool final) const
        {
          if (final)
            offsets(i) = sum;
          if (i < counts.extent(0))
            sum += counts(i);
        }
      };

      /// @brief Stable duplicate removal shared by transpose and composition.
      struct UniqueEntries
      {
        struct Mark
        {}; ///< Mark first occurrences.
        struct Fill
        {}; ///< Fill compact indices and count rows.
        DeviceArray<Entry> entries; ///< Sorted candidate incidences.
        DeviceArray<Index> flags, positions, indices, counts; ///< Compaction arrays.
        bool ordered; ///< Whether to retain original traversal order.

        /// @param i Sorted candidate index. @returns Whether it is the first incidence of its pair.
        KOKKOS_INLINE_FUNCTION bool first(Index i) const
        {
          return i == 0 || entries(i).row != entries(i - 1).row ||
            entries(i).target != entries(i - 1).target;
        }

        /// @param tag Mark pass. @param i Sorted candidate index.
        KOKKOS_INLINE_FUNCTION void operator()([[maybe_unused]] Mark tag, Index i) const
        {
          if (first(i))
            flags(ordered ? entries(i).order : i) = 1;
        }

        /// @param tag Fill pass. @param i Sorted candidate index.
        KOKKOS_INLINE_FUNCTION void operator()([[maybe_unused]] Fill tag, Index i) const
        {
          if (first(i))
          {
            indices(positions(ordered ? entries(i).order : i)) = entries(i).target;
            ::Kokkos::atomic_inc(&counts(entries(i).row));
          }
        }
      };

      /// @brief Emission of transpose pairs, parallel over original source rows.
      struct TransposeEntries
      {
        IncidenceView source; ///< Original CSR relation.
        DeviceArray<Entry> entries; ///< Transposed pairs.

        /// @param j Original source index.
        KOKKOS_INLINE_FUNCTION void operator()(Index j) const
        {
          for (Index k = source.offsets[j]; k < source.offsets[j + 1]; ++k)
            entries(k) = {source.indices[k], j, k};
        }
      };

      /// @brief Counts or emits Rodin's ordered, subset-filtered composition paths.
      struct ComposeEntries
      {
        DeviceView topology; ///< Resident ordered vertex lists.
        IncidenceView left, right; ///< Two relations through the pivot dimension.
        DeviceArray<Index> counts, offsets; ///< Candidate row lengths and offsets.
        DeviceArray<Entry> entries; ///< Emitted pairs.
        size_t dimension, targetDimension; ///< Source and target dimensions.
        bool fill; ///< Whether this is the emission pass.

        /// @param i Source entity index.
        KOKKOS_INLINE_FUNCTION void operator()(Index i) const
        {
          Index n = 0;
          const auto vertices = topology.getPolytope(dimension, i);
          const auto pivots = left[i];
          for (Index k = 0; k < pivots.size(); ++k)
          {
            const auto candidates = right[pivots[k]];
            for (Index a = 0; a < candidates.size(); ++a)
            {
              const Index j = candidates[a];
              if (dimension == targetDimension && i == j)
                continue;
              bool subset = true;
              if (dimension > targetDimension)
              {
                const auto targetVertices = topology.getPolytope(targetDimension, j);
                for (Index v = 0; v < targetVertices.size(); ++v)
                {
                  bool found = false;
                  for (Index u = 0; u < vertices.size(); ++u)
                    found = found || vertices[u] == targetVertices[v];
                  subset = subset && found;
                }
              }
              if (subset)
              {
                if (fill)
                  entries(offsets(i) + n) = {i, j, offsets(i) + n};
                ++n;
              }
            }
          }
          if (!fill)
            counts(i) = n;
        }
      };

      /// @brief Resident discovery stages and their shared workspace.
      struct Discovery
      {
        struct Count
        {}; ///< Count local reference subentities.
        struct Existing
        {}; ///< Emit existing identities first.
        struct Emit
        {}; ///< Emit cell-local candidates.
        struct Groups
        {}; ///< Mark canonical-key groups.
        struct Leaders
        {}; ///< Choose first occurrences.
        struct Assign
        {}; ///< Assign stable entity IDs.
        struct Resolve
        {}; ///< Resolve candidate identities.
        struct NewEntities
        {}; ///< Compact newly discovered entities.
        struct CountRows
        {}; ///< Count retained cell-local incidences.
        struct FillRows
        {}; ///< Fill cell-local incidence in reference order.
        DeviceView topology; ///< Input resident topology.
        DeviceArray<Reference> reference; ///< Rodin reference-cell constants.
        DeviceArray<Index> counts, offsets; ///< Emission row lengths and offsets.
        DeviceArray<Entity> records, added; ///< Oriented candidates and new entities.
        DeviceArray<Candidate> candidates; ///< Canonical identity keys.
        DeviceArray<Index> flags, groupOffsets, groups, leaders, leaderFlags,
          leaderOffsets, assigned, ids; ///< Grouping and numbering workspace.
        DeviceArray<unsigned char> matched; ///< Explicit membership of each occurrence.
        DeviceArray<Index> rowOffsets, rowIndices; ///< Output cell-to-entity relation.
        size_t dimension; ///< Intermediate entity dimension.
        Index oldCount; ///< Existing entities, which retain their indices.
        bool discover; ///< Whether absent keys create new entities.

        /// @param entity Ordered entity. @param order Original occurrence. @returns Symmetric key with deterministic tie-breaking.
        KOKKOS_INLINE_FUNCTION static Candidate key(const Entity& entity, Index order)
        {
          Candidate result{entity.vertices, entity.size, order};
          for (Index a = 1; a < result.size; ++a)
          {
            const Index value = result.vertices[a];
            Index b = a;
            while (b && value < result.vertices[b - 1])
            {
              result.vertices[b] = result.vertices[b - 1];
              --b;
            }
            result.vertices[b] = value;
          }
          return result;
        }

        /// @param a First candidate. @param b Second candidate. @returns Whether their vertex multisets agree.
        KOKKOS_INLINE_FUNCTION static bool equal(const Candidate& a, const Candidate& b)
        {
          if (a.size != b.size)
            return false;
          for (Index k = 0; k < a.size; ++k)
          {
            if (a.vertices[k] != b.vertices[k])
              return false;
          }
          return true;
        }

        /// @param tag Count pass. @param cell Cell index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Count tag, Index cell) const
        {
          counts(cell) = reference(
            static_cast<Index>(topology.getGeometry(topology.dimension, cell)) * 4 +
            dimension)
                           .size;
        }

        /// @param tag Existing-entity pass. @param i Existing index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Existing tag, Index i) const
        {
          const auto vertices = topology.getPolytope(dimension, i);
          Entity entity;
          entity.size = vertices.size();
          entity.geometry = topology.getGeometry(dimension, i);
          for (Index a = 0; a < entity.size; ++a)
            entity.vertices[a] = vertices[a];
          records(i) = entity;
          candidates(i) = key(entity, i);
        }

        /// @param tag Candidate emission. @param cell Cell index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Emit tag, Index cell) const
        {
          const auto& local = reference(
            static_cast<Index>(topology.getGeometry(topology.dimension, cell)) * 4 +
            dimension);
          const auto vertices = topology.getPolytope(topology.dimension, cell);
          for (Index a = 0; a < local.size; ++a)
          {
            Entity entity = local.entities[a];
            for (Index v = 0; v < entity.size; ++v)
              entity.vertices[v] = vertices[entity.vertices[v]];
            const Index order = oldCount + offsets(cell) + a;
            records(order) = entity;
            candidates(order) = key(entity, order);
          }
        }

        /// @param tag Group detection. @param i Sorted candidate index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Groups tag, Index i) const
        {
          flags(i) = i == 0 || !equal(candidates(i), candidates(i - 1));
        }

        /// @param tag First-occurrence selection. @param i Sorted candidate index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Leaders tag, Index i) const
        {
          const Index group = groupOffsets(i + 1) - 1;
          const Index order = candidates(i).order;
          groups(order) = group;
          if (flags(i))
          {
            leaders(group) = order;
            if (order >= oldCount)
              leaderFlags(order) = 1;
          }
        }

        /// @param tag Stable numbering. @param i Original candidate index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Assign tag, Index i) const
        {
          if (i < oldCount)
            assigned(i) = i;
          else if (leaderFlags(i))
            assigned(i) = oldCount + leaderOffsets(i);
        }

        /// @param tag Identity resolution. @param i Original candidate index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] Resolve tag, Index i) const
        {
          const Index leader = leaders(groups(i));
          matched(i) = discover || leader < oldCount;
          if (matched(i))
            ids(i) = assigned(leader);
        }

        /// @param tag New-entity compaction. @param i Original candidate index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] NewEntities tag, Index i) const
        {
          if (leaderFlags(i))
            added(assigned(i) - oldCount) = records(i);
        }

        /// @param tag Restricted row count. @param cell Cell index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] CountRows tag, Index cell) const
        {
          Index n = 0;
          for (Index a = offsets(cell); a < offsets(cell + 1); ++a)
            n += matched(oldCount + a) != 0;
          counts(cell) = n;
        }

        /// @param tag Reference-ordered row fill. @param cell Cell index.
        KOKKOS_INLINE_FUNCTION void operator()(
          [[maybe_unused]] FillRows tag, Index cell) const
        {
          Index out = rowOffsets(cell);
          for (Index a = offsets(cell); a < offsets(cell + 1); ++a)
          {
            if (matched(oldCount + a))
              rowIndices(out++) = ids(oldCount + a);
          }
        }
      };

      /// @brief Resident entity-table insertion while preserving all existing local IDs.
      struct PackEntities
      {
        DeviceView source; ///< Old dimension-major tables.
        DeviceArray<Entity> added; ///< Newly discovered entities.
        DeviceArray<Polytope::Type> geometries; ///< New geometry table.
        DeviceArray<Index> counts, offsets, vertices; ///< New vertex lists.
        Index insertion; ///< Old end of the growing dimension.
        bool fill; ///< Whether vertex storage is allocated.

        /// @param i New flattened entity index.
        KOKKOS_INLINE_FUNCTION void operator()(Index i) const
        {
          const bool inserted = i >= insertion && i < insertion + added.extent(0);
          const Index old = i < insertion ? i : i - added.extent(0);
          if (!fill)
          {
            counts(i) = inserted
              ? added(i - insertion).size
              : source.vertexOffsets[old + 1] - source.vertexOffsets[old];
            geometries(i) =
              inserted ? added(i - insertion).geometry : source.geometries[old];
          }
          else
          {
            for (Index a = 0; a < counts(i); ++a)
            {
              vertices(offsets(i) + a) = inserted
                ? added(i - insertion).vertices[a]
                : source.vertexIndices[source.vertexOffsets[old] + a];
            }
          }
        }
      };

      /// @brief Counts the geometry families in newly discovered entities.
      struct CountGeometry
      {
        DeviceArray<Entity> entities; ///< New entities.
        Polytope::Type geometry; ///< Selected family.

        /// @param i Entity index. @param count Geometry count.
        KOKKOS_INLINE_FUNCTION void operator()(Index i, Index& count) const
        {
          count += entities(i).geometry == geometry;
        }
      };

      /// @brief Counts and copies ordered entity vertex lists into an independent CSR relation.
      struct VertexEntries
      {
        DeviceView topology; ///< Ordered keys.
        size_t dimension; ///< Source dimension.
        DeviceArray<Index> counts, offsets, indices; ///< Relation arrays.
        bool fill; ///< Whether offsets have been computed.

        /// @param i Entity index.
        KOKKOS_INLINE_FUNCTION void operator()(Index i) const
        {
          const auto vertices = topology.getPolytope(dimension, i);
          if (!fill)
            counts(i) = vertices.size();
          else
          {
            for (Index a = 0; a < vertices.size(); ++a)
              indices(offsets(i) + a) = vertices[a];
          }
        }
      };

      /// @param count Number of work items. @returns Execution-instance policy.
      auto policy(Index count) const
      {
        return ::Kokkos::RangePolicy<ExecutionSpace>(
          m_context.getExecutionSpace(), 0, count);
      }

      /// @param count Number of entries. @param label Allocation label. @returns Resident array.
      template <class T>
      DeviceArray<T> allocate(Index count, const char* label) const
      {
        return DeviceArray<T>(
          ::Kokkos::view_alloc(m_context.getExecutionSpace(), std::string(label)), count);
      }

      /// @tparam T Entry type. @param values Construction values. @param label Allocation label. @returns Resident array.
      template <class T>
      DeviceArray<T> upload(const std::vector<T>& values, const char* label) const
      {
        return Rodin::Kokkos::Internal::upload(
          m_context.getExecutionSpace(), values, label);
      }

      /// @param d First dimension. @param dp Second dimension.
      void checkDimensions(size_t d, size_t dp) const
      {
        if (d > m_dimension || dp > m_dimension)
        {
          Alert::MemberFunctionException(*this, __func__)
            << "Invalid incidence dimensions." << Alert::Raise;
        }
      }

      /// @brief Publishes the small resident descriptor table after relation replacement.
      void publishRelations()
      {
        std::vector<Descriptor> descriptors;
        for (size_t d = 0; d <= m_dimension; ++d)
        {
          for (size_t dp = 0; dp <= m_dimension; ++dp)
          {
            const auto& relation = m_incidence[d][dp];
            descriptors.push_back({static_cast<std::uint8_t>(relation.present),
              relation.offsets.data(), relation.indices.data(),
              relation.offsets.extent(0) ? relation.offsets.extent(0) - 1 : 0});
          }
        }
        m_descriptors = upload(descriptors, "incidence descriptors");
      }

      /// @param d Source dimension. @param dp Target dimension. @returns Materialized resident relation.
      IncidenceView relation(size_t d, size_t dp) const
      {
        const auto& r = m_incidence[d][dp];
        assert(r.present);
        return {r.offsets.data(), r.indices.data(), r.offsets.extent(0) - 1};
      }

      /// @param counts Resident row counts. @param offsets Output offsets including the terminal entry. @returns Total entry count.
      Index scan(const DeviceArray<Index>& counts, DeviceArray<Index>& offsets) const;

      /// @param d Source dimension. @param dp Target dimension. @param entries Candidate incidences. @param ordered Preserve original traversal order, rather than sorted target order.
      void finishRelation(
        size_t d, size_t dp, const DeviceArray<Entry>& entries, bool ordered);

      /// @brief Initializes reference-cell constants without examining resident mesh data.
      void initializeReference();

      /// @param d Growing dimension. @param entities New entities in discovery order.
      void appendEntities(size_t d, const DeviceArray<Entity>& entities);

      /// @param d Entity dimension. @brief Materializes ordered entity-to-vertex incidence.
      void vertexRelation(size_t d);

      ContextType m_context;
      size_t m_dimension = 0;
      std::array<Index, 5> m_offsets{};
      std::array<size_t, 8> m_geometryCounts{};
      DeviceArray<Index> m_entityOffsets, m_vertexOffsets, m_vertices;
      DeviceArray<Polytope::Type> m_geometries;
      DeviceArray<Descriptor> m_descriptors;
      std::array<std::array<Relation, 4>, 4> m_incidence;
      DeviceArray<Reference> m_reference;
      mutable std::mutex m_mutex;
      mutable std::map<Index, Polytope::Key> m_keys;
      mutable std::map<std::pair<size_t, size_t>, Incidence> m_relations;
  };

  template <class Backend>
  Index Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::scan(
    const DeviceArray<Index>& counts, DeviceArray<Index>& offsets) const
  {
    offsets = allocate<Index>(counts.extent(0) + 1, "CSR offsets");
    Index total = 0;
    ::Kokkos::parallel_scan(
      "Rodin CSR prefix", policy(counts.extent(0) + 1), Prefix{counts, offsets}, total);
    m_context.getExecutionSpace().fence("Rodin CSR sizes");
    return total;
  }

  template <class Backend>
  void
  Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::finishRelation(
    size_t d, size_t dp, const DeviceArray<Entry>& entries, bool ordered)
  {
    const auto& execution = m_context.getExecutionSpace();
    ::Kokkos::sort(execution, entries, EntryLess{});
    UniqueEntries work{entries, allocate<Index>(entries.extent(0), "unique flags"), {},
      {}, allocate<Index>(getCount(d), "CSR counts"), ordered};
    ::Kokkos::parallel_for("Rodin incidence duplicates",
      ::Kokkos::RangePolicy<ExecutionSpace, typename UniqueEntries::Mark>(
        execution, 0, entries.extent(0)),
      work);
    const Index total = scan(work.flags, work.positions);
    work.indices = allocate<Index>(total, "CSR indices");
    ::Kokkos::parallel_for("Rodin incidence fill",
      ::Kokkos::RangePolicy<ExecutionSpace, typename UniqueEntries::Fill>(
        execution, 0, entries.extent(0)),
      work);
    Relation next;
    scan(work.counts, next.offsets);
    next.indices = std::move(work.indices);
    next.present = true;
    next.dirty = false;
    m_incidence[d][dp] = std::move(next);
    m_relations.erase({d, dp});
    publishRelations();
  }

  template <class Backend>
  auto Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::transpose(
    size_t d, size_t dp) -> Connectivity&
  {
    checkDimensions(d, dp);
    assert(d < dp);
    auto entries =
      allocate<Entry>(m_incidence[dp][d].indices.extent(0), "transpose pairs");
    ::Kokkos::parallel_for("Rodin transpose emit", policy(getCount(dp)),
      TransposeEntries{relation(dp, d), entries});
    finishRelation(d, dp, entries, false);
    return *this;
  }

  template <class Backend>
  auto Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::intersection(
    size_t d, size_t dp, size_t pivot) -> Connectivity&
  {
    checkDimensions(d, dp);
    checkDimensions(d, pivot);
    assert(d >= dp);
    ComposeEntries work{getDeviceView(), relation(d, pivot), relation(pivot, dp),
      allocate<Index>(getCount(d), "composition counts"), {}, {}, d, dp, false};
    ::Kokkos::parallel_for("Rodin composition count", policy(getCount(d)), work);
    const Index total = scan(work.counts, work.offsets);
    work.entries = allocate<Entry>(total, "composition pairs");
    work.fill = true;
    ::Kokkos::parallel_for("Rodin composition emit", policy(getCount(d)), work);
    finishRelation(d, dp, work.entries, true);
    return *this;
  }

  template <class Backend>
  void Connectivity<
    Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::initializeReference()
  {
    if (m_reference.extent(0))
      return;
    // These are reference-cell constants, not a host mirror of the actual mesh.
    static const auto reference = [] {
      std::vector<Reference> result(8 * 4);
      for (size_t g = 0; g < 8; ++g)
      {
        const auto geometry = static_cast<Polytope::Type>(g);
        const Polytope::Traits traits(geometry);
        if (traits.getDimension() < 2)
          continue;
        Connectivity<Rodin::Context::Local> local;
        Polytope::Key vertices(traits.getVertexCount());
        for (Index a = 0; a < vertices.size(); ++a)
          vertices[a] = a;
        local.initialize(traits.getDimension())
          .nodes(vertices.size())
          .polytope(geometry, vertices);
        for (size_t d = 1; d < traits.getDimension(); ++d)
        {
          std::vector<typename Connectivity<Rodin::Context::Local>::SubPolytope> entities;
          local.getSubPolytopes(entities, 0, d);
          auto& r = result[g * 4 + d];
          r.size = entities.size();
          assert(r.size <= r.entities.size());
          for (Index a = 0; a < r.size; ++a)
          {
            auto& out = r.entities[a];
            out.geometry = entities[a].geometry;
            const auto& key = entities[a].vertices;
            out.size = key.size();
            assert(out.size <= out.vertices.size());
            for (Index v = 0; v < out.size; ++v)
              out.vertices[v] = key[v];
          }
        }
      }
      return result;
    }();
    m_reference = upload(reference, "Rodin reference subentities");
  }

  template <class Backend>
  auto Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::build(
    size_t d, Mode mode) -> Connectivity&
  {
    assert(d > 0 && d < m_dimension);
    initializeReference();
    const auto& execution = m_context.getExecutionSpace();
    Discovery work;
    work.topology = getDeviceView();
    work.reference = m_reference;
    work.dimension = d;
    work.oldCount = getCount(d);
    work.discover = mode == Mode::Discover;
    work.counts = allocate<Index>(getCount(m_dimension), "local entity counts");
    ::Kokkos::parallel_for("Rodin discovery count",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Count>(
        execution, 0, getCount(m_dimension)),
      work);
    const Index emitted = scan(work.counts, work.offsets);
    const Index total = work.oldCount + emitted;
    work.records = allocate<Entity>(total, "oriented entity candidates");
    work.candidates = allocate<Candidate>(total, "symmetric entity keys");
    ::Kokkos::parallel_for("Rodin discovery existing",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Existing>(
        execution, 0, work.oldCount),
      work);
    ::Kokkos::parallel_for("Rodin discovery emit",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Emit>(
        execution, 0, getCount(m_dimension)),
      work);
    ::Kokkos::sort(execution, work.candidates, CandidateLess{});
    work.flags = allocate<Index>(total, "entity group flags");
    ::Kokkos::parallel_for("Rodin discovery groups",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Groups>(
        execution, 0, total),
      work);
    const Index groupCount = scan(work.flags, work.groupOffsets);
    work.groups = allocate<Index>(total, "candidate group IDs");
    work.leaders = allocate<Index>(groupCount, "entity first occurrences");
    work.leaderFlags = allocate<Index>(total, "new entity flags");
    ::Kokkos::parallel_for("Rodin discovery leaders",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Leaders>(
        execution, 0, total),
      work);
    const Index newCount = scan(work.leaderFlags, work.leaderOffsets);
    work.assigned = allocate<Index>(total, "stable entity IDs");
    work.ids = allocate<Index>(total, "resolved entity IDs");
    work.matched = allocate<unsigned char>(total, "entity membership");
    ::Kokkos::parallel_for("Rodin discovery assign",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Assign>(
        execution, 0, total),
      work);
    ::Kokkos::parallel_for("Rodin discovery resolve",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::Resolve>(
        execution, 0, total),
      work);
    if (work.discover && newCount)
    {
      work.added = allocate<Entity>(newCount, "new ordered entities");
      ::Kokkos::parallel_for("Rodin discovery compact",
        ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::NewEntities>(
          execution, 0, total),
        work);
      appendEntities(d, work.added);
    }
    ::Kokkos::parallel_for("Rodin discovery row counts",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::CountRows>(
        execution, 0, getCount(m_dimension)),
      work);
    const Index incidenceCount = scan(work.counts, work.rowOffsets);
    work.rowIndices = allocate<Index>(incidenceCount, "cell subentity indices");
    ::Kokkos::parallel_for("Rodin discovery row fill",
      ::Kokkos::RangePolicy<ExecutionSpace, typename Discovery::FillRows>(
        execution, 0, getCount(m_dimension)),
      work);
    execution.fence("Rodin discovered incidence");
    auto& r = m_incidence[m_dimension][d];
    r = Relation{work.rowOffsets, work.rowIndices, true, false};
    m_relations.erase({m_dimension, d});
    vertexRelation(d);
    publishRelations();
    return *this;
  }

  template <class Backend>
  void
  Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::appendEntities(
    size_t d, const DeviceArray<Entity>& entities)
  {
    const auto& execution = m_context.getExecutionSpace();
    const Index count = m_offsets[m_dimension + 1] + entities.extent(0);
    PackEntities work{getDeviceView(), entities,
      allocate<Polytope::Type>(count, "entity geometries"),
      allocate<Index>(count, "entity vertex counts"), {}, {}, m_offsets[d + 1], false};
    ::Kokkos::parallel_for("Rodin entity table insertion", policy(count), work);
    const Index vertices = scan(work.counts, work.offsets);
    work.vertices = allocate<Index>(vertices, "ordered entity vertices");
    work.fill = true;
    ::Kokkos::parallel_for("Rodin entity vertex insertion", policy(count), work);
    for (size_t g = 0; g < m_geometryCounts.size(); ++g)
    {
      Index added = 0;
      ::Kokkos::parallel_reduce("Rodin discovered geometry counts",
        policy(entities.extent(0)),
        CountGeometry{entities, static_cast<Polytope::Type>(g)}, added);
      m_geometryCounts[g] += added;
    }
    execution.fence("Rodin topology growth");
    for (size_t k = d + 1; k <= m_dimension + 1; ++k)
      m_offsets[k] += entities.extent(0);
    m_entityOffsets =
      upload(IndexVector(m_offsets.begin(), m_offsets.begin() + m_dimension + 2),
        "entity offsets");
    m_geometries = std::move(work.geometries);
    m_vertexOffsets = std::move(work.offsets);
    m_vertices = std::move(work.vertices);
    // Any relation touching the grown dimension must acquire the new rows/targets.
    for (size_t k = 0; k <= m_dimension; ++k)
    {
      m_incidence[d][k] = Relation{};
      m_incidence[k][d] = Relation{};
      m_relations.erase({d, k});
      m_relations.erase({k, d});
    }
    m_keys.clear();
  }

  template <class Backend>
  void
  Connectivity<Rodin::Context::Kokkos<Rodin::Context::Local, Backend>>::vertexRelation(
    size_t d)
  {
    VertexEntries work{getDeviceView(), d,
      allocate<Index>(getCount(d), "vertex incidence counts"), {}, {}, false};
    ::Kokkos::parallel_for("Rodin vertex incidence count", policy(getCount(d)), work);
    const Index total = scan(work.counts, work.offsets);
    work.indices = allocate<Index>(total, "vertex incidence indices");
    work.fill = true;
    ::Kokkos::parallel_for("Rodin vertex incidence fill", policy(getCount(d)), work);
    m_context.getExecutionSpace().fence("Rodin vertex incidence");
    m_incidence[d][0] = Relation{work.offsets, work.indices, true, false};
    m_relations.erase({d, 0});
  }

}

#endif
