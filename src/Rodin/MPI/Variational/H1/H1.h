/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_VARIATIONAL_H1_H1_H
#define RODIN_MPI_VARIATIONAL_H1_H1_H

/**
 * @file
 * @brief Distributed H1 finite element space specializations on MPI meshes.
 *
 * Provides @ref Rodin::Variational::H1 specializations for
 * @ref Rodin::Geometry::Mesh<Rodin::Context::MPI> for both scalar and
 * vector-valued ranges. Each specialization wraps a rank-local H1 space on
 * the mesh shard and maintains local/global degree-of-freedom mappings
 * consistent across ownership and ghost layers for every mesh entity
 * dimension (vertices, edges, faces, cells).
 */

#include <array>
#include <limits>
#include <numeric>
#include <vector>

#include <boost/mpi/collectives.hpp>

#include "Rodin/Array.h"
#include "Rodin/Serialization/Optional.h"
#include "Rodin/Serialization/Vector.h"
#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Variational/FiniteElementSpace.h"
#include "Rodin/Variational/H1/H1.h"

namespace Rodin::Variational
{
  // =========================================================================
  // H1<K, Scalar, Mesh<MPI>>  —  scalar/complex distributed specialization
  // =========================================================================

  /**
   * @brief Distributed H1 finite element space for MPI meshes (scalar range).
   *
   * Wraps a local shard H1 space and augments it with distributed global
   * indexing. DOFs are generated at every mesh entity dimension:
   *
   *  - d = 0 (vertices):  1 DOF per vertex
   *  - d = 1 (edges, K≥2): K−1 interior DOFs per edge
   *  - d = 2 (faces, K≥3): (K−1)(K−2)/2 interior DOFs per triangle,
   *                         (K−1)² per quadrilateral
   *  - d = 3 (cells, K≥4): tetrahedron/wedge/hexahedron interior DOFs
   *
   * Owned degrees of freedom receive contiguous global indices per rank,
   * and ghost/shared DOFs are synchronized from the owning rank.
   *
   * @tparam K    Polynomial degree (≥ 1).
   * @tparam Scalar Scalar type (Real or Complex).
   */
  template <size_t K, class Scalar>
  class H1<K, Scalar, Geometry::Mesh<Context::MPI>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>,
        H1<K, Scalar, Geometry::Mesh<Context::MPI>>>
  {
    public:
      static_assert(K > 0, "Polynomial degree K must be greater than 0.");

      /**
       * @brief Bidirectional map between local and global DOF indices.
       */
      struct IndexBimap
      {
        /// Local-to-global: `left[localDof] = globalDof`.
        std::vector<Index> left;
        /// Global-to-local: `right[globalDof] = localDof`.
        FlatMap<Index, Index> right;
      };

      /// Context type.
      using ContextType = Context::MPI;

      /// Underlying local shard finite element space type.
      using FESType = H1<K, Scalar, Geometry::Mesh<Context::Local>>;

      /// Scalar type carried by DOF values.
      using ScalarType = typename FESType::ScalarType;

      /// Range type of the space.
      using RangeType = ScalarType;

      /// Distributed mesh type.
      using MeshType = Geometry::Mesh<ContextType>;

      /// Finite element type.
      using ElementType = H1Element<K, RangeType>;

      /// CRTP parent.
      using Parent = FiniteElementSpace<MeshType, H1<K, Scalar, MeshType>>;

      using Parent::getGlobalIndex;

      /**
       * @brief Pullback of a function to the reference polytope.
       */
      template <class FunctionDerived>
      class Pullback
        : public FiniteElementSpacePullbackBase<Pullback<FunctionDerived>>
      {
        public:
          using FunctionType = FunctionBase<FunctionDerived>;

          Pullback(const Geometry::Polytope& polytope, const FunctionType& v)
            : m_polytope(polytope), m_v(v.copy())
          {}

          Pullback(const Pullback&) = default;

          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return getFunction()(p);
          }

          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return getFunction()(res, p);
          }

          constexpr
          const FunctionType& getFunction() const
          {
            assert(m_v);
            return *m_v;
          }

        private:
          Geometry::Polytope m_polytope;
          std::unique_ptr<FunctionType> m_v;
      };

      /**
       * @brief Pushforward of a function from the reference polytope.
       */
      template <class CallableType>
      class Pushforward
        : public FiniteElementSpacePushforwardBase<Pushforward<CallableType>>
      {
        public:
          using FunctionType = CallableType;

          Pushforward(const FunctionType& v)
            : m_v(v)
          {}

          Pushforward(const Pushforward&) = default;

          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          template <class T>
          auto operator()(T& res, const Geometry::Point& p) const
          {
            return getFunction()(res, p.getReferenceCoordinates());
          }

          constexpr
          const FunctionType& getFunction() const
          {
            return m_v;
          }

        private:
          const FunctionType m_v;
      };

      /**
       * @brief Constructs the distributed scalar H1 space on the given mesh.
       *
       * The constructor builds a local shard H1 space, then assigns contiguous
       * global DOF ranges to each rank for every entity dimension d = 0..D.
       * Ghost/shared DOFs are populated by communicating from the owning rank.
       *
       * @param[in] mesh Distributed mesh on which the space is defined.
       */
      H1(std::integral_constant<size_t, K>, const MeshType& mesh)
        : m_mesh(mesh),
          m_fes(std::integral_constant<size_t, K>{}, mesh.getShard())
      {
        build(mesh);
      }

      H1(const H1&) = default;
      H1(H1&&) = default;
      H1& operator=(H1&&) = default;

      /**
       * @brief Returns the local shard finite element space.
       */
      const FESType& getShard() const
      {
        return m_fes;
      }

      /**
       * @brief Returns the global ownership range [begin, end).
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        begin = m_offset;
        end   = m_offset + m_owned;
      }

      /**
       * @brief Returns the global DOF index for a local shard DOF.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        return m_local_to_global.left.at(localIdx);
      }

      /**
       * @brief Returns the local DOF index for a global DOF, or nullopt if not present.
       */
      Optional<Index> getLocalIndex(Index globalIdx) const
      {
        auto it = m_local_to_global.right.find(globalIdx);
        if (it == m_local_to_global.right.end())
          return std::nullopt;
        return it->second;
      }

      /**
       * @brief Returns the finite element on local polytope @f$(d, i)@f$.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_fes.getFiniteElement(d, i);
      }

      /**
       * @brief Returns the global number of DOFs of the distributed space.
       */
      size_t getSize() const override
      {
        return m_globalSize;
      }

      /**
       * @brief Returns the vector dimension (1 for scalar spaces).
       */
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /**
       * @brief Returns the distributed mesh.
       */
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Returns global DOF indices for local polytope @f$(d, i)@f$.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        static thread_local IndexArray s_dofs;
        s_dofs = m_fes.getDOFs(d, i);
        for (auto& dof : s_dofs)
          dof = getGlobalIndex(dof);
        return s_dofs;
      }

      /**
       * @brief Returns the global DOF index for a local basis function.
       */
      Index getGlobalIndex(
          const std::pair<size_t, Index>& p, Index localDof) const override
      {
        const auto& [d, i] = p;
        const Index local = m_fes.getGlobalIndex({ d, i }, localDof);
        return getGlobalIndex(local);
      }

      /**
       * @brief Returns a pullback wrapper for a function on local polytope @f$(d, i)@f$.
       */
      template <class FunctionDerived>
      auto getPullback(
          const std::pair<size_t, Index>& p,
          const FunctionBase<FunctionDerived>& v) const
      {
        const auto& [d, i] = p;
        const auto& mesh = getMesh();
        return Pullback<FunctionDerived>(*mesh.getPolytope(d, i), v);
      }

      /**
       * @brief Returns a pushforward wrapper on local polytope @f$(d, i)@f$.
       */
      template <class CallableType>
      auto getPushforward(
          const std::pair<size_t, Index>&, const CallableType& v) const
      {
        return Pushforward<CallableType>(v);
      }

      /**
       * @brief Returns a pushforward wrapper for an explicit polytope object.
       */
      template <class CallableType>
      auto getPushforward(const Geometry::Polytope&, const CallableType& v) const
      {
        return Pushforward<CallableType>(v);
      }

    private:
      /**
       * @brief Returns the positions in getDOFs(d, i) that correspond to interior
       * (non-shared-with-boundary) DOFs for geometry @p g.
       *
       * These positions are deterministic for a given geometry type and degree K.
       * For example, for a Segment the interior positions are [1, ..., K-1];
       * for a Triangle they are the Fekete positions with i2 > 0, j2 > 0,
       * i2 + j2 < K.
       *
       * @param[in] g Geometry type of the mesh entity.
       * @return Sorted list of interior DOF positions within the entity's closure.
       */
      static std::vector<size_t> interiorPositions(Geometry::Polytope::Type g)
      {
        std::vector<size_t> res;
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            res.push_back(0);
            break;
          }

          case Geometry::Polytope::Type::Segment:
          {
            // Interior positions: 1 .. K-1  (K-1 interior DOFs)
            for (size_t k = 1; k + 1 < K + 1; ++k)
              res.push_back(k);
            break;
          }

          case Geometry::Polytope::Type::Triangle:
          {
            // Fekete ordering: idx(i2, j2) = j2*(K+1) - j2*(j2-1)/2 + i2
            // Interior: i2 > 0, j2 > 0, i2 + j2 < K
            for (size_t j2 = 1; j2 < K; ++j2)
            {
              const size_t rowStart = j2 * (K + 1) - j2 * (j2 - 1) / 2;
              for (size_t i2 = 1; i2 < K - j2; ++i2)
                res.push_back(rowStart + i2);
            }
            break;
          }

          case Geometry::Polytope::Type::Quadrilateral:
          {
            // Tensor GLL grid: idx(i, j) = j*(K+1) + i
            // Interior: 1 <= i <= K-1, 1 <= j <= K-1
            const size_t N1 = K + 1;
            for (size_t j = 1; j < K; ++j)
              for (size_t i = 1; i < K; ++i)
                res.push_back(j * N1 + i);
            break;
          }

          case Geometry::Polytope::Type::Tetrahedron:
          {
            // Fekete tet ordering.
            // tetraTotal = (K+1)(K+2)(K+3)/6
            // For (i,j,k) with i+j+k <= K, i,j,k >= 0:
            //   offset_k = tetraTotal - (K-k+1)(K-k+2)(K-k+3)/6
            //   offset_j = j*(K-k+1) - j*(j-1)/2
            //   idx      = offset_k + offset_j + i
            // Interior: i > 0, j > 0, k > 0, i+j+k < K
            const size_t tetraTotal = (K + 1) * (K + 2) * (K + 3) / 6;
            for (size_t k = 1; k + 2 < K; ++k)          // k <= K-3
            {
              const size_t m_tail   = K - k;
              const size_t tetraTail =
                  (m_tail + 1) * (m_tail + 2) * (m_tail + 3) / 6;
              const size_t offset_k = tetraTotal - tetraTail;
              for (size_t j = 1; j + k + 1 < K; ++j)    // j <= K-k-2
              {
                const size_t offset_j = j * (K - k + 1) - j * (j - 1) / 2;
                for (size_t i = 1; i + j + k < K; ++i)  // i <= K-j-k-1
                  res.push_back(offset_k + offset_j + i);
              }
            }
            break;
          }

          case Geometry::Polytope::Type::Wedge:
          {
            // Wedge DOF ordering: wedgeIdx = s * TriCount + triIdx
            // where s is the z-layer (0..K) and triIdx is the Fekete triangle index.
            // Interior wedge: s in [1, K-1] AND triIdx is interior to the triangle.
            const size_t TriCount = (K + 1) * (K + 2) / 2;

            // Compute interior triangle positions
            std::vector<size_t> triInterior;
            for (size_t j2 = 1; j2 < K; ++j2)
            {
              const size_t rowStart = j2 * (K + 1) - j2 * (j2 - 1) / 2;
              for (size_t i2 = 1; i2 < K - j2; ++i2)
                triInterior.push_back(rowStart + i2);
            }

            for (size_t s = 1; s < K; ++s)
              for (const size_t triIdx : triInterior)
                res.push_back(s * TriCount + triIdx);
            break;
          }

          case Geometry::Polytope::Type::Hexahedron:
          {
            // Tensor GLL grid: hexIdx(i,j,k) = k*N1^2 + j*N1 + i
            // Interior: 1 <= i,j,k <= K-1
            const size_t N1 = K + 1;
            for (size_t k = 1; k < K; ++k)
              for (size_t j = 1; j < K; ++j)
                for (size_t i = 1; i < K; ++i)
                  res.push_back(k * N1 * N1 + j * N1 + i);
            break;
          }

          default:
            break;
        }
        return res;
      }

      static size_t triangleIndex(size_t i, size_t j)
      {
        return j * (K + 1) - j * (j - 1) / 2 + i;
      }

      static Optional<size_t> findOrdinal(
          const std::vector<size_t>& positions,
          size_t position)
      {
        for (size_t i = 0; i < positions.size(); ++i)
        {
          if (positions[i] == position)
            return i;
        }
        return std::nullopt;
      }

      static Optional<std::vector<size_t>> vertexPermutation(
          const std::vector<Index>& localVertexIDs,
          const std::vector<Index>& ownerVertexIDs)
      {
        if (localVertexIDs.size() != ownerVertexIDs.size())
          return std::nullopt;

        std::vector<size_t> localToOwner(localVertexIDs.size());
        for (size_t li = 0; li < localVertexIDs.size(); ++li)
        {
          bool found = false;
          for (size_t oi = 0; oi < ownerVertexIDs.size(); ++oi)
          {
            if (localVertexIDs[li] == ownerVertexIDs[oi])
            {
              localToOwner[li] = oi;
              found = true;
              break;
            }
          }
          if (!found)
            return std::nullopt;
        }
        return localToOwner;
      }

      static Optional<size_t> orientedOwnerOrdinal(
          Geometry::Polytope::Type g,
          const std::vector<size_t>& positions,
          size_t localOrdinal,
          const std::vector<Index>& localVertexIDs,
          const std::vector<Index>& ownerVertexIDs)
      {
        if (localOrdinal >= positions.size())
          return std::nullopt;

        const auto permOpt = vertexPermutation(localVertexIDs, ownerVertexIDs);
        if (!permOpt)
          return std::nullopt;
        const auto& localToOwner = *permOpt;

        switch (g)
        {
          case Geometry::Polytope::Type::Segment:
          {
            if (localToOwner.size() != 2)
              return std::nullopt;
            if (localToOwner[0] == 1 && localToOwner[1] == 0)
              return positions.size() - 1 - localOrdinal;
            return localOrdinal;
          }

          case Geometry::Polytope::Type::Triangle:
          {
            if (localToOwner.size() != 3)
              return std::nullopt;

            const size_t raw = positions[localOrdinal];
            Optional<size_t> li;
            Optional<size_t> lj;
            for (size_t j = 1; j < K; ++j)
            {
              for (size_t i = 1; i < K - j; ++i)
              {
                if (triangleIndex(i, j) == raw)
                {
                  li = i;
                  lj = j;
                  break;
                }
              }
              if (li)
                break;
            }
            if (!li || !lj)
              return std::nullopt;

            std::array<size_t, 3> localBary = {
              K - *li - *lj,
              *li,
              *lj
            };
            std::array<size_t, 3> ownerBary = { 0, 0, 0 };
            for (size_t lv = 0; lv < 3; ++lv)
              ownerBary[localToOwner[lv]] = localBary[lv];

            const size_t ownerRaw = triangleIndex(ownerBary[1], ownerBary[2]);
            return findOrdinal(positions, ownerRaw);
          }

          case Geometry::Polytope::Type::Quadrilateral:
          {
            if (localToOwner.size() != 4)
              return std::nullopt;

            const size_t N1 = K + 1;
            const size_t raw = positions[localOrdinal];
            const size_t i = raw % N1;
            const size_t j = raw / N1;

            using Coord = std::array<size_t, 2>;
            const std::array<Coord, 4> corners = {{
              Coord{ 0, 0 },
              Coord{ K, 0 },
              Coord{ K, K },
              Coord{ 0, K }
            }};

            using Transform = Coord (*)(size_t, size_t);
            const std::array<Transform, 8> transforms = {{
              +[](size_t x, size_t y) -> Coord { return { x, y }; },
              +[](size_t x, size_t y) -> Coord { return { K - x, y }; },
              +[](size_t x, size_t y) -> Coord { return { x, K - y }; },
              +[](size_t x, size_t y) -> Coord { return { K - x, K - y }; },
              +[](size_t x, size_t y) -> Coord { return { y, x }; },
              +[](size_t x, size_t y) -> Coord { return { K - y, x }; },
              +[](size_t x, size_t y) -> Coord { return { y, K - x }; },
              +[](size_t x, size_t y) -> Coord { return { K - y, K - x }; }
            }};

            for (const auto& transform : transforms)
            {
              bool matches = true;
              for (size_t lv = 0; lv < 4; ++lv)
              {
                const Coord mapped =
                    transform(corners[lv][0], corners[lv][1]);
                if (mapped != corners[localToOwner[lv]])
                {
                  matches = false;
                  break;
                }
              }

              if (matches)
              {
                const Coord ownerCoord = transform(i, j);
                const size_t ownerRaw = ownerCoord[1] * N1 + ownerCoord[0];
                return findOrdinal(positions, ownerRaw);
              }
            }
            return std::nullopt;
          }

          default:
            return localOrdinal;
        }
      }

      /**
       * @brief Core DOF-numbering algorithm: builds m_local_to_global.
       *
       * For each entity dimension d = 0..D, owned entities receive contiguous
       * global DOF indices.  Non-owned entities obtain their global DOF
       * indices via a direct owner-push exchange:
       *
       * - The owner iterates @c shard.getHalo(d) and sends
       *   @c (gid, firstGlobalDOF, vertexIDs) to every rank listed as a holder.
       *   The halo is guaranteed to be complete because @c reconcile()
       *   runs an iterative holder-set propagation phase that reaches
       *   non-adjacent ranks sharing the same entity (e.g. a K=2 edge at
       *   the junction of four partitions arranged in a ring).
       *
       * - Each non-owner posts a matching irecv from its owner rank
       *   (recorded in @c shard.getOwner(d)).  It then uses the ordered vertex
       *   IDs to map local edge/face interior DOFs to the owner's orientation.
       *
       * @param[in] mesh The distributed mesh.
       */
      void build(const MeshType& mesh)
      {
        const auto& ctx   = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();
        const auto& shard = mesh.getShard();

        const int P    = comm.size();
        const int rank = comm.rank();
        const size_t D = shard.getDimension();

        // ------------------------------------------------------------------
        // Step 1: count owned DOFs per dimension to compute ownership range.
        // ------------------------------------------------------------------
        std::vector<size_t> ownedPerDim(D + 1, 0);
        for (size_t d = 0; d <= D; ++d)
        {
          const size_t count = shard.getPolytopeCount(d);
          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            if (!shard.isOwned(d, i))
              continue;
            const auto g = shard.getGeometry(d, i);
            ownedPerDim[d] += interiorPositions(g).size();
          }
        }

        m_owned = 0;
        for (size_t d = 0; d <= D; ++d)
          m_owned += ownedPerDim[d];

        const size_t inclusive = boost::mpi::scan(comm, m_owned, std::plus<size_t>());
        m_offset = inclusive - m_owned;

        // ------------------------------------------------------------------
        // Step 2: compute global size via all_reduce.
        // ------------------------------------------------------------------
        boost::mpi::all_reduce(comm, m_owned, m_globalSize, std::plus<size_t>());

        // ------------------------------------------------------------------
        // Step 3: pre-allocate local_to_global.left.
        // ------------------------------------------------------------------
        const size_t localDofCount = m_fes.getSize();
        m_local_to_global.left.assign(
            localDofCount, std::numeric_limits<Index>::max());

        // ------------------------------------------------------------------
        // Step 4: for each dimension, assign owned global indices and
        //         exchange with non-owned counterparts via two-phase
        //         owner-push (see method documentation above).
        //
        // Interior DOFs are numbered contiguously in the owner's orientation.
        // Receivers map their local edge/face orientation to the owner's
        // ordinal before adding firstGlobalDOF.
        // ------------------------------------------------------------------
        Index dimOffset = m_offset;

        // Owner-push message:
        //   (globalEntityID, (firstGlobalDOF, ordered global vertex IDs)).
        //
        // The ordered vertex IDs allow non-owner ranks to install shared
        // high-order edge/face DOFs with the same physical orientation as the
        // owner.
        using EntityMsg = std::pair<Index, std::pair<Index, std::vector<Index>>>;

        for (size_t d = 0; d <= D; ++d)
        {
          const size_t count = shard.getPolytopeCount(d);
          const auto& owner  = shard.getOwner(d);
          const auto* d20 =
              d > 0 ? &shard.getConnectivity().getIncidence(d, 0) : nullptr;

          const auto getOrderedVertexIDs =
            [&](Index entity) -> std::vector<Index>
            {
              std::vector<Index> res;
              if (d20 == nullptr)
                return res;
              const auto& vertices = (*d20)[entity];
              res.reserve(vertices.size());
              for (const Index v : vertices)
                res.push_back(shard.getPolytopeMap(0).left.at(v));
              return res;
            };

          Index dofIdx = dimOffset;

          // ----------------------------------------------------------------
          // 4a. Assign global DOF indices to every owned entity.
          // ----------------------------------------------------------------
          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            if (!shard.isOwned(d, i))
              continue;

            const auto g   = shard.getGeometry(d, i);
            const auto pos = interiorPositions(g);
            if (pos.empty())
              continue;

            const auto& entityDOFs = m_fes.getDOFs(d, i);
            for (const size_t p : pos)
              m_local_to_global.left[entityDOFs(p)] = dofIdx++;
          }

          // ----------------------------------------------------------------
          // 4b. Owner-push: owners send (gid, firstGlobalDOF, vertexIDs) to every
          //     rank listed in halo[d].  Non-owners post irecv from their
          //     owner rank.  The halo is complete after reconcile() runs
          //     the holder-set propagation, so no all_to_all is required.
          // ----------------------------------------------------------------
          const auto& halo_d = shard.getHalo(d);

          std::vector<std::vector<EntityMsg>> push_send(static_cast<size_t>(P));

          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            if (!shard.isOwned(d, i))
              continue;

            const auto g   = shard.getGeometry(d, i);
            const auto pos = interiorPositions(g);
            if (pos.empty())
              continue;

            auto hit = halo_d.find(i);
            if (hit == halo_d.end())
              continue;

            const auto& entityDOFs = m_fes.getDOFs(d, i);
            const Index  firstGDOF = m_local_to_global.left[entityDOFs(pos[0])];
            const Index  gid       = mesh.getGlobalIndex(d, i);
            const auto   orderedVertexIDs = getOrderedVertexIDs(i);

            for (const Index h : hit->second)
            {
              const int rh = static_cast<int>(h);
              if (rh == rank)
                continue;
              push_send[static_cast<size_t>(rh)].push_back(
                  { gid, { firstGDOF, orderedVertexIDs } });
            }
          }

          std::vector<int> need_recv(static_cast<size_t>(P), 0);
          for (Index i = 0; i < static_cast<Index>(count); ++i)
          {
            if (shard.isOwned(d, i))
              continue;
            const auto g = shard.getGeometry(d, i);
            if (interiorPositions(g).empty())
              continue;
            auto oit = owner.find(i);
            if (oit == owner.end())
              continue;
            const int ro = static_cast<int>(oit->second);
            if (ro != rank)
              need_recv[static_cast<size_t>(ro)] = 1;
          }

          std::vector<std::vector<EntityMsg>> push_recv(static_cast<size_t>(P));

          {
            std::vector<boost::mpi::request> reqs;
            reqs.reserve(static_cast<size_t>(2 * P));

            const int tag_push = static_cast<int>(d);

            for (int r = 0; r < P; ++r)
            {
              if (need_recv[static_cast<size_t>(r)])
                reqs.push_back(comm.irecv(r, tag_push,
                    push_recv[static_cast<size_t>(r)]));
            }
            for (int r = 0; r < P; ++r)
            {
              if (!push_send[static_cast<size_t>(r)].empty())
                reqs.push_back(comm.isend(r, tag_push,
                    push_send[static_cast<size_t>(r)]));
            }
            boost::mpi::wait_all(reqs.begin(), reqs.end());
          }

          // ----------------------------------------------------------------
          // 4c. Install received global DOF numbering for non-owned entities.
          // ----------------------------------------------------------------
          for (int r = 0; r < P; ++r)
          {
            for (const auto& [gid, payload] : push_recv[static_cast<size_t>(r)])
            {
              const auto& [firstGDOF, ownerVertexIDs] = payload;
              const auto liOpt = mesh.getLocalIndex(d, gid);
              assert(liOpt);
              const Index li = *liOpt;
              assert(!shard.isOwned(d, li));

              const auto& entityDOFs = m_fes.getDOFs(d, li);
              const auto  g          = shard.getGeometry(d, li);
              const auto  pos        = interiorPositions(g);
              const auto  localVertexIDs = getOrderedVertexIDs(li);

              for (size_t k = 0; k < pos.size(); ++k)
              {
                const auto ownerK = orientedOwnerOrdinal(
                    g, pos, k, localVertexIDs, ownerVertexIDs);
                m_local_to_global.left[entityDOFs(pos[k])] =
                    firstGDOF + static_cast<Index>(ownerK.value_or(k));
              }
            }
          }

          dimOffset += ownedPerDim[d];
        } // end for each dimension d

        // ------------------------------------------------------------------
        // Step 5: build global-to-local reverse map.
        // ------------------------------------------------------------------
#ifndef NDEBUG
        for (size_t local = 0; local < localDofCount; ++local)
          assert(m_local_to_global.left[local] != std::numeric_limits<Index>::max());
#endif

        for (size_t local = 0; local < localDofCount; ++local)
        {
          const Index global = m_local_to_global.left[local];
          m_local_to_global.right.emplace(global, static_cast<Index>(local));
        }
      }

      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;

      size_t m_offset;
      size_t m_owned;
      size_t m_globalSize;
      IndexBimap m_local_to_global;
  };


  // =========================================================================
  // H1<K, SpatialVector<Scalar>, Mesh<MPI>>  —  vector distributed specialization
  // =========================================================================

  /**
   * @brief Distributed vector-valued H1 finite element space for MPI meshes.
   *
   * Each scalar interior DOF node is replicated @p vdim times (once per
   * component). Internally this builds a scalar MPI H1 space and derives
   * the vector DOF mappings as:
   *
   * @f[
   *   \text{globalVecDOF}(q, c) =
   *       \text{vdim} \cdot \text{globalScalarDOF}(q) + c
   * @f]
   *
   * @tparam K    Polynomial degree (≥ 1).
   * @tparam Scalar Underlying scalar type (Real or Complex).
   */
  template <size_t K, class Scalar>
  class H1<K, Math::SpatialVector<Scalar>, Geometry::Mesh<Context::MPI>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>,
        H1<K, Math::SpatialVector<Scalar>, Geometry::Mesh<Context::MPI>>>
  {
    public:
      static_assert(K > 0, "Polynomial degree K must be greater than 0.");

      /**
       * @brief Bidirectional map between local and global DOF indices.
       */
      struct IndexBimap
      {
        std::vector<Index> left;
        FlatMap<Index, Index> right;
      };

      using ContextType = Context::MPI;

      using FESType =
          H1<K, Math::SpatialVector<Scalar>, Geometry::Mesh<Context::Local>>;

      using ScalarFESType = H1<K, Scalar, Geometry::Mesh<Context::MPI>>;

      using ScalarType = Scalar;

      using RangeType = Math::SpatialVector<ScalarType>;

      using MeshType = Geometry::Mesh<ContextType>;

      using ElementType = H1Element<K, RangeType>;

      using Parent =
          FiniteElementSpace<MeshType, H1<K, Math::SpatialVector<Scalar>, MeshType>>;

      using Parent::getGlobalIndex;

      /**
       * @brief Pullback of a function to the reference polytope.
       */
      template <class FunctionDerived>
      class Pullback
        : public FiniteElementSpacePullbackBase<Pullback<FunctionDerived>>
      {
        public:
          using FunctionType = FunctionBase<FunctionDerived>;

          Pullback(const Geometry::Polytope& polytope, const FunctionType& v)
            : m_polytope(polytope), m_v(v.copy())
          {}

          Pullback(const Pullback&) = default;

          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return getFunction()(p);
          }

          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return getFunction()(res, p);
          }

          constexpr
          const FunctionType& getFunction() const
          {
            assert(m_v);
            return *m_v;
          }

        private:
          Geometry::Polytope m_polytope;
          std::unique_ptr<FunctionType> m_v;
      };

      /**
       * @brief Pushforward of a function from the reference polytope.
       */
      template <class CallableType>
      class Pushforward
        : public FiniteElementSpacePushforwardBase<Pushforward<CallableType>>
      {
        public:
          using FunctionType = CallableType;

          Pushforward(const FunctionType& v)
            : m_v(v)
          {}

          Pushforward(const Pushforward&) = default;

          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return getFunction()(p.getReferenceCoordinates());
          }

          template <class T>
          auto operator()(T& res, const Geometry::Point& p) const
          {
            return getFunction()(res, p.getReferenceCoordinates());
          }

          constexpr
          const FunctionType& getFunction() const
          {
            return m_v;
          }

        private:
          const FunctionType m_v;
      };

      /**
       * @brief Constructs a distributed vector-valued H1 space.
       *
       * @param[in] mesh  Distributed mesh on which the space is defined.
       * @param[in] vdim  Number of components per DOF node.
       */
      H1(std::integral_constant<size_t, K>,
         const MeshType& mesh,
         size_t vdim)
        : m_mesh(mesh),
          m_fes(std::integral_constant<size_t, K>{}, mesh.getShard(), vdim),
          m_vdim(vdim)
      {
        // Build the scalar MPI H1 space first.
        ScalarFESType scalarSpace(std::integral_constant<size_t, K>{}, mesh);

        const size_t scalarLocalSize  = scalarSpace.getShard().getSize();
        const size_t scalarGlobalSize = scalarSpace.getSize();

        m_globalSize = scalarGlobalSize * vdim;
        m_owned      = 0;

        {
          Index begin, end;
          scalarSpace.getOwnershipRange(begin, end);
          m_offset = begin * static_cast<Index>(vdim);
          m_owned  = (end - begin) * vdim;
        }

        // ------------------------------------------------------------------
        // Build local_to_global from scalar space.
        //
        // Local vector H1 layout (from H1<K, SpatialVector, Local>):
        //   vecLocal[q*vdim + c] = scalarLocal[q] + c * scalarLocalSize
        //
        // Global vector layout interleaves components per scalar DOF:
        //   vecGlobal = scalarGlobal * vdim + c
        //
        // This keeps each rank's PETSc ownership range contiguous:
        //   [scalarBegin * vdim, scalarEnd * vdim)
        // ------------------------------------------------------------------
        const size_t localVecSize = m_fes.getSize();
        m_local_to_global.left.assign(
            localVecSize, std::numeric_limits<Index>::max());

        for (size_t localScalar = 0; localScalar < scalarLocalSize; ++localScalar)
        {
          const Index globalScalar = scalarSpace.getGlobalIndex(
              static_cast<Index>(localScalar));

          for (size_t c = 0; c < vdim; ++c)
          {
            const Index localVec  =
                static_cast<Index>(localScalar + c * scalarLocalSize);
            const Index globalVec =
                static_cast<Index>(globalScalar * vdim + c);

            m_local_to_global.left[localVec] = globalVec;
          }
        }

#ifndef NDEBUG
        for (size_t local = 0; local < localVecSize; ++local)
          assert(m_local_to_global.left[local] != std::numeric_limits<Index>::max());
#endif

        for (size_t local = 0; local < localVecSize; ++local)
        {
          const Index global = m_local_to_global.left[local];
          m_local_to_global.right.emplace(global, static_cast<Index>(local));
        }
      }

      H1(const H1&) = default;
      H1(H1&&) = default;
      H1& operator=(H1&&) = default;

      /**
       * @brief Returns the local shard finite element space.
       */
      const FESType& getShard() const
      {
        return m_fes;
      }

      /**
       * @brief Returns the global ownership range [begin, end).
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        begin = m_offset;
        end   = m_offset + m_owned;
      }

      /**
       * @brief Returns the global DOF index for a local shard DOF.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        return m_local_to_global.left.at(localIdx);
      }

      /**
       * @brief Returns the local DOF index for a global DOF, or nullopt if not present.
       */
      Optional<Index> getLocalIndex(Index globalIdx) const
      {
        auto it = m_local_to_global.right.find(globalIdx);
        if (it == m_local_to_global.right.end())
          return std::nullopt;
        return it->second;
      }

      /**
       * @brief Returns the finite element on local polytope @f$(d, i)@f$.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_fes.getFiniteElement(d, i);
      }

      /**
       * @brief Returns the global number of DOFs.
       */
      size_t getSize() const override
      {
        return m_globalSize;
      }

      /**
       * @brief Returns the vector dimension of the space.
       */
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      /**
       * @brief Returns the distributed mesh.
       */
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Returns global DOF indices for local polytope @f$(d, i)@f$.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        static thread_local IndexArray s_dofs;
        s_dofs = m_fes.getDOFs(d, i);
        for (auto& dof : s_dofs)
          dof = getGlobalIndex(dof);
        return s_dofs;
      }

      /**
       * @brief Returns the global DOF index for a local basis function.
       */
      Index getGlobalIndex(
          const std::pair<size_t, Index>& p, Index localDof) const override
      {
        const auto& [d, i] = p;
        const Index local = m_fes.getGlobalIndex({ d, i }, localDof);
        return getGlobalIndex(local);
      }

      /**
       * @brief Returns a pullback wrapper for a function on local polytope @f$(d, i)@f$.
       */
      template <class FunctionDerived>
      auto getPullback(
          const std::pair<size_t, Index>& p,
          const FunctionBase<FunctionDerived>& v) const
      {
        const auto& [d, i] = p;
        const auto& mesh = getMesh();
        return Pullback<FunctionDerived>(*mesh.getPolytope(d, i), v);
      }

      /**
       * @brief Returns a pushforward wrapper on local polytope @f$(d, i)@f$.
       */
      template <class CallableType>
      auto getPushforward(
          const std::pair<size_t, Index>&, const CallableType& v) const
      {
        return Pushforward<CallableType>(v);
      }

      /**
       * @brief Returns a pushforward wrapper for an explicit polytope object.
       */
      template <class CallableType>
      auto getPushforward(const Geometry::Polytope&, const CallableType& v) const
      {
        return Pushforward<CallableType>(v);
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;
      size_t m_vdim;

      size_t m_offset;
      size_t m_owned;
      size_t m_globalSize;
      IndexBimap m_local_to_global;
  };
}

namespace Rodin::MPI
{
  /**
   * @brief Convenience alias for the default distributed scalar H1 space.
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  using H1 = Variational::H1<K, Real, Geometry::Mesh<Context::MPI>>;

  /**
   * @brief Convenience alias for the default distributed vector H1 space.
   * @tparam K Polynomial degree.
   */
  template <size_t K>
  using VectorH1 =
      Variational::H1<K, Math::SpatialVector<Real>, Geometry::Mesh<Context::MPI>>;
}

#endif // RODIN_MPI_VARIATIONAL_H1_H1_H
