/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_VARIATIONAL_P0_P0_H
#define RODIN_MPI_VARIATIONAL_P0_P0_H

/**
 * @file
 * @brief Distributed P0 finite element space specialization on MPI meshes.
 *
 * This file provides the @ref Rodin::Variational::P0 specialization for
 * @ref Rodin::Geometry::Mesh<Rodin::Context::MPI>. The implementation builds
 * a rank-local P0 space on each shard and maintains local/global degree-of-
 * freedom mappings consistent across ownership and ghost layers.
 *
 * P0 degrees of freedom are associated with mesh cells: each cell carries
 * one DOF per field component. Cells are partitioned across ranks; owned cells receive
 * contiguous global indices, and ghost cells synchronize their global DOF
 * index from the owning rank via a direct owner-push exchange.
 */

#include <map>
#include <limits>
#include <algorithm>
#include <type_traits>
#include <vector>
#include <utility>

#include <boost/mpi/collectives.hpp>

#include "Rodin/Array.h"

#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Variational/FiniteElementSpace.h"

#include "Rodin/Variational/P0/P0.h"

namespace Rodin::Variational
{
  /**
   * @brief Distributed P0 finite element space for MPI meshes.
   *
   * This specialization wraps a local shard P0 space and augments it with
   * distributed global indexing. Each cell of the distributed mesh carries
   * exactly one degree of freedom per field component. Owned cells receive contiguous global
   * indices per rank, and ghost cells synchronize their global DOF index
   * from the owning rank.
   *
   * @tparam Range Scalar or vector range descriptor used by P0.
   */
  template <class Range>
  class P0<Range, Geometry::Mesh<Context::MPI>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>, P0<Range, Geometry::Mesh<Context::MPI>>>
  {
    public:
      /// @brief Bidirectional map between local and global DOF indices.
      struct IndexBimap
      {
        /**
         * @note The names `left` and `right` mirror bidirectional-map
         * terminology:
         * - `left[local] = global`
         * - `right[global] = local`
         */
        /// Local-to-global mapping (indexed by local DOF = local cell index).
        std::vector<Index> left;
        /// Global-to-local mapping (keyed by global DOF).
        FlatMap<Index, Index> right;
      };

      /// Represents the Context of the P0 space.
      using ContextType = Context::MPI;

      /// Underlying local shard finite element space type.
      using FESType = P0<Range, Geometry::Mesh<Context::Local>>;

      /// Scalar type carried by DOF values.
      using ScalarType = typename FESType::ScalarType;

      /// Range type of value.
      using RangeType = typename FESType::RangeType;

      /// Type of mesh on which the finite element space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element.
      using ElementType = P0Element<RangeType>;

      /// Parent class.
      using Parent = FiniteElementSpace<MeshType, P0<RangeType, MeshType>>;

      using Parent::getGlobalIndex;

      /**
       * @brief Pullback of a function to the reference polytope.
       *
       * Given a function @f$ v @f$ defined on the physical polytope @f$ K @f$,
       * this class defines the function @f$ \hat v @f$ on the reference
       * polytope @f$ \hat K @f$ by
       *
       * @f[
       *   \hat v(\hat x) = v(x(\hat x)),
       * @f]
       *
       * where @f$ x(\hat x) @f$ is the point of @f$ K @f$ with reference
       * coordinates @f$ \hat x @f$.
       */
      template <class Callable>
      class Pullback : public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          /// @brief Callable type evaluated on physical points.
          using CallableType = Callable;

          /**
           * @brief Constructs a pullback on a physical polytope.
           * @tparam Function Callable type accepted by the wrapper.
           * @param polytope Physical polytope.
           * @param v Callable evaluated on physical points.
           */
          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          /**
           * @brief Copy constructor.
           * @param other Object to copy from.
           */
          Pullback(const Pullback& other) = default;

          /**
           * @brief Evaluates the pulled-back callable.
           * @param r Reference coordinates at which to evaluate the basis.
           * @return Callable value at the mapped physical point.
           */
          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

          /**
           * @brief Evaluates the pulled-back callable into storage.
           * @tparam T Result storage type.
           * @param res Result storage.
           * @param r Reference coordinates at which to evaluate the basis.
           * @return Value returned by the wrapped callable.
           */
          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(res, p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      /**
       * @brief Pushforward of a function from the reference polytope.
       *
       * Given a function @f$ \hat v @f$ defined on the reference polytope
       * @f$ \hat K @f$, this class defines the function @f$ v @f$ on the
       * physical polytope @f$ K @f$ by
       *
       * @f[
       *   v(x) = \hat v(\hat x),
       * @f]
       *
       * where @f$ \hat x @f$ are the reference coordinates of the physical
       * point @f$ x @f$.
       */
      template <class CallableType>
      class Pushforward : public FiniteElementSpacePushforwardBase<Pushforward<CallableType>>
      {
        public:
          /// @brief Callable type evaluated on reference coordinates.
          using FunctionType = CallableType;

          /**
           * @brief Constructs a pushforward from a reference-space callable.
           * @param v Callable evaluated on reference coordinates.
           */
          Pushforward(const FunctionType& v)
            : m_v(v)
          {}

          /**
           * @brief Copy constructor.
           * @param other Object to copy from.
           */
          Pushforward(const Pushforward& other) = default;

          /**
           * @brief Evaluates the pushed-forward callable.
           * @param p Physical point.
           * @return Callable value at the point reference coordinates.
           */
          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

          /**
           * @brief Evaluates the pushed-forward callable into storage.
           * @tparam T Result storage type.
           * @param res Result storage.
           * @param p Physical point.
           * @return Value returned by the wrapped callable.
           */
          template <class T>
          auto operator()(T& res, const Geometry::Point& p) const
          {
            return m_v(res, p.getReferenceCoordinates());
          }

          /**
           * @brief Returns the wrapped reference-space callable.
           * @return Callable reference.
           */
          constexpr
          const FunctionType& getFunction() const
          {
            return m_v;
          }

        private:
          const FunctionType m_v;
      };

      /**
       * @brief Constructs the distributed P0 space on the given mesh.
       *
       * Builds the local shard P0 space, then assigns contiguous global
       * DOF indices to owned cells. Ghost cells obtain their global DOF
       * index via a direct owner-push exchange:
       *
       * - Each rank iterates its owned cells and assigns a global DOF index
       *   @f$ g = \text{offset} + k @f$ where @f$ k @f$ is the local owned-
       *   cell counter.
       * - For each owned cell that also appears on neighboring ranks as a ghost
       *   (recorded in @c shard.getHalo(D)), the owner sends
       *   @c (globalCellID, globalDOF) to all holder ranks.
       * - Each non-owner rank posts a matching @c irecv from its owner
       *   (recorded in @c shard.getOwner(D)) and installs the received global
       *   DOF index in the local-to-global map.
       *
       * @note Collective over the mesh communicator, including empty ranks.
       * @param[in] mesh Distributed mesh on which the space is defined.
       */
      P0(const MeshType& mesh)
        : P0(mesh, 1)
      {}

      /**
       * @brief Collectively constructs one constant DOF per cell and component.
       *
       * All ranks in the mesh communicator participate, including empty ranks.
       * Ownership numbering and global size are established for this fixed
       * layout; subsequent size queries and copies require no communication.
       * @param[in] mesh Mesh on which the object is defined.
       * @param[in] vdim Number of components in the value range.
       */
      P0(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_fes(makeShardFES(mesh, vdim))
      {
        const auto& ctx   = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();
        const auto& shard = mesh.getShard();

        const size_t D = mesh.getDimension();

        // halo(D):  owned local cell -> set of remote ranks that have it as ghost
        const auto& halo  = shard.getHalo(D);
        // owner(D): ghost/shared local cell -> owner rank
        const auto& owner = shard.getOwner(D);

        const int rank = comm.rank();

        // Count owned local cells.
        m_owned = 0;
        const size_t localCellCount = shard.getPolytopeCount(D);
        for (size_t i = 0; i < localCellCount; ++i)
        {
          if (shard.isOwned(D, i))
            m_owned += vdim;
        }

        // Assign contiguous global DOF range for owned cells via prefix scan.
        const size_t inclusive = boost::mpi::scan(comm, m_owned, std::plus<size_t>());
        m_offset = inclusive - m_owned;
        boost::mpi::all_reduce(comm, m_owned, m_globalSize, std::plus<size_t>());

        // The local DOFs of cell i occupy [i * vdim, (i + 1) * vdim).
        // Pre-allocate the left map with an invalid sentinel.
        m_localToGlobal.left.assign(
          localCellCount * vdim, std::numeric_limits<Index>::max());

        // send[r]: messages to rank r — pairs (globalCellID, globalDOF).
        // Keyed by rank so every neighbor gets an entry (possibly empty).
        UnorderedMap<int, std::vector<std::pair<Index, Index>>> send;

        // Build the symmetric neighbor set from halo ∪ owner.
        // Both sides of every partition interface appear, so every isend is
        // matched by a corresponding irecv and no messages are orphaned.
        std::vector<int> neighbors;
        {
          UnorderedSet<int> nbrs;
          for (const auto& [i, peers] : halo)
          {
            for (const Index r : peers)
            {
              if (static_cast<int>(r) != rank)
                nbrs.insert(static_cast<int>(r));
            }
          }
          for (const auto& entry : owner)
          {
            if (static_cast<int>(entry.second) != rank)
              nbrs.insert(static_cast<int>(entry.second));
          }
          neighbors.assign(nbrs.begin(), nbrs.end());
        }
        for (int r : neighbors)
          send[r]; // default-construct empty vector for every neighbor

        // Assign global DOFs to owned cells and enqueue ghost notifications.
        Index dofIdx = 0;
        for (size_t i = 0; i < localCellCount; ++i)
        {
          if (!shard.isOwned(D, i))
            continue;

          const Index gid    = mesh.getGlobalIndex(D, i);
          const Index global = m_offset + dofIdx;
          for (size_t c = 0; c < vdim; ++c)
          {
            const Index local = static_cast<Index>(i * vdim + c);
            const Index globalComponent = global + c;
            m_localToGlobal.left[local] = globalComponent;
          }
          dofIdx += vdim;

          // Notify all neighbors that have this cell as a ghost.
          auto hit = halo.find(i);
          if (hit != halo.end())
          {
            for (const Index peer : hit->second)
            {
              const int rpeer = static_cast<int>(peer);
              if (rpeer == rank)
                continue;
              send[rpeer].push_back({ gid, global });
            }
          }
        }
        assert(dofIdx == static_cast<Index>(m_owned));

        // Symmetric exchange: every neighbor always sends and receives,
        // even if the payload is empty.  Tag 50 is reserved for P0 DOF
        // exchange and does not overlap with other FES or SubMesh tags.
        static constexpr int kTagP0Dof = 50;

        UnorderedMap<int, std::vector<std::pair<Index, Index>>> recv;
        std::vector<boost::mpi::request> reqs;
        reqs.reserve(2 * neighbors.size());

        for (int r : neighbors)
        {
          recv[r]; // default-construct
          reqs.push_back(comm.irecv(r, kTagP0Dof, recv[r]));
          reqs.push_back(comm.isend(r, kTagP0Dof, send[r]));
        }

        boost::mpi::wait_all(reqs.begin(), reqs.end());

        // Install global DOFs for ghost cells received from owners.
        for (auto& [r, msgs] : recv)
        {
          for (const auto& [gid, global] : msgs)
          {
            const auto liOpt = mesh.getLocalIndex(D, gid);
            // The owner's halo is derived from the parent shard, which may
            // list ranks that have this entity only in ghost cells.  When the
            // mesh is a SubMesh that excludes those ghost cells the entity is
            // absent here, so simply skip the entry.
            if (!liOpt)
              continue;
            const Index li = *liOpt;
            assert(!shard.isOwned(D, li));

            for (size_t c = 0; c < vdim; ++c)
            {
              const Index local = li * vdim + c;
              const Index globalComponent = global + c;
              m_localToGlobal.left[local] = globalComponent;
            }
          }
        }

#ifndef NDEBUG
        for (size_t i = 0; i < localCellCount * vdim; ++i)
          assert(m_localToGlobal.left[i] != std::numeric_limits<Index>::max());
#endif

        // Ghost numbering arrives in neighbor order, not global-index order.
        // Sorting once avoids repeatedly shifting the flat map's storage.
        std::vector<std::pair<Index, Index>> globalLocalPairs;
        globalLocalPairs.reserve(m_localToGlobal.left.size());
        for (Index local = 0; local < m_localToGlobal.left.size(); ++local)
          globalLocalPairs.emplace_back(m_localToGlobal.left[local], local);
        std::sort(globalLocalPairs.begin(), globalLocalPairs.end());
        m_localToGlobal.right.reserve(globalLocalPairs.size());
        for (const auto& [global, local] : globalLocalPairs)
          m_localToGlobal.right.emplace_hint(m_localToGlobal.right.end(), global, local);
        assert(m_localToGlobal.right.size() == globalLocalPairs.size());
      }

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      P0(const P0& other) = default;

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      P0(P0&& other) = default;

      /**
       * @brief Move assignment operator.
       * @returns Reference to this object after the operation.
       * @param other Object to move from.
       */
      P0& operator=(P0&& other) = default;

      /**
       * @brief Returns the local shard finite element space.
       *
       * The distributed P0 space is constructed from the P0 space defined on
       * the local mesh shard. This method provides access to that underlying
       * local finite element space.
       * @returns The local shard finite element space.
       */
      const FESType& getShard() const
      {
        return m_fes;
      }

      /**
       * @brief Returns the global ownership range of this rank.
       *
       * The interval @f$ [\text{begin}, \text{end}) @f$ identifies the global
       * degrees of freedom owned by the current MPI rank.
       *
       * @param[out] begin First global DOF owned by this rank.
       * @param[out] end   One past the last global DOF owned by this rank.
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        begin = m_offset;
        end   = m_offset + m_owned;
      }

      /**
       * @brief Returns the global distributed DOF index for a local shard DOF.
       *
       * @param[in] localIdx Local shard degree-of-freedom index (= local cell index).
       * @return Global distributed DOF index.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        return m_localToGlobal.left.at(localIdx);
      }

      /**
       * @brief Returns the local shard DOF index for a global distributed DOF.
       *
       * If the global DOF is present on this rank (owned or ghost), its local
       * shard index (= local cell index) is returned.
       *
       * @param[in] globalIdx Global distributed DOF index.
       * @return Local shard index if present, otherwise @c std::nullopt.
       */
      Optional<Index> getLocalIndex(Index globalIdx) const
      {
        auto it = m_localToGlobal.right.find(globalIdx);
        if (it == m_localToGlobal.right.end())
          return std::nullopt;
        return it->second;
      }

      /**
       * @brief Returns the finite element attached to local polytope @f$(d, i)@f$.
       *
       * @param[in] d Topological dimension of the polytope.
       * @param[in] i Local index of the polytope within the shard.
       * @return Finite element bound to the local polytope.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_fes.getFiniteElement(d, i);
      }

      /**
       * @brief Returns the global number of degrees of freedom.
       *
       * For a P0 space, degrees of freedom are associated with mesh cells.
       * The global dimension of the space equals the total number of cells
       * in the distributed mesh:
       *
       * @f[
       *   \dim(V_h) = N_c \, k
       * @f]
       *
       * where @f$ N_c @f$ is the global number of mesh cells and @f$k@f$
       * is the number of field components.
       *
       * @note Non-collective: returns the size established during construction.
       * @return Global number of degrees of freedom.
       */
      size_t getSize() const override
      {
        return m_globalSize;
      }

      /**
       * @brief Returns the vector dimension (1 for scalar P0).
       *
       * @return Vector dimension of the space.
       */
      size_t getVectorDimension() const override
      {
        return m_fes.getVectorDimension();
      }

      /**
       * @brief Returns the distributed mesh on which the finite element space
       * is defined.
       *
       * @return Reference to the distributed mesh.
       */
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Returns global DOF indices attached to local polytope @f$(d, i)@f$.
       *
       * The local shard DOF list (= @f$\{i\}@f$ for P0) is retrieved and
       * mapped to the distributed global DOF index.
       *
       * @param[in] d Topological dimension of the polytope (must equal mesh dimension).
       * @param[in] i Local polytope index in the shard.
       * @return Thread-local array containing the single global DOF index.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        static thread_local IndexArray s_dofs;
        assert(i < static_cast<Index>(getMesh().getShard().getPolytopeCount(d)));

        s_dofs = m_fes.getDOFs(d, i);
        for (auto& dof : s_dofs)
          dof = this->getGlobalIndex(dof);

        return s_dofs;
      }

      /**
       * @brief Returns the global DOF index for local basis function @p localDof
       * on polytope @f$(d, i)@f$.
       *
       * @param[in] p       Pair @f$(d, i)@f$ identifying a local shard polytope.
       * @param[in] localDof Local basis-function index on the polytope (always 0 for P0).
       * @return Global distributed DOF index.
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index localDof) const override
      {
        const auto& [d, i] = p;
        assert(i < static_cast<Index>(getMesh().getShard().getPolytopeCount(d)));

        const Index local = m_fes.getGlobalIndex({ d, i }, localDof);
        return getGlobalIndex(local);
      }

      /**
       * @brief Returns a pullback wrapper on local polytope @f$(d, i)@f$.
       *
       * @tparam Callable Callable type (FunctionBase or plain callable).
       * @param[in] p Local polytope pair @f$(d, i)@f$.
       * @param[in] v Function defined on physical coordinates.
       * @return Pullback wrapper mapping reference points to physical evaluation.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& p, Callable&& v) const
      {
        const auto& [d, i] = p;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /**
       * @brief Returns a pushforward wrapper on local polytope @f$(d, i)@f$.
       *
       * @tparam CallableType Callable defined on reference coordinates.
       * @param[in] v Reference-space callable.
       * @return Pushforward wrapper mapping physical points to reference evaluation.
       * @param entity Mesh entity; this pushforward does not require entity-dependent data.
       */
      template <class CallableType>
      auto getPushforward(
        [[maybe_unused]] const std::pair<size_t, Index>& entity, CallableType&& v) const
      {
        return Pushforward<CallableType>(std::forward<CallableType>(v));
      }

      /**
       * @brief Returns a pushforward wrapper for an explicit polytope object.
       *
       * @tparam CallableType Callable defined on reference coordinates.
       * @param[in] v Reference-space callable.
       * @return Pushforward wrapper mapping physical points to reference evaluation.
       * @param entity Mesh entity; this pushforward does not require entity-dependent data.
       */
      template <class CallableType>
      auto getPushforward(
        [[maybe_unused]] const Geometry::Polytope& entity, CallableType&& v) const
      {
        return Pushforward<CallableType>(std::forward<CallableType>(v));
      }

    private:
      /**
       * @brief Creates the local finite element space used by a mesh shard.
       * @param mesh Mesh whose entities are traversed.
       * @param vdim Number of components in the value range.
       * @returns Local finite element space with the requested vector dimension.
       */
      static FESType makeShardFES(const MeshType& mesh, size_t vdim)
      {
        if constexpr (std::is_same_v<Range, Real> || std::is_same_v<Range, Complex>)
        {
          assert(vdim == 1);
          return FESType(mesh.getShard());
        }
        else
        {
          return FESType(mesh.getShard(), vdim);
        }
      }

      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;

      size_t m_offset;
      size_t m_owned;
      size_t m_globalSize;
      IndexBimap m_localToGlobal;
  };
}

namespace Rodin::MPI
{
  /// @brief Convenience alias for the default distributed scalar P0 space.
  using P0 = Variational::P0<Real, Geometry::Mesh<Context::MPI>>;
}

namespace Rodin::Variational
{
  /// @brief Matrix-range finite element or expression specialization.
  template <class Scalar>
  class P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::MPI>,
        P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::MPI>>>
  {
    public:
      /// @brief Local matrix space on the mesh shard.
      using FESType = P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>>;
      /// @brief Mesh and execution context of this specialization.
      using MeshType = Geometry::Mesh<Context::MPI>;
      /// @brief Scalar space supplying this family's DOF numbering.
      using ScalarSpace = P0<Scalar, Geometry::Mesh<Context::MPI>>;
      /// @brief Matrix reference element for this family.
      using ElementType = P0Element<Math::SpatialMatrix<Scalar>>;
      /// @brief Existing finite-element-space interface.
      using Parent = FiniteElementSpace<MeshType,
        P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::MPI>>>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename ScalarSpace::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Local or distributed execution context.
      using ContextType = typename ScalarSpace::ContextType;
      /**
       * @brief Expands scalar DOF maps into interleaved row-major matrix components.
       * @param mesh Mesh on which the object is defined.
       * @param rows Number of rows.
       * @param cols Number of columns.
       */
      P0(const MeshType& mesh, size_t rows, size_t cols)
        : m_scalar(ScalarSpace(mesh)),
          m_rows(rows),
          m_cols(cols),
          m_shard(mesh.getShard(), rows, cols),
          m_size(0)
      {
        if (rows == 0 || cols == 0 || rows > RODIN_MAXIMAL_SPACE_DIMENSION ||
          cols > RODIN_MAXIMAL_SPACE_DIMENSION)
          Alert::Exception() << "SpatialMatrix ranges require 1 to 3 rows and columns."
                             << Alert::Raise;
        // Distributed scalar spaces may compute their size collectively.
        // Resolve it at construction, never from a local evaluation loop.
        m_size = m_scalar.getSize() * rows * cols;
        m_dofs.resize(mesh.getDimension() + 1);
        for (size_t d = mesh.getDimension(); d <= mesh.getDimension(); ++d)
        {
          const size_t count = mesh.getConnectivity().getCount(d);
          m_dofs[d].reserve(count);
          for (size_t i = 0; i < count; ++i)
          {
            const auto& scalarDOFs = m_scalar.getDOFs(d, i);
            auto& dofs = m_dofs[d].emplace_back(scalarDOFs.size() * rows * cols);
            for (size_t a = 0; a < static_cast<size_t>(scalarDOFs.size()); ++a)
            {
              for (size_t c = 0; c < rows * cols; ++c)
                dofs[a * rows * cols + c] = scalarDOFs[a] * rows * cols + c;
            }
            const auto& scalarFE = m_scalar.getFiniteElement(d, i);
            m_elements.try_emplace(scalarFE.getGeometry(), scalarFE, rows, cols);
          }
        }
        for (size_t i = 0; i < m_shard.getSize(); ++i)
          m_globalToLocal.emplace(getGlobalIndex(i), i);
      }

      /**
       * @brief Copies the space and its DOF maps.
       * @param other Object to copy from.
       */
      P0(const P0& other) = default;
      /**
       * @brief Moves the space and its DOF maps.
       * @param other Object to move from.
       */
      P0(P0&& other) = default;
      /**
       * @brief Copies the space and its DOF maps.
       * @returns Reference to this object after the operation.
       * @param other Object to copy from.
       */
      P0& operator=(const P0& other) = default;
      /**
       * @brief Moves the space and its DOF maps.
       * @returns Reference to this object after the operation.
       * @param other Object to move from.
       */
      P0& operator=(P0&& other) = default;

      size_t getSize() const override
      {
        return m_size;
      }
      size_t getVectorDimension() const override
      {
        return m_rows * m_cols;
      }
      /**
       * @brief Returns the number of matrix rows.
       * @returns The number of matrix rows.
       */
      size_t getRows() const
      {
        return m_rows;
      }
      /**
       * @brief Returns the number of matrix columns.
       * @returns The number of matrix columns.
       */
      size_t getColumns() const
      {
        return m_cols;
      }
      const MeshType& getMesh() const override
      {
        return m_scalar.getMesh();
      }
      /**
       * @brief Returns the scalar space supplying topology and component-independent maps.
       * @returns The scalar space supplying topology and component-independent maps.
       */
      const ScalarSpace& getScalarSpace() const
      {
        return m_scalar;
      }

      /**
       * @brief Returns the matrix reference element of a mesh entity.
       * @param i Index of the requested entry.
       * @returns The matrix reference element of a mesh entity.
       * @param d Topological dimension of the entity.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_elements.at(m_scalar.getFiniteElement(d, i).getGeometry());
      }

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_dofs.at(d).at(i);
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index local) const override
      {
        const Index components = m_rows * m_cols;
        return m_scalar.getGlobalIndex(p, local / components) * components +
          local % components;
      }

      /**
       * @brief Pulls a physical callable back to a reference element.
       * @param value Value to store or assign.
       * @returns The pullback.
       * @param p Entity identified by its dimension and local index.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& p, Callable&& value) const
      {
        return m_scalar.getPullback(p, std::forward<Callable>(value));
      }

      /**
       * @brief Pushes a reference callable forward to the physical mesh.
       * @param value Value to store or assign.
       * @returns The pushforward.
       * @param p Entity identified by its dimension and local index.
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& p, Callable&& value) const
      {
        return m_scalar.getPushforward(p, std::forward<Callable>(value));
      }

      /**
       * @brief Returns the matrix space on the local mesh shard.
       * @returns The matrix space on the local mesh shard.
       */
      const FESType& getShard() const
      {
        return m_shard;
      }

      /**
       * @brief Returns the half-open range of owned global component DOFs.
       * @param begin Iterator to the first entry.
       * @param end Iterator past the last entry.
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        this->getScalarSpace().getOwnershipRange(begin, end);
        begin *= this->getVectorDimension();
        end *= this->getVectorDimension();
      }

      /**
       * @brief Maps a local component DOF to its global coefficient index.
       * @param local Index in the local numbering.
       * @returns The global index.
       */
      Index getGlobalIndex(Index local) const
      {
        assert(static_cast<size_t>(local) < m_shard.getSize());
        const Index components = this->getVectorDimension();
        return this->getScalarSpace().getGlobalIndex(local / components) * components +
          local % components;
      }

      /**
       * @brief Returns the shard-local index of a global DOF, if present.
       * @param global Index in the global numbering.
       * @returns The shard-local index of a global DOF, if present.
       */
      Optional<Index> getLocalIndex(Index global) const
      {
        const auto it = m_globalToLocal.find(global);
        if (it == m_globalToLocal.end())
          return std::nullopt;
        return it->second;
      }

    private:
      ScalarSpace m_scalar;
      size_t m_rows, m_cols;
      FESType m_shard;
      std::map<Index, Index> m_globalToLocal;
      size_t m_size;
      std::vector<std::vector<IndexArray>> m_dofs;
      std::map<Geometry::Polytope::Type, ElementType> m_elements;
  };
}

#endif
