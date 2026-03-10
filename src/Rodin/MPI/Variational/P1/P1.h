#ifndef RODIN_MPI_VARIATIONAL_P1_P1_H
#define RODIN_MPI_VARIATIONAL_P1_P1_H

#include <mpi.h>
#include <sys/mman.h>

#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>
#include "Rodin/Serialization/Optional.h"

#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Variational/FiniteElementSpace.h"

#include "Rodin/Variational/P1/P1.h"

#include "Rodin/Array.h"

namespace Rodin::Variational
{
  template <class Range>
  class P1<Range, Geometry::Mesh<Context::MPI>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>, P1<Range, Geometry::Mesh<Context::MPI>>>
  {
    public:
      struct IndexBimap
      {
        std::vector<Index> left;
        FlatMap<Index, Index> right;
      };

      /// Represents the Context of the P1 space
      using ContextType = Context::MPI;

      /// Type of the local finite element space
      using FESType = P1<Range, Geometry::Mesh<Context::Local>>;

      using ScalarType = typename FESType::ScalarType;

      /// Range type of value
      using RangeType = typename FESType::RangeType;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

      using Parent::getGlobalIndex;

      /**
       * @brief Pullback of a function to the reference polytope.
       *
       * Given a function @f$ v @f$ defined on the physical polytope @f$ K @f$, this
       * class defines the function @f$ \hat v @f$ on the reference polytope
       * @f$ \hat K @f$ by
       *
       * @f[
       *   \hat v(\hat x) = v(x(\hat x)),
       * @f]
       *
       * where @f$ x(\hat x) @f$ is the point of @f$ K @f$ with reference coordinates
       * @f$ \hat x @f$.
       */
      template <class FunctionDerived>
      class Pullback : public FiniteElementSpacePullbackBase<Pullback<FunctionDerived>>
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
       *
       * Given a function @f$ \hat v @f$ defined on the reference polytope
       * @f$ \hat K @f$, this class defines the function @f$ v @f$ on the
       * physical polytope @f$ K @f$ by
       *
       * @f[
       *   v(x) = \hat v(\hat x)
       * @f]
       *
       * where @f$ \hat x @f$ are the reference coordinates of the physical
       * point @f$ x @f$.
       */
      template <class CallableType>
      class Pushforward : public FiniteElementSpacePushforwardBase<Pushforward<CallableType>>
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
       * @brief Constructs a distributed scalar P1 finite element space on the given mesh.
       *
       * This constructor builds the local shard P1 space and assigns a global
       * distributed numbering to its degrees of freedom. For scalar P1 elements,
       * degrees of freedom are associated with mesh vertices.
       *
       * Each rank numbers its owned vertex degrees of freedom contiguously, then
       * exchanges this numbering with neighboring ranks so that ghost vertices are
       * assigned the global indices of their owners. The resulting maps store the
       * correspondence between local shard degrees of freedom and global distributed
       * degrees of freedom.
       *
       * @param[in] mesh Distributed mesh on which the space is defined.
       */
      P1(const MeshType& mesh)
        : m_mesh(mesh),
          m_fes(mesh.getShard())
      {
        const auto& ctx   = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();
        const auto& shard = mesh.getShard();

        // halo: owned local vertex -> peers that need it
        const auto& halo  = shard.getHalo(0);
        // owner: ghost local vertex -> owning rank
        const auto& owner = shard.getOwner(0);

        const int P    = comm.size();
        const int rank = comm.rank();

        // Owned vertices are exactly the keys in halo(0) (by construction).
        m_owned = halo.size();

        const size_t inclusive = boost::mpi::scan(comm, m_owned, std::plus<size_t>());
        m_offset = inclusive - m_owned;

        // 1) Use rank-indexed buffers (dense ranks), no FlatMap in the hot path.
        std::vector<std::vector<std::pair<Index, Index>>> send(P); // to peer: (globalVertexId, globalDof)
        std::vector<std::vector<std::pair<Index, Index>>> recv(P); // from owner: (globalVertexId, globalDof)
        std::vector<char> need_recv(P, 0);

        // 2) Collect (globalDof -> localDof) as a flat vector, then bulk-build FlatMap once.
        //    Estimate: owned + ghosts (owner map size is a good proxy for ghosts in this shard).
        std::vector<std::pair<Index, Index>> gl_pairs;
        gl_pairs.reserve(m_owned + owner.size());

        // ----- Owned part: iterate halo entries (do NOT scan all vertices) -----
        Index dofIdx = 0;
        for (const auto& [lv, peers] : halo)
        {
          // lv is local vertex index (owned)
          const Index gid   = mesh.getGlobalIndex(0, lv);   // global vertex id
          const Index local = m_fes.getDOFs(0, lv)[0];      // scalar P1 => 1 dof (often == lv)
          const Index global = m_offset + dofIdx++;

          gl_pairs.push_back({global, local});

          for (const Index& peer : peers)
          {
            const int rpeer = static_cast<int>(peer);
            if (rpeer == rank) continue;
            send[rpeer].push_back({gid, global});
          }
        }
        assert(dofIdx == static_cast<Index>(m_owned));

        // ----- Ghost part: we only need to post receives to owners that exist in this shard -----
        for (const auto& [lv, own] : owner)
        {
          const int ro = static_cast<int>(own);
          if (ro != rank)
            need_recv[ro] = 1;
        }

        // 3) Post receives and sends.
        std::vector<boost::mpi::request> reqs;
        reqs.reserve(static_cast<size_t>(P) * 2);

        for (int r = 0; r < P; ++r)
          if (need_recv[r])
            reqs.push_back(comm.irecv(r, 0, recv[r]));

        for (int r = 0; r < P; ++r)
          if (!send[r].empty())
            reqs.push_back(comm.isend(r, 0, send[r]));

        boost::mpi::wait_all(reqs.begin(), reqs.end());

        // 4) Consume received (gid, globalDof) pairs and map them to local DOFs.
        for (int r = 0; r < P; ++r)
        {
          for (const auto& [gid, global] : recv[r])
          {
            const auto lvOpt = mesh.getLocalIndex(0, gid);
            assert(lvOpt);
            const Index lv = *lvOpt;

            const Index local = m_fes.getDOFs(0, lv)[0];
            gl_pairs.push_back({global, local});
          }
        }

        // 5) Bulk-build the FlatMap ONCE (avoid O(n^2) incremental inserts).
        std::sort(gl_pairs.begin(), gl_pairs.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        gl_pairs.erase(std::unique(gl_pairs.begin(), gl_pairs.end(),
                                   [](const auto& a, const auto& b) { return a.first == b.first; }),
                       gl_pairs.end());

        m_local_to_global.right = FlatMap<Index, Index>(gl_pairs.begin(), gl_pairs.end());

        // 6) Build left in one linear pass.
        // For scalar P1 on the local shard, local DOF count should be shard vertex count.
        const size_t localDofCount = shard.getVertexCount();
        m_local_to_global.left.assign(localDofCount, std::numeric_limits<Index>::max());

        for (const auto& [global, local] : m_local_to_global.right)
        {
          assert(local < localDofCount);
          m_local_to_global.left[local] = global;
        }
      }

      P1(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_fes(mesh.getShard(), vdim)
      {
        static thread_local std::vector<Index> s_send;

        const auto& ctx = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();
        const auto& shard = mesh.getShard();
        const auto& halo = shard.getHalo(0);
        const auto& owner = shard.getOwner(0);

        const size_t nv = halo.size();
        m_owned = nv * vdim;
        const size_t inclusive = boost::mpi::scan(comm, m_owned, std::plus<size_t>());
        m_offset = inclusive - m_owned;

        FlatMap<Index, std::vector<std::pair<Index, std::vector<Index>>>> push, pull;
        Index dofIdx = 0;
        for (size_t i = 0; i < shard.getVertexCount(); ++i)
        {
          if (shard.isOwned(0, i))
          {
            const Index id = mesh.getGlobalIndex(0, i);
            const auto& dofs = m_fes.getDOFs(0, i);

            s_send.clear();
            for (const Index& local : dofs)
            {
              const Index global = dofIdx + m_offset;
              s_send.push_back(global);

              const auto [it, inserted] = m_local_to_global.right.emplace(global, local);
              assert(inserted);

              dofIdx++;
            }

            for (const Index& peer : halo.at(i))
            {
              assert(comm.rank() >= 0);
              assert(peer != static_cast<Index>(comm.rank()));
              push[peer].push_back({ id, s_send });
            }
          }
          else
          {
            pull.try_emplace(owner.at(i));
          }
        }

        std::vector<boost::mpi::request> irecv;
        for (auto& [owner, requested] : pull)
          irecv.push_back(comm.irecv(owner, 0, pull[owner]));

        std::vector<boost::mpi::request> isend;
        for (const auto& [peer, requested] : push)
          isend.push_back(comm.isend(peer, 0, push[peer]));

        boost::mpi::wait_all(isend.begin(), isend.end());

        boost::mpi::wait_all(irecv.begin(), irecv.end());

        for (const auto& [owner, requested] : pull)
        {
          for (const auto& [id, global] : requested)
          {
            const auto i = mesh.getLocalIndex(0, id);
            assert(i);
            const auto& dofs = m_fes.getDOFs(0, *i);
            assert(dofs.size() >= 0);
            assert(dofs.size() == global.size());
            for (size_t k = 0; k < global.size(); k++)
            {
              const auto [it, inserted] = m_local_to_global.right.emplace(global[k], dofs[k]);
              assert(inserted);
            }
          }
        }

        const size_t localDofCount = m_fes.getSize();
        m_local_to_global.left.assign(localDofCount, 0);

        for (const auto& [global, local] : m_local_to_global.right)
        {
          assert(local < localDofCount);
          m_local_to_global.left[local] = global;
        }
      }

      P1(const P1& other) = default;

      P1(P1&& other) = default;

      P1& operator=(P1&& other) = default;

      /**
       * @brief Returns the local shard finite element space.
       *
       * The distributed P1 space is constructed from the P1 space defined on the
       * local mesh shard. This method provides access to that underlying local
       * finite element space.
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
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        begin = m_offset;
        end = m_offset + m_owned;
      }

      /**
       * @brief Returns the global distributed index corresponding to a local
       * shard degree of freedom.
       *
       * @param[in] localIdx Local shard degree of freedom index.
       * @return Global distributed degree of freedom index.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        return m_local_to_global.left.at(localIdx);
      }

      /**
       * @brief Returns the local shard index of a global distributed degree of freedom.
       *
       * If the global degree of freedom is present on this rank (owned or ghost),
       * its local shard index is returned.
       *
       * @param[in] globalIdx Global distributed degree of freedom index.
       * @return Local shard index if present, otherwise `std::nullopt`.
       */
      Optional<Index> getLocalIndex(Index globalIdx) const
      {
        auto find = m_local_to_global.right.find(globalIdx);
        if (find == m_local_to_global.right.end())
          return std::nullopt;
        else
          return *find;
      }

      /**
       * @brief Returns the finite element associated with the local shard polytope @f$(d,i)@f$.
       *
       * The key @f$(d,i)@f$ refers to a polytope in the local mesh shard of this MPI rank.
       * Finite elements are therefore retrieved directly from the underlying shard
       * finite element space.
       *
       * @param[in] d Topological dimension of the polytope.
       * @param[in] i Local index of the polytope within the shard.
       *
       * @pre @f$0 \le i < |\mathcal{T}_d^{(loc)}|@f$.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        assert(i < this->getMesh().getShard();.getPolytopeCount(d));
        return m_fes.getFiniteElement(d, i);
      }

      /**
       * @brief Returns the global number of degrees of freedom of the distributed P1 space.
       *
       * For a P1 space, degrees of freedom are associated with mesh vertices. If the
       * vector dimension is @f$k@f$, the global dimension of the space is
       *
       * @f[
       *   \dim(V_h) = N_v \, k
       * @f]
       *
       * where @f$N_v@f$ is the global number of mesh vertices.
       *
       * @return Global number of degrees of freedom.
       */
      size_t getSize() const override
      {
        const auto& mesh = getMesh();
        return mesh.getVertexCount() * getVectorDimension();
      }

      /**
       * @brief Returns the vector dimension of the finite element space.
       *
       * For scalar P1 spaces this value is @f$ 1 @f$. For vector-valued spaces
       * it corresponds to the number of components carried by each degree of
       * freedom.
       *
       * @return Vector dimension of the space.
       */
      size_t getVectorDimension() const override
      {
        const auto& fes = this->getShard();
        return fes.getVectorDimension();
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

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        static thread_local IndexArray s_dofs;
        const auto& shard = getMesh().getShard();
        assert(i < shard.getPolytopeCount(d));

        s_dofs = m_fes.getDOFs(d, i);
        for (auto& dof : s_dofs)
          dof = this->getGlobalIndex(dof);

        return s_dofs;
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index localDof) const override
      {
        const auto& [d, i] = p;
        const auto& shard = getMesh().getShard();
        assert(i < shard.getPolytopeCount(d));

        const Index local = m_fes.getGlobalIndex({ d, i }, localDof);
        return getGlobalIndex(local);
      }

      template <class FunctionDerived>
      auto getPullback(const std::pair<size_t, Index>& p, const FunctionBase<FunctionDerived>& v) const
      {
        const auto& [d, i] = p;
        const auto& mesh = getMesh();
        return Pullback<FunctionDerived>(*mesh.getPolytope(d, i), v);
      }

      template <class CallableType>
      auto getPushforward(const std::pair<size_t, Index>&, const CallableType& v) const
      {
        return typename FESType::template Pushforward<CallableType>(v);
      }

      template <class CallableType>
      auto getPushforward(const Geometry::Polytope&, const CallableType& v) const
      {
        return typename FESType::Pushforward(v);
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;

      size_t m_offset;
      size_t m_owned;
      IndexBimap m_local_to_global;
  };
}

namespace Rodin::MPI
{
  using P1 = Variational::P1<Real, Geometry::Mesh<Context::MPI>>;
}

#endif
