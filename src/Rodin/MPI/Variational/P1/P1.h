#ifndef RODIN_MPI_VARIATIONAL_P1_P1_H
#define RODIN_MPI_VARIATIONAL_P1_P1_H

#include <mpi.h>
#include <sys/mman.h>

#include <boost/serialization/version.hpp>
#include <boost/serialization/split_free.hpp>

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
       * @brief Mapping for the scalar/complex P1 space.
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
       * @brief Inverse mapping for the scalar/complex P1 space.
       */
      template <class CallableType>
      class Pushforward : public FiniteElementSpacePushforwardBase<Pushforward<CallableType>>
      {
        public:
          using FunctionType = CallableType;

          /**
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
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

      P1(const MeshType& mesh)
        : m_mesh(mesh),
          m_fes(mesh.getShard())
      {
        const auto& ctx = mesh.getContext();
        const auto& comm  = ctx.getCommunicator();
        const auto& shard = mesh.getShard();
        const auto& halo = shard.getHalo(0);
        const auto& owner = shard.getOwner(0);

        m_owned = halo.size();
        const size_t inclusive = boost::mpi::scan(comm, m_owned, std::plus<size_t>());
        m_offset = inclusive - m_owned;

        FlatMap<Index, std::vector<std::pair<Index, Index>>> push, pull;
        Index dofIdx = 0;
        for (size_t i = 0; i < shard.getVertexCount(); ++i)
        {
          if (shard.isOwned(0, i))
          {
            const Index id = mesh.getGlobalIndex(0, i);
            const auto& dofs = m_fes.getDOFs(0, i);
            assert(dofs.size() == 1);
            const Index local = dofs[0];
            const Index global = dofIdx + m_offset;
            const auto [it, inserted] = m_local_to_global.right.emplace(global, local);
            assert(inserted);
            dofIdx++;

            for (const Index& peer : halo.at(i))
            {
              assert(comm.rank() >= 0);
              assert(peer != static_cast<Index>(comm.rank()));
              push[peer].push_back({ id, global });
            }
          }
          else
          {
            pull.try_emplace(owner.at(i));
          }
        }

        assert(m_owned == dofIdx);

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
            assert(dofs.size() == 1);
            const auto [it, inserted] = m_local_to_global.right.emplace(global, dofs[0]);
            assert(inserted);
          }
        }

        m_local_to_global.left.resize(m_local_to_global.right.size());
        for (const auto& [global, local] : m_local_to_global.right)
          m_local_to_global.left[local] = global;
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

        m_owned = halo.size();
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

        assert(m_owned == dofIdx);

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
      }

      P1(const P1& other) = default;

      P1(P1&& other) = default;

      P1& operator=(P1&& other) = default;

      const FESType& getShard() const
      {
        return m_fes;
      }

      void getOwnershipRange(Index& begin, Index& end) const
      {
        begin = m_offset;
        end = m_offset + m_owned;
      }

      /**
       * @brief Returns the global distributed index of the given local
       * distributed index.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        return m_local_to_global.left.at(localIdx);
      }

      /**
       * @brief Returns the local distributed index of the given global
       * distributed index.
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
       * @brief Returns the global finite element for the given global index.
       * @param[in] d Dimension of the element
       * @param[in] globalIdx Global distributed index of the polytope in the
       * distributed mesh
       */
      const auto& getFiniteElement(size_t d, Index globalIdx) const
      {
        static thread_local ElementType s_element;
        const auto& mesh = getMesh();
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        const auto& shard = mesh.getShard();
        const auto localIdx = mesh.getLocalIndex(d, globalIdx);
        Optional<ElementType> send;
        if (localIdx)
        {
          const Index owner = shard.getOwner(d).at(*localIdx);
          if (owner == comm.rank())
            send = m_fes.getFiniteElement(*localIdx, globalIdx);
        }
        auto recv =
          boost::mpi::all_reduce(
              comm, send, [](const auto& a, const auto& b) { return a ? a : b; });
        assert(recv);
        s_element = *recv;
        return s_element;
      }

      size_t getSize() const override
      {
        const auto& mesh = getMesh();
        return mesh.getVertexCount() * getVectorDimension();
      }

      size_t getVectorDimension() const override
      {
        const auto& fes = this->getShard();
        return fes.getVectorDimension();
      }

      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      const IndexArray& getDOFs(size_t d, Index globalIdx) const override
      {
        static thread_local IndexArray s_dofs;
        const auto& mesh = getMesh();
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        const auto& shard = mesh.getShard();
        const auto localIdx = mesh.getLocalIndex(d, globalIdx);
        Optional<IndexArray> send;
        if (localIdx)
        {
          const Index& owner = shard.getOwner(d).at(*localIdx);
          assert(comm.rank() >= 0);
          if (owner == static_cast<Index>(comm.rank()))
            send = m_fes.getDOFs(*localIdx, globalIdx);
        }
        auto recv =
          boost::mpi::all_reduce(
              comm, send, [](const auto& a, const auto& b) { return a ? a : b; });
        assert(recv);
        s_dofs = *recv;
        for (auto& dof : s_dofs)
          dof = getGlobalIndex(dof);
        return s_dofs;
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index localDof) const override
      {
        const auto& [d, globalIdx] = p;
        const auto& mesh = getMesh();
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        const auto& shard = mesh.getShard();
        const auto& fes = getShard();
        const auto localIdx = mesh.getLocalIndex(d, globalIdx);
        Optional<Index> send;
        if (localIdx)
        {
          const Index& owner = shard.getOwner(d).at(*localIdx);
          assert(comm.rank() >= 0);
          if (owner == static_cast<Index>(comm.rank()))
            send = fes.getGlobalIndex({ d, *localIdx }, localDof);
        }
        auto recv =
          boost::mpi::all_reduce(
              comm, send, [](const auto& a, const auto& b) { return a ? a : b; });
        assert(recv);
        return *recv;
      }

      template <class FunctionDerived>
      auto getPullback(const std::pair<size_t, Index>& p, const FunctionBase<FunctionDerived>& v) const
      {
        const auto& [d, globalIdx] = p;
        const auto& mesh = getMesh();
        return Pullback<FunctionDerived>(*mesh.getPolytope(d, globalIdx), v);
      }

      template <class CallableType>
      auto getPushforward(const std::pair<size_t, Index>& idx, const CallableType& v) const
      {
        return typename FESType::template Pushforward<CallableType>(v);
      }

      template <class CallableType>
      auto getPushforward(const Geometry::Polytope& polytope, const CallableType& v) const
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
