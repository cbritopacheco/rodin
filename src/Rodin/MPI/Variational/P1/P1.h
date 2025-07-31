#ifndef RODIN_MPI_VARIATIONAL_P1_P1_H
#define RODIN_MPI_VARIATIONAL_P1_P1_H

#include <mpi.h>
#include <boost/serialization/optional.hpp>

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
      using Bimap =
        boost::bimap<
          boost::bimaps::vector_of<Index>,
          boost::bimaps::unordered_set_of<Index>>;

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
      class Mapping : public FiniteElementSpaceMappingBase<Mapping<FunctionDerived>>
      {
        public:
          using FunctionType = FunctionBase<FunctionDerived>;

          Mapping(const Geometry::Polytope& polytope, const FunctionType& v)
            : m_polytope(polytope), m_v(v.copy())
          {}

          Mapping(const Mapping&) = default;

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
      class InverseMapping : public FiniteElementSpaceInverseMappingBase<InverseMapping<CallableType>>
      {
        public:
          using FunctionType = CallableType;

          /**
           * @param[in] polytope Reference to polytope on the mesh.
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          InverseMapping(const FunctionType& v)
            : m_v(v)
          {}

          InverseMapping(const InverseMapping&) = default;

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
        const auto& comm  = mesh.getContext().getCommunicator();
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
            for (const Index& peer : halo.at(i))
            {
              assert(comm.rank() >= 0);
              assert(peer != static_cast<Index>(comm.rank()));
              const Index& global = mesh.getGlobalIndex(0, i);
              push[peer].push_back({ global, dofIdx + m_offset });
            }
            m_loc2glob.right.insert({ dofIdx + m_offset, i });
            dofIdx++;
          }
          else
          {
            pull.try_emplace(owner.at(i));
          }
        }

        assert(m_owned == dofIdx);

        std::vector<boost::mpi::request> requests;
        for (auto& [owner, requested] : pull)
          requests.push_back(comm.irecv(owner, 0, pull[owner]));
        for (const auto& [peer, requested] : push)
          requests.push_back(comm.isend(peer, 0, push[peer]));
        boost::mpi::wait_all(requests.begin(), requests.end());

        for (const auto& [owner, requested] : pull)
        {
          for (size_t i = 0; i < requested.size(); ++i)
          {
            const auto local = mesh.getLocalIndex(0, requested[i].first);
            assert(local);
            const Index& global = requested[i].second;
            m_loc2glob.right.insert({ global, *local});
          }
        }

for (int r = 0; r < comm.size(); ++r) {
  if (comm.rank() == r) {
    std::cout << "=== rank " << r << " ===\n";
    for (auto it = m_loc2glob.left.begin(); it != m_loc2glob.left.end(); ++it) {
      std::cout
        << " local=" << it->first
        << " -> global="  << it->second
        << "\n";
    }
    std::cout << std::flush;
  }
  comm.barrier();    // wait so next rank prints cleanly
}


      }

      P1(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_fes(mesh.getShard(), vdim)
      {}

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
      Index getGlobalIndex(Index globalIdx) const
      {
        return m_loc2glob.left.at(globalIdx).get_right();
      }

      /**
       * @brief Returns the local distributed index of the given global
       * distributed index.
       */
      Optional<Index> getLocalIndex(Index globalIdx) const
      {
        auto find = m_loc2glob.right.find(globalIdx);
        if (find == m_loc2glob.right.end())
          return std::nullopt;
        else
          return find->get_left();
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
        boost::optional<ElementType> send;
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
        boost::optional<IndexArray> send;
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
        boost::optional<Index> send;
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
      auto getMapping(const std::pair<size_t, Index>& p, const FunctionBase<FunctionDerived>& v) const
      {
        const auto& [d, globalIdx] = p;
        const auto& mesh = getMesh();
        return Mapping<FunctionDerived>(*mesh.getPolytope(d, globalIdx), v);
      }

      template <class CallableType>
      auto getInverseMapping(const std::pair<size_t, Index>& idx, const CallableType& v) const
      {
        return typename FESType::template InverseMapping<CallableType>(v);
      }

      template <class CallableType>
      auto getInverseMapping(const Geometry::Polytope& polytope, const CallableType& v) const
      {
        return typename FESType::InverseMapping(v);
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;

      size_t m_offset;
      size_t m_owned;
      Bimap m_loc2glob;
  };
}

namespace Rodin::MPI
{
  using P1 = Variational::P1<Real, Geometry::Mesh<Context::MPI>>;
}

#endif
