#ifndef RODIN_MPI_VARIATIONAL_P1_P1_H
#define RODIN_MPI_VARIATIONAL_P1_P1_H

#include <boost/serialization/optional.hpp>

#include "Rodin/Variational/P1/P1.h"

#include "Rodin/MPI/Variational/FiniteElementSpace.h"
#include "Rodin/MPI/Geometry/Mesh.h"

namespace Rodin::Variational
{
  template <class Range>
  class P1<Range, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>, P1<Real, Geometry::Mesh<Context::MPI>>>
  {
    public:
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
            : m_polytope(polytope), m_trans(m_polytope.getTransformation()), m_v(v.copy())
          {}

          Mapping(const Mapping&) = default;

          auto operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
            return getFunction()(p);
          }

          template <class T>
          auto operator()(T& res, const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, m_trans.get(), r);
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
          std::reference_wrapper<const Geometry::PolytopeTransformation> m_trans;
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
            return m_v.get();
          }

        private:
          std::reference_wrapper<const FunctionType> m_v;
      };

      P1(const MeshType& mesh)
        : m_mesh(mesh),
          m_fes(mesh.getShard())
      {
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        size_t local = m_fes.getSize();
        Index inclusiveSum = 0;
        boost::mpi::scan(comm, local, inclusiveSum, std::plus<Index>());
        m_offset = inclusiveSum - local;
      }

      P1(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_fes(mesh.getShard(), vdim)
      {
        const auto& ctx = mesh.getContext();
        const auto& comm = ctx.getCommunicator();
        size_t local = m_fes.getSize();
        Index inclusiveSum = 0;
        boost::mpi::scan(comm, local, inclusiveSum, std::plus<Index>());
        m_offset = inclusiveSum - local;
      }

      P1(const P1& other) = default;

      P1(P1&& other) = default;

      P1& operator=(P1&& other) = default;

      const FESType& getShard() const
      {
        return m_fes;
      }

      /**
       * @brief Returns the global distributed index of the given local
       * distributed index.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        return m_offset + localIdx;
      }

      /**
       * @brief Returns the local distributed index of the given global
       * distributed index.
       */
      std::optional<Index> getLocalIndex(Index globalIdx) const
      {
        const size_t sz = m_fes.getSize();
        if (globalIdx >= m_offset && globalIdx <  m_offset + sz)
          return globalIdx - m_offset;
        else
          return std::nullopt;
      }

      /**
       * @brief Returns the global finite element for the given global index.
       * @param[in] d Dimension of the element
       * @param[in] globalIdx Global distributed index of the polytope in the
       * distributed mesh
       */
      const auto& getFiniteElement(size_t d, Index globalIdx) const
      {
        const auto& mesh = getMesh();
        const auto& shard = mesh.getShard();
        const auto idx = mesh.getLocalIndex(d, globalIdx);
        boost::optional<Geometry::Polytope::Type> local;
        if (idx)
        {
          if (!shard.isGhost(d, *idx))
            local = shard.getGeometry(d, *idx);
        }
        auto res = boost::mpi::all_reduce(
            mesh.getContext().getCommunicator(), local,
            [](auto const& a, auto const& b) { return a ? a : b; });
        assert(res);
        const auto& index = m_fes.getElementIndex();
        return index[*res];
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
        const auto& mesh = getMesh();
        const auto& shard = mesh.getShard();
        const auto& fes = this->getShard();
        const auto idx = mesh.getLocalIndex(d, globalIdx);
        IndexArray local;
        if (idx)
        {
          if (!shard.isGhost(d, *idx))
            local = fes.getDOFs(d, *idx);
        }
        auto dofs = boost::mpi::all_reduce(
            mesh.getContext().getCommunicator(), local,
            [](auto const& a, auto const& b) { return a.size() > 0 ? a : b; });
        assert(dofs.size() > 0);
        for (auto& dof : dofs)
          dof = getGlobalIndex(dof);
        return dofs;
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index localDof) const override
      {
        const auto& [d, globalIdx] = p;
        const auto& mesh = getMesh();
        const auto& shard = mesh.getShard();
        const auto& fes = this->getShard();
        const auto idx = mesh.getLocalIndex(d, globalIdx);
        boost::optional<Index> local;
        if (idx)
        {
          if (!shard.isGhost(d, *idx))
            local = fes.getGlobalIndex({ d, *idx }, localDof);
        }
        auto res = boost::mpi::all_reduce(
            mesh.getContext().getCommunicator(), local,
            [](auto const& a, auto const& b) { return a ? a : b; });
        assert(res);
        return *res;
      }

      template <class FunctionDerived>
      auto getMapping(const std::pair<size_t, Index>& p, const FunctionBase<FunctionDerived>& v) const
      {
        const auto& [d, globalIdx] = p;
        const auto& mesh = getMesh();
        return Mapping<FunctionDerived>(*mesh.getPolytope(d, globalIdx), v);
      }

      template <class FunctionDerived>
      auto getMapping(const Geometry::Polytope& polytope, const FunctionBase<FunctionDerived>& v) const
      {
        return Mapping<FunctionDerived>(polytope, v);
      }

      template <class CallableType>
      auto getInverseMapping(const std::pair<size_t, Index>& idx, const CallableType& v) const
      {
        return typename FESType::InverseMapping(v);
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

      mutable std::vector<FlatMap<Index, IndexArray>> m_dofs;
  };
}

#endif
