/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MPI_VARIATIONAL_H1_H
#define RODIN_KOKKOS_MPI_VARIATIONAL_H1_H

/**
 * @file
 * @brief H1 Lagrange spaces on a distributed resident Kokkos mesh.
 */

#include <functional>
#include <type_traits>

#include "Rodin/MPI/Variational/H1.h"

#include "Rodin/Kokkos/Variational/H1.h"

#include "Rodin/Kokkos/MPI/Geometry/Mesh.h"

namespace Rodin::Variational
{
  /**
   * @ingroup H1Specializations
   * @brief H1 Lagrange space on a distributed resident Kokkos mesh.
   *
   * The space pairs the H1 space of the resident shard, which numbers the
   * local degrees of freedom of the kernels, with the distributed host H1
   * space, which numbers them across processes. Both spaces number the
   * shard alike, since the resident shard keeps the polytopes of the host
   * shard in the same order. The host members answer as the distributed
   * host space does: getSize() is the global size, and getDOFs() and
   * getGlobalIndex() return global indices. The device members,
   * getDeviceDOFs() and getLocalSize(), refer to the local numbering.
   *
   * @tparam K Polynomial degree.
   * @tparam Range Real for a scalar space, Math::SpatialVector<Real> for a
   * vector space.
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <size_t K, class Range, class Backend>
  class H1<K, Range, Geometry::Mesh<Context::Kokkos<Context::MPI, Backend>>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::MPI, Backend>>,
        H1<K, Range, Geometry::Mesh<Context::Kokkos<Context::MPI, Backend>>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = Range;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::MPI, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Distributed host space.
      using HostType = H1<K, RangeType, Geometry::Mesh<Context::MPI>>;

      /// @brief Resident space of the shard.
      using LocalType =
        H1<K, RangeType, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>;

      /// @brief Type of finite element.
      using ElementType = typename LocalType::ElementType;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType = typename LocalType::DOFsType;

      using Parent::getGlobalIndex;

      /**
       * @brief Constructs the scalar space; collective over the processes.
       * @param[in] degree Polynomial degree.
       * @param[in] mesh Mesh without pyramids supporting the space, which must
       * outlive it.
       */
      H1(std::integral_constant<size_t, K> degree, const MeshType& mesh)
      requires std::is_same_v<RangeType, Real>
        : m_mesh(mesh),
          m_host(degree, mesh.getHost()),
          m_local(degree, mesh.getLocal())
      {}

      /**
       * @brief Constructs the vector space; collective over the processes.
       * @param[in] degree Polynomial degree.
       * @param[in] mesh Mesh without pyramids supporting the space, which must
       * outlive it.
       * @param[in] vdim Vector dimension.
       */
      H1(std::integral_constant<size_t, K> degree, const MeshType& mesh, size_t vdim)
      requires(!std::is_same_v<RangeType, Real>)
        : m_mesh(mesh),
          m_host(degree, mesh.getHost(), vdim),
          m_local(degree, mesh.getLocal(), vdim)
      {}

      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_host(other.m_host),
          m_local(other.m_local)
      {}

      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_host(std::move(other.m_host)),
          m_local(std::move(other.m_local))
      {}

      /// @param[in] d Dimension of the polytope. @param[in] i Local index of the polytope. @return Element of its reference geometry.
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_local.getFiniteElement(d, i);
      }

      /// @brief Gets the global number of degrees of freedom. @return Global size.
      size_t getSize() const override
      {
        return m_host.getSize();
      }

      /// @brief Gets the number of local degrees of freedom, ghosts included. @return Local size.
      size_t getLocalSize() const
      {
        return m_local.getSize();
      }

      size_t getVectorDimension() const override
      {
        return m_local.getVectorDimension();
      }

      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /// @returns Distributed host space.
      const HostType& getHost() const
      {
        return m_host;
      }

      /// @returns Resident space of the shard.
      const LocalType& getLocal() const
      {
        return m_local;
      }

      /// @returns Host space of the shard, numbering the local degrees of freedom.
      const auto& getShard() const
      {
        return m_host.getShard();
      }

      /// @param[out] begin First owned global index. @param[out] end One past the last.
      void getOwnershipRange(Index& begin, Index& end) const
      {
        m_host.getOwnershipRange(begin, end);
      }

      /// @param[in] local Local degree of freedom. @return Global index.
      Index getGlobalIndex(Index local) const
      {
        return m_host.getGlobalIndex(local);
      }

      /// @param[in] global Global degree of freedom. @return Local index, if local.
      Optional<Index> getLocalIndex(Index global) const
      {
        return m_host.getLocalIndex(global);
      }

      /// @param[in] d Dimension of the polytope. @param[in] i Local index of the polytope. @return Global indices of its degrees of freedom.
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_host.getDOFs(d, i);
      }

      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        return m_host.getGlobalIndex(idx, local);
      }

      /**
       * @brief Gets the pullback of a callable on a polytope, as on the host.
       * @param[in] idx Pair of dimension and local polytope index.
       * @param[in] v Callable of the physical points.
       * @return Callable of the reference coordinates.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return m_local.getPullback(idx, std::forward<Callable>(v));
      }

      /// @brief Gets the device degrees of freedom. @return Device numbering and basis.
      const DOFsType& getDeviceDOFs() const
      {
        return m_local.getDeviceDOFs();
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      HostType m_host;
      LocalType m_local;
  };
}

#endif
