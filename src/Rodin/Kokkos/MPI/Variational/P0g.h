/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_MPI_VARIATIONAL_P0G_H
#define RODIN_KOKKOS_MPI_VARIATIONAL_P0G_H

/**
 * @file
 * @brief Global constant space on a distributed resident Kokkos mesh.
 */

#include <functional>
#include <utility>

#include "Rodin/MPI/Variational/P0g.h"

#include "Rodin/Kokkos/Variational/P0g.h"

#include "Rodin/Kokkos/MPI/Geometry/Mesh.h"

namespace Rodin::Variational
{
  /**
   * @ingroup P0gSpecializations
   * @brief Real valued global constant space on a distributed resident
   * Kokkos mesh.
   *
   * The space has a single degree of freedom, held by every process and
   * owned by the first one, as the distributed host P0g space. Its local and
   * global indices are both zero. It is typically the Lagrange multiplier
   * fixing the mean of a pressure, each process integrating its owned
   * polytopes into the shared row.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class P0g<Real, Geometry::Mesh<Context::Kokkos<Context::MPI, Backend>>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::MPI, Backend>>,
        P0g<Real, Geometry::Mesh<Context::Kokkos<Context::MPI, Backend>>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = ScalarType;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::MPI, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Distributed host space.
      using HostType = P0g<RangeType, Geometry::Mesh<Context::MPI>>;

      /// @brief Resident space of the shard.
      using LocalType =
        P0g<RangeType, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>;

      /// @brief Type of finite element.
      using ElementType = typename LocalType::ElementType;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType = typename LocalType::DOFsType;

      using Parent::getGlobalIndex;

      /**
       * @brief Constructs the space; no communication is needed.
       * @param[in] mesh Mesh supporting the space, which must outlive it.
       */
      explicit P0g(const MeshType& mesh)
        : m_mesh(mesh),
          m_host(mesh.getHost()),
          m_local(mesh.getLocal())
      {}

      /**
       * @brief Copies the space, sharing its mesh.
       * @param other Space to copy.
       */
      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_host(other.m_host),
          m_local(other.m_local)
      {}

      /**
       * @brief Moves the space.
       * @param other Space to move.
       */
      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_host(std::move(other.m_host)),
          m_local(std::move(other.m_local))
      {}

      /// @param[in] d Dimension of the polytope. @param[in] i Local index of the polytope. @return Constant element.
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_local.getFiniteElement(d, i);
      }

      /// @brief Gets the global number of degrees of freedom. @return One.
      size_t getSize() const override
      {
        return 1;
      }

      /// @brief Gets the number of local degrees of freedom. @return One.
      size_t getLocalSize() const
      {
        return 1;
      }

      /// @brief Gets the vector dimension of the range. @return One.
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /// @brief Gets the mesh supporting the space. @return Distributed resident mesh.
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

      /// @returns Host space of the shard.
      const auto& getShard() const
      {
        return m_host.getShard();
      }

      /// @param[out] begin First owned global index. @param[out] end One past the last.
      void getOwnershipRange(Index& begin, Index& end) const
      {
        m_host.getOwnershipRange(begin, end);
      }

      /// @param[in] local Local degree of freedom. @return Zero.
      Index getGlobalIndex([[maybe_unused]] Index local) const
      {
        return 0;
      }

      /// @param[in] global Global degree of freedom. @return Zero, the local index of the global one.
      Optional<Index> getLocalIndex(Index global) const
      {
        if (global == 0)
          return Index(0);
        return {};
      }

      /// @param[in] d Dimension of the polytope. @param[in] i Local index of the polytope. @return The single global degree of freedom.
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_host.getDOFs(d, i);
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of (dimension, local index) of the polytope.
       * @param[in] local Local degree of freedom on the polytope.
       * @return Zero.
       */
      Index getGlobalIndex([[maybe_unused]] const std::pair<size_t, Index>& idx,
        [[maybe_unused]] Index local) const override
      {
        return 0;
      }

      /// @brief Gets the device degrees of freedom. @return Single constant basis function.
      DOFsType getDeviceDOFs() const
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
