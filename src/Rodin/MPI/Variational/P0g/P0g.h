/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_VARIATIONAL_P0G_P0G_H
#define RODIN_MPI_VARIATIONAL_P0G_P0G_H

/**
 * @file
 * @brief Distributed P0g (global constant) finite element space for MPI meshes.
 *
 * This file provides @ref Rodin::Variational::P0g specializations for
 * @ref Rodin::Geometry::Mesh<Rodin::Context::MPI>. The implementation wraps
 * the local shard P0g space and exposes a globally consistent interface.
 *
 * P0g degrees of freedom are globally constant:
 * - Scalar P0g: exactly 1 global DOF (index 0) shared across all ranks.
 * - Vector P0g with @p vdim components: @p vdim global DOFs (indices 0 to
 *   vdim-1) shared across all ranks.
 *
 * Ownership convention for PETSc compatibility: rank 0 owns all global DOFs
 * (@f$[0, \text{vdim})@f$); all other ranks own an empty range
 * (@f$[\text{vdim}, \text{vdim})@f$).  This allows PETSc's ADD_MODE assembly
 * to accumulate contributions from every rank into the single owner's storage.
 *
 * No MPI communication is required during construction.
 */

#include <cassert>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "Rodin/Types.h"
#include "Rodin/Array.h"

#include "Rodin/MPI/Geometry/Mesh.h"
#include "Rodin/MPI/Variational/FiniteElementSpace.h"

#include "Rodin/Variational/P0g/P0g.h"

namespace Rodin::Variational
{
  // --------------------------------------------------------------------------
  // Scalar P0g<Real, Mesh<Context::MPI>>
  // --------------------------------------------------------------------------

  /**
   * @brief Distributed scalar P0g finite element space for MPI meshes.
   *
   * Wraps a local shard scalar P0g space and provides the distributed
   * interface.  All ranks see the same single global DOF (index 0).
   */
  template <>
  class P0g<Real, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>,
        P0g<Real, Geometry::Mesh<Context::MPI>>>
  {
    public:
      using ScalarType  = Real;
      using RangeType   = ScalarType;
      /// @brief Execution context type.
      using ContextType = Context::MPI;
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<RangeType>;

      /// Underlying local shard finite element space type.
      using FESType = P0g<Real, Geometry::Mesh<Context::Local>>;

      /// Parent class.
      using Parent = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

      using Parent::getGlobalIndex;

      // Reuse the local shard Pullback and Pushforward directly.
      template <class Callable>
      using Pullback = typename FESType::template Pullback<Callable>;

      template <class Callable>
      using Pushforward = typename FESType::template Pushforward<Callable>;

      /**
       * @brief Constructs the distributed scalar P0g space on the given mesh.
       *
       * The underlying shard P0g space is created from the local mesh shard.
       * No MPI communication is required.
       *
       * @param[in] mesh Distributed mesh on which the space is defined.
       */
      explicit P0g(const MeshType& mesh)
        : m_mesh(mesh),
          m_fes(mesh.getShard())
      {}

      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_fes(other.m_fes)
      {}

      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh),
          m_fes(std::move(other.m_fes))
      {}

      P0g& operator=(const P0g& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_fes = other.m_fes;
        }
        return *this;
      }

      P0g& operator=(P0g&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = other.m_mesh;
          m_fes = std::move(other.m_fes);
        }
        return *this;
      }

      ~P0g() override = default;

      /**
       * @brief Returns the underlying local shard P0g space.
       */
      const FESType& getShard() const
      {
        return m_fes;
      }

      /**
       * @brief Returns the ownership range @f$[\text{begin}, \text{end})@f$.
       *
       * Rank 0 owns the single global DOF: @p begin = 0, @p end = 1.
       * All other ranks own nothing: @p begin = @p end = 1.
       *
       * @param[out] begin First global DOF owned by this rank.
       * @param[out] end   One past the last global DOF owned by this rank.
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        const int rank = getMesh().getContext().getCommunicator().rank();
        if (rank == 0)
        {
          begin = 0;
          end   = 1;
        }
        else
        {
          begin = 1;
          end   = 1;
        }
      }

      /**
       * @brief Returns the global number of degrees of freedom (always 1).
       */
      size_t getSize() const override
      {
        return 1;
      }

      /**
       * @brief Returns the vector dimension (always 1 for scalar P0g).
       */
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /**
       * @brief Returns the distributed mesh on which the space is defined.
       */
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Returns the finite element attached to local polytope @f$(d, i)@f$.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_fes.getFiniteElement(d, i);
      }

      /**
       * @brief Returns the global DOF array for local polytope @f$(d, i)@f$.
       *
       * For scalar P0g every polytope maps to the single global DOF 0.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_fes.getDOFs(d, i);
      }

      /**
       * @brief Returns the global DOF index for local shard DOF @p localIdx.
       *
       * For scalar P0g there is only 1 global DOF (index 0), so this always
       * returns 0 regardless of @p localIdx.
       *
       * @param[in] localIdx Local shard DOF index (must be 0 for scalar P0g).
       * @return Global distributed DOF index (always 0).
       */
      Index getGlobalIndex(Index) const
      {
        return Index(0);
      }

      /**
       * @brief Returns the global DOF index for local basis function @p localDof
       * on polytope @f$(d, i)@f$.
       *
       * For scalar P0g this always returns 0.
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index localDof) const override
      {
        return m_fes.getGlobalIndex(p, localDof);
      }

      /**
       * @brief Returns a pullback wrapper on local polytope @f$(d, i)@f$.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& p, Callable&& v) const
      {
        return m_fes.getPullback(p, std::forward<Callable>(v));
      }

      /**
       * @brief Returns a pushforward wrapper on local polytope @f$(d, i)@f$.
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& p, Callable&& v) const
      {
        return m_fes.getPushforward(p, std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;
  };

  // --------------------------------------------------------------------------
  // Vector P0g<Math::SpatialVector<Real>, Mesh<Context::MPI>>
  // --------------------------------------------------------------------------

  /**
   * @brief Distributed vector P0g finite element space for MPI meshes.
   *
   * Wraps a local shard vector P0g space and provides the distributed
   * interface.  All ranks see the same @p vdim global DOFs (indices 0 to
   * @p vdim - 1).
   */
  template <>
  class P0g<Math::SpatialVector<Real>, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::MPI>,
        P0g<Math::SpatialVector<Real>, Geometry::Mesh<Context::MPI>>>
  {
    public:
      using ScalarType  = Real;
      using RangeType   = Math::SpatialVector<Real>;
      /// @brief Execution context type.
      using ContextType = Context::MPI;
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<Math::SpatialVector<ScalarType>>;

      /// Underlying local shard finite element space type.
      using FESType = P0g<Math::SpatialVector<Real>, Geometry::Mesh<Context::Local>>;

      /// Parent class.
      using Parent = FiniteElementSpace<MeshType, P0g<Math::SpatialVector<Real>, MeshType>>;

      using Parent::getGlobalIndex;

      template <class Callable>
      using Pullback = typename FESType::template Pullback<Callable>;

      template <class Callable>
      using Pushforward = typename FESType::template Pushforward<Callable>;

      /**
       * @brief Constructs the distributed vector P0g space on the given mesh.
       *
       * @param[in] mesh  Distributed mesh on which the space is defined.
       * @param[in] vdim  Number of vector components.
       */
      P0g(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_fes(mesh.getShard(), vdim)
      {
        assert(vdim > 0);
      }

      /**
       * @brief Integral-constant constructor for compile-time vdim.
       *
       * @tparam VDim  Compile-time vector dimension.
       * @param[in] mesh  Distributed mesh.
       */
      template <size_t VDim>
      P0g(std::integral_constant<size_t, VDim>, const MeshType& mesh)
        : P0g(mesh, VDim)
      {}

      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_fes(other.m_fes)
      {}

      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh),
          m_fes(std::move(other.m_fes))
      {}

      P0g& operator=(const P0g& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_fes = other.m_fes;
        }
        return *this;
      }

      P0g& operator=(P0g&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_fes = std::move(other.m_fes);
        }
        return *this;
      }

      ~P0g() override = default;

      /**
       * @brief Returns the underlying local shard P0g space.
       */
      const FESType& getShard() const
      {
        return m_fes;
      }

      /**
       * @brief Returns the ownership range @f$[\text{begin}, \text{end})@f$.
       *
       * Rank 0 owns all @p vdim global DOFs: @p begin = 0, @p end = vdim.
       * All other ranks own nothing: @p begin = @p end = vdim.
       *
       * @param[out] begin First global DOF owned by this rank.
       * @param[out] end   One past the last global DOF owned by this rank.
       */
      void getOwnershipRange(Index& begin, Index& end) const
      {
        const size_t vdim = getVectorDimension();
        const int rank = getMesh().getContext().getCommunicator().rank();
        if (rank == 0)
        {
          begin = 0;
          end   = static_cast<Index>(vdim);
        }
        else
        {
          begin = static_cast<Index>(vdim);
          end   = static_cast<Index>(vdim);
        }
      }

      /**
       * @brief Returns the global number of degrees of freedom (@p vdim).
       */
      size_t getSize() const override
      {
        return m_fes.getSize();
      }

      /**
       * @brief Returns the vector dimension.
       */
      size_t getVectorDimension() const override
      {
        return m_fes.getVectorDimension();
      }

      /**
       * @brief Returns the distributed mesh on which the space is defined.
       */
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Returns the finite element attached to local polytope @f$(d, i)@f$.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        return m_fes.getFiniteElement(d, i);
      }

      /**
       * @brief Returns the global DOF array for local polytope @f$(d, i)@f$.
       *
       * For vector P0g every polytope maps to DOFs @f$\{0, 1, \ldots, \text{vdim}-1\}@f$.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_fes.getDOFs(d, i);
      }

      /**
       * @brief Returns the global DOF index for local shard DOF @p localIdx.
       *
       * For vector P0g the global DOF index equals the local DOF index
       * @f$(0, 1, \ldots, \text{vdim}-1)@f$, which are the same on all ranks.
       *
       * @param[in] localIdx Local shard DOF index in @f$[0, \text{vdim})@f$.
       * @return Global distributed DOF index.
       */
      Index getGlobalIndex(Index localIdx) const
      {
        assert(localIdx < static_cast<Index>(getVectorDimension()));
        return localIdx;
      }

      /**
       * @brief Returns the global DOF index for local basis function @p localDof.
       *
       * For vector P0g this returns @p localDof (global index equals local index).
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& p, Index localDof) const override
      {
        return m_fes.getGlobalIndex(p, localDof);
      }

      /**
       * @brief Returns a pullback wrapper on local polytope @f$(d, i)@f$.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& p, Callable&& v) const
      {
        return m_fes.getPullback(p, std::forward<Callable>(v));
      }

      /**
       * @brief Returns a pushforward wrapper on local polytope @f$(d, i)@f$.
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& p, Callable&& v) const
      {
        return m_fes.getPushforward(p, std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      FESType m_fes;
  };

} // namespace Rodin::Variational

namespace Rodin::MPI
{
  /**
   * @brief Convenience alias for the default distributed scalar P0g space.
   */
  using P0g = Variational::P0g<Real, Geometry::Mesh<Context::MPI>>;

  /**
   * @brief Convenience alias for the default distributed vector P0g space.
   */
  using VectorP0g = Variational::P0g<Math::SpatialVector<Real>, Geometry::Mesh<Context::MPI>>;
}

#endif
