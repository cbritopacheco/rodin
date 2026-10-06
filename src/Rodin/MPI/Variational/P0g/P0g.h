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

#include <map>
#include <cassert>
#include <cstddef>
#include <functional>
#include <type_traits>
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
  // Scalar P0g<Real or Complex, Mesh<Context::MPI>>
  // --------------------------------------------------------------------------

  /**
   * @brief Distributed scalar P0g finite element space for MPI meshes.
   *
   * Wraps a local shard scalar P0g space and provides the distributed
   * interface.  All ranks see the same single global DOF (index 0).
   */
  template <class Scalar>
    requires(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>)
  class P0g<Scalar, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::MPI>,
        P0g<Scalar, Geometry::Mesh<Context::MPI>>>
  {
    public:
      /// @brief Scalar coefficient type.
      using ScalarType = Scalar;
      /// @brief Value type represented by the finite element space.
      using RangeType   = ScalarType;
      /// @brief Execution context type.
      using ContextType = Context::MPI;
      /// @brief Distributed mesh type.
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<RangeType>;

      /// Underlying local shard finite element space type.
      using FESType = P0g<Scalar, Geometry::Mesh<Context::Local>>;

      /// Parent class.
      using Parent = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

      using Parent::getGlobalIndex;

      /// @brief Pullback type reused from the local shard finite element space.
      template <class Callable>
      using Pullback = typename FESType::template Pullback<Callable>;

      /// @brief Pushforward type reused from the local shard finite element space.
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

      /// @brief Copy-constructs the distributed scalar P0g space.
      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_fes(other.m_fes)
      {}

      /// @brief Move-constructs the distributed scalar P0g space.
      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh),
          m_fes(std::move(other.m_fes))
      {}

      /// @brief Copy-assigns the distributed scalar P0g space.
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

      /// @brief Move-assigns the distributed scalar P0g space.
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
       * returns 0 regardless of the local shard DOF index.
       *
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
        // Reuse the shared pullback implementation, but retain the MPI mesh
        // identity. A shard-attached point loses SubMesh ancestry and cannot
        // be included in a parent GridFunction's mesh.
        const auto& [d, i] = p;
        return Pullback<Callable>(*getMesh().getPolytope(d, i), std::forward<Callable>(v));
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
  // Vector P0g<Math::SpatialVector<Real or Complex>, Mesh<Context::MPI>>
  // --------------------------------------------------------------------------

  /**
   * @brief Distributed vector P0g finite element space for MPI meshes.
   *
   * Wraps a local shard vector P0g space and provides the distributed
   * interface.  All ranks see the same @p vdim global DOFs (indices 0 to
   * @p vdim - 1).
   */
  template <class Scalar>
    requires(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>)
  class P0g<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::MPI>,
        P0g<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::MPI>>>
  {
    public:
      /// @brief Scalar coefficient type.
      using ScalarType = Scalar;
      /// @brief Vector value type represented by the finite element space.
      using RangeType = Math::SpatialVector<Scalar>;
      /// @brief Execution context type.
      using ContextType = Context::MPI;
      /// @brief Distributed mesh type.
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<Math::SpatialVector<ScalarType>>;

      /// Underlying local shard finite element space type.
      using FESType = P0g<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::Local>>;

      /// Parent class.
      using Parent =
        FiniteElementSpace<MeshType, P0g<Math::SpatialVector<Scalar>, MeshType>>;

      using Parent::getGlobalIndex;

      /// @brief Pullback type reused from the local shard finite element space.
      template <class Callable>
      using Pullback = typename FESType::template Pullback<Callable>;

      /// @brief Pushforward type reused from the local shard finite element space.
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

      /// @brief Copy-constructs the distributed vector P0g space.
      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_fes(other.m_fes)
      {}

      /// @brief Move-constructs the distributed vector P0g space.
      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh),
          m_fes(std::move(other.m_fes))
      {}

      /// @brief Copy-assigns the distributed vector P0g space.
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

      /// @brief Move-assigns the distributed vector P0g space.
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
        // The mathematical pullback is shared; point provenance belongs to
        // the MPI mesh, not its rank-local shard (as for MPI P0 and H1).
        const auto& [d, i] = p;
        return Pullback<Callable>(*getMesh().getPolytope(d, i), std::forward<Callable>(v));
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

namespace Rodin::Variational
{
  /// @brief Matrix-range finite element or expression specialization.
  template <class Scalar>
  class P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::MPI>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::MPI>,
        P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::MPI>>>
  {
    public:
      /// @brief Local matrix space on the mesh shard.
      using FESType = P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>>;
      /// @brief Mesh and execution context of this specialization.
      using MeshType = Geometry::Mesh<Context::MPI>;
      /// @brief Scalar space supplying this family's DOF numbering.
      using ScalarSpace = P0g<Scalar, Geometry::Mesh<Context::MPI>>;
      /// @brief Matrix reference element for this family.
      using ElementType = P0gElement<Math::SpatialMatrix<Scalar>>;
      /// @brief Existing finite-element-space interface.
      using Parent = FiniteElementSpace<MeshType,
        P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::MPI>>>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename ScalarSpace::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Local or distributed execution context.
      using ContextType = typename ScalarSpace::ContextType;
      /// @brief Expands scalar DOF maps into interleaved row-major matrix components.
      P0g(const MeshType& mesh, size_t rows, size_t cols)
        : m_scalar(ScalarSpace(mesh)),
          m_rows(rows),
          m_cols(cols),
          m_shard(mesh.getShard(), rows, cols)
      {
        if (rows == 0 || cols == 0 || rows > RODIN_MAXIMAL_SPACE_DIMENSION ||
          cols > RODIN_MAXIMAL_SPACE_DIMENSION)
          Alert::Exception() << "SpatialMatrix ranges require 1 to 3 rows and columns."
                             << Alert::Raise;
        // Distributed scalar spaces may compute their size collectively.
        // Resolve it at construction, never from a local evaluation loop.
        m_size = m_scalar.getSize() * rows * cols;
        m_dofs.resize(mesh.getDimension() + 1);
        for (size_t d = 0; d <= mesh.getDimension(); ++d)
        {
          const size_t count = mesh.getConnectivity().getCount(d);
          m_dofs[d].reserve(count);
          for (size_t i = 0; i < count; ++i)
          {
            const auto& scalarDOFs = m_scalar.getDOFs(d, i);
            auto& dofs = m_dofs[d].emplace_back(scalarDOFs.size() * rows * cols);
            for (size_t a = 0; a < static_cast<size_t>(scalarDOFs.size()); ++a)
              for (size_t c = 0; c < rows * cols; ++c)
                dofs[a * rows * cols + c] = scalarDOFs[a] * rows * cols + c;
            const auto& scalarFE = m_scalar.getFiniteElement(d, i);
            m_elements.try_emplace(scalarFE.getGeometry(), scalarFE, rows, cols);
          }
        }
        for (size_t i = 0; i < m_shard.getSize(); ++i)
          m_globalToLocal.emplace(getGlobalIndex(i), i);
      }

      /// @brief Copies the space and its DOF maps.
      P0g(const P0g&) = default;
      /// @brief Moves the space and its DOF maps.
      P0g(P0g&&) = default;
      /// @brief Copies the space and its DOF maps.
      P0g& operator=(const P0g&) = default;
      /// @brief Moves the space and its DOF maps.
      P0g& operator=(P0g&&) = default;

      size_t getSize() const override
      {
        return m_size;
      }
      size_t getVectorDimension() const override
      {
        return m_rows * m_cols;
      }
      /// @brief Returns the number of matrix rows.
      size_t getRows() const
      {
        return m_rows;
      }
      /// @brief Returns the number of matrix columns.
      size_t getColumns() const
      {
        return m_cols;
      }
      const MeshType& getMesh() const override
      {
        return m_scalar.getMesh();
      }
      /// @brief Returns the scalar space supplying topology and component-independent maps.
      const ScalarSpace& getScalarSpace() const
      {
        return m_scalar;
      }

      /// @brief Returns the matrix reference element of a mesh entity.
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

      /// @brief Pulls a physical callable back to a reference element.
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& p, Callable&& value) const
      {
        return m_scalar.getPullback(p, std::forward<Callable>(value));
      }

      /// @brief Pushes a reference callable forward to the physical mesh.
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& p, Callable&& value) const
      {
        return m_scalar.getPushforward(p, std::forward<Callable>(value));
      }

      /// @brief Returns the matrix space on the local mesh shard.
      const FESType& getShard() const
      {
        return m_shard;
      }

      /// @brief Returns the half-open range of owned global component DOFs.
      void getOwnershipRange(Index& begin, Index& end) const
      {
        this->getScalarSpace().getOwnershipRange(begin, end);
        begin *= this->getVectorDimension();
        end *= this->getVectorDimension();
      }

      /// @brief Maps a local component DOF to its global coefficient index.
      Index getGlobalIndex(Index local) const
      {
        assert(static_cast<size_t>(local) < m_shard.getSize());
        const Index components = this->getVectorDimension();
        return this->getScalarSpace().getGlobalIndex(local / components) * components +
          local % components;
      }

      /// @brief Returns the shard-local index of a global DOF, if present.
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
      size_t m_size = 0;
      std::vector<std::vector<IndexArray>> m_dofs;
      std::map<Geometry::Polytope::Type, ElementType> m_elements;
  };
}

#endif
