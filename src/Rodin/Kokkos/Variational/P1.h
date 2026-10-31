/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_P1_H
#define RODIN_KOKKOS_VARIATIONAL_P1_H

/**
 * @file
 * @brief P1 Lagrange spaces on a resident Kokkos mesh.
 */

#include <algorithm>
#include <array>
#include <functional>
#include <vector>

#include "Rodin/Variational/P1/P1.h"

#include "Rodin/Kokkos/Geometry/Mesh.h"

#include "DOFs.h"

namespace Rodin::Variational
{
  /**
   * @ingroup P1Specializations
   * @brief Real valued P1 Lagrange space on a resident Kokkos mesh.
   *
   * The space numbers one degree of freedom per mesh vertex, as the host P1
   * space does, so a local mesh and the resident mesh built from it share
   * their global numbering. Local degrees of freedom follow the vertex order
   * of each polytope.
   *
   * The host members answer the form language queries made while building
   * and copying forms. Device evaluation of the basis goes through
   * Rodin::Kokkos::Geometry::PointView::evaluateBasis(), which is the P1
   * Lagrange basis of each reference geometry.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class P1<Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        P1<Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = ScalarType;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::Local, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Type of finite element.
      using ElementType = P1Element<RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType = Rodin::Kokkos::Variational::VertexDOFs;

      /**
       * @brief Constructs the P1 space on a resident mesh.
       * @param[in] mesh Mesh supporting the space, which must outlive it.
       */
      P1(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      P1& operator=(const P1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Space to move from.
       * @return Reference to this space.
       */
      P1& operator=(P1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
        }
        return *this;
      }

      /**
       * @brief Gets the finite element of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return P1 element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        static const auto s_elements = [] {
          std::array<ElementType, Geometry::Polytope::Types.size()> elements;
          for (size_t g = 0; g < elements.size(); ++g)
            elements[g] = ElementType(Geometry::Polytope::Types[g]);
          return elements;
        }();
        return s_elements[static_cast<size_t>(getMesh().getGeometry(d, i))];
      }

      /// @brief Gets the number of degrees of freedom. @return Vertex count.
      size_t getSize() const override
      {
        return getMesh().getVertexCount();
      }

      /// @brief Gets the vector dimension. @return One, the space being scalar.
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /// @brief Gets the underlying mesh. @return Resident mesh.
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return Scalar P1 space on @p mesh.
       */
      P1<RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return P1<RangeType, Geometry::Mesh<Context::Local>>(mesh);
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Ordered vertex indices of the polytope.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        static thread_local IndexArray s_dofs;
        const auto& key = getMesh().getConnectivity().getPolytope(d, i);
        s_dofs.resize(key.size());
        std::copy(key.begin(), key.end(), s_dofs.begin());
        return s_dofs;
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] local Local vertex index within the polytope.
       * @return Global vertex index.
       */
      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        const auto& key = getMesh().getConnectivity().getPolytope(d, i);
        assert(local < static_cast<size_t>(key.size()));
        return key(local);
      }

      /// @brief Gets the device degrees of freedom. @return One per vertex.
      DOFsType getDeviceDOFs() const
      {
        return {static_cast<Index>(getMesh().getVertexCount())};
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
  };

  /**
   * @ingroup P1Specializations
   * @brief Vector valued P1 Lagrange space on a resident Kokkos mesh.
   *
   * The space numbers its degrees of freedom as the host vector P1 space
   * does: component @f$ r @f$ of vertex @f$ q @f$ has the global index
   * @f$ q + r n_v @f$, where @f$ n_v @f$ is the vertex count, and the local
   * index @f$ a\,d + r @f$ when @f$ q @f$ is the local vertex @f$ a @f$ of a
   * polytope and @f$ d @f$ the vector dimension.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class P1<Math::SpatialVector<Real>,
    Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>
    final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        P1<Math::SpatialVector<Real>,
          Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;

      /// @brief Range type of value.
      using RangeType = Math::SpatialVector<ScalarType>;

      /// @brief Context of the mesh.
      using ContextType = Context::Kokkos<Context::Local, Backend>;

      /// @brief Type of mesh on which the space is built.
      using MeshType = Geometry::Mesh<ContextType>;

      /// @brief Type of finite element.
      using ElementType = P1Element<RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType = Rodin::Kokkos::Variational::VertexDOFs;

      /**
       * @brief Constructs the vector P1 space on a resident mesh.
       * @param[in] mesh Mesh supporting the space, which must outlive it.
       * @param[in] vdim Vector dimension, at most three.
       */
      P1(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_vdim(vdim)
      {
        assert(vdim <= RODIN_MAXIMAL_SPACE_DIMENSION);
      }

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(other.m_vdim)
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      P1& operator=(const P1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_vdim = other.m_vdim;
        }
        return *this;
      }

      /**
       * @brief Move assignment operator.
       * @param[in] other Space to move from.
       * @return Reference to this space.
       */
      P1& operator=(P1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = other.m_vdim;
        }
        return *this;
      }

      /**
       * @brief Gets the finite element of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Vector P1 element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        constexpr size_t Dimensions = RODIN_MAXIMAL_SPACE_DIMENSION + 1;
        static const auto s_elements = [] {
          std::vector<ElementType> elements;
          elements.reserve(Geometry::Polytope::Types.size() * Dimensions);
          for (const auto g : Geometry::Polytope::Types)
          {
            for (size_t vdim = 0; vdim < Dimensions; ++vdim)
              elements.emplace_back(g, vdim);
          }
          return elements;
        }();
        const auto g = static_cast<size_t>(getMesh().getGeometry(d, i));
        return s_elements[g * Dimensions + m_vdim];
      }

      /// @brief Gets the number of degrees of freedom. @return Vertex count times vector dimension.
      size_t getSize() const override
      {
        return getMesh().getVertexCount() * m_vdim;
      }

      /// @brief Gets the vector dimension. @return Number of components.
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      /// @brief Gets the underlying mesh. @return Resident mesh.
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return Vector P1 space on @p mesh.
       */
      P1<RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return P1<RangeType, Geometry::Mesh<Context::Local>>(mesh, m_vdim);
      }

      /**
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return Global indices in local order.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        static thread_local IndexArray s_dofs;
        const auto& key = getMesh().getConnectivity().getPolytope(d, i);
        const size_t count = key.size();
        const size_t vn = getMesh().getVertexCount();
        s_dofs.resize(count * m_vdim);
        for (size_t local = 0; local < count * m_vdim; ++local)
          s_dofs.coeffRef(local) = key(local / m_vdim) + (local % m_vdim) * vn;
        return s_dofs;
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of dimension and polytope index.
       * @param[in] local Local index within the polytope.
       * @return Global index.
       */
      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        const auto& key = getMesh().getConnectivity().getPolytope(d, i);
        assert(local / m_vdim < static_cast<size_t>(key.size()));
        return key(local / m_vdim) + (local % m_vdim) * getMesh().getVertexCount();
      }

      /// @brief Gets the device degrees of freedom. @return One per vertex.
      DOFsType getDeviceDOFs() const
      {
        return {static_cast<Index>(getMesh().getVertexCount())};
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_vdim;
  };
}

#endif
