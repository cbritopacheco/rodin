/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_KOKKOS_VARIATIONAL_P0G_H
#define RODIN_KOKKOS_VARIATIONAL_P0G_H

/**
 * @file
 * @brief Global constant space on a resident Kokkos mesh.
 */

#include <functional>
#include <vector>

#include "Rodin/Variational/P0g/P0g.h"
#include "Rodin/Variational/P0g/Div.h"
#include "Rodin/Variational/P0g/Grad.h"
#include "Rodin/Variational/P0g/Jacobian.h"
#include "Rodin/Variational/P0g/ShapeFunction.h"

#include "Rodin/Kokkos/Geometry/Mesh.h"

#include "DOFs.h"

namespace Rodin::Variational
{
  /**
   * @ingroup P0gSpecializations
   * @brief Real valued global constant space on a resident Kokkos mesh.
   *
   * The space has a single degree of freedom, whose basis function is one
   * on the whole mesh, as the host P0g space. It is typically the Lagrange
   * multiplier fixing the mean of a pressure.
   *
   * @tparam Backend Kokkos backend of the mesh context.
   */
  template <class Backend>
  class P0g<Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>,
        P0g<Real, Geometry::Mesh<Context::Kokkos<Context::Local, Backend>>>>
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
      using ElementType = P0gElement<RangeType>;

      /// @brief Parent class.
      using Parent = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

      /// @brief Device degrees of freedom.
      using DOFsType = Rodin::Kokkos::Variational::ConstantDOFs;

      /**
       * @brief Constructs the space on a resident mesh.
       * @param[in] mesh Mesh supporting the space, which must outlive it.
       */
      explicit P0g(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Space to copy.
       */
      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Space to move from.
       */
      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh))
      {}

      /**
       * @brief Copy assignment operator.
       * @param[in] other Space to copy.
       * @return Reference to this space.
       */
      P0g& operator=(const P0g& other)
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
      P0g& operator=(P0g&& other)
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
       * @return Constant element of the polytope's reference geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        static const auto s_elements = [] {
          std::vector<ElementType> elements;
          elements.reserve(Geometry::Polytope::Types.size());
          for (const auto g : Geometry::Polytope::Types)
            elements.emplace_back(g);
          return elements;
        }();
        return s_elements[static_cast<size_t>(getMesh().getGeometry(d, i))];
      }

      /// @brief Gets the number of degrees of freedom. @return One.
      size_t getSize() const override
      {
        return 1;
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
       * @brief Gets the global degrees of freedom of a polytope.
       * @param[in] d Dimension of the polytope.
       * @param[in] i Index of the polytope.
       * @return The single degree of freedom.
       */
      const IndexArray& getDOFs(
        [[maybe_unused]] size_t d, [[maybe_unused]] Index i) const override
      {
        static const IndexArray s_dofs{{0}};
        return s_dofs;
      }

      /**
       * @brief Converts a local degree of freedom to its global index.
       * @param[in] idx Pair of (dimension, index) of the polytope.
       * @param[in] local Local degree of freedom on the polytope.
       * @return Zero.
       */
      Index getGlobalIndex([[maybe_unused]] const std::pair<size_t, Index>& idx,
        [[maybe_unused]] Index local) const override
      {
        return 0;
      }

      /**
       * @brief Builds the host space with the same numbering.
       * @param[in] mesh Host copy of the mesh, see Geometry::Mesh::toHost().
       * @return P0g space on @p mesh.
       */
      P0g<RangeType, Geometry::Mesh<Context::Local>> toHost(
        const Geometry::Mesh<Context::Local>& mesh) const
      {
        return P0g<RangeType, Geometry::Mesh<Context::Local>>(mesh);
      }

      /// @brief Gets the device degrees of freedom. @return Single constant basis function.
      DOFsType getDeviceDOFs() const
      {
        return {};
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
  };
}

namespace Rodin::Kokkos::Variational
{
  /// @brief Global constant space on a resident mesh.
  template <class Backend>
  using DeviceP0g = ::Rodin::Variational::P0g<Real,
    ::Rodin::Geometry::Mesh<::Rodin::Context::Kokkos<::Rodin::Context::Local, Backend>>>;
}

#endif
