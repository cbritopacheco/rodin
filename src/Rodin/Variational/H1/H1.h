/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1_H
#define RODIN_VARIATIONAL_H1_H1_H

#include <cstdint>
#include <functional>
#include <type_traits>

#include <boost/multi_array.hpp>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"

#include "Rodin/Utility/ForConstexpr.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "H1Element.h"
#include "Rodin/Variational/H1/Fekete.h"
#include "Rodin/Variational/H1/GLL.h"

#include "Rodin/Utility/DependentValue.h"

namespace Rodin::FormLanguage
{
  template <size_t K, class Scalar, class Mesh>
  struct Traits<Variational::H1<K, Scalar, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = ScalarType;
    using ElementType = Variational::H1Element<K, RangeType>;
  };

  template <size_t K, class Scalar, class Mesh>
  struct Traits<Variational::H1<K, Math::Vector<Scalar>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = Math::Vector<ScalarType>;
    using ElementType = Variational::H1Element<K, RangeType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup H1Specializations H1 Template Specializations
   * @brief Template specializations of the H1 class.
   * @see H1
   */

  template <size_t K, class Range, class Mesh = Geometry::Mesh<Context::Local>>
  class H1;

  /**
   * @ingroup H1Specializations
   * @brief Degree K H1-conforming Lagrange finite element space
   *
   * Represents the finite element space composed of scalar valued continuous,
   * piecewise polynomial functions of degree K:
   * @f[
   *  \mathbb{P}_K (\mathcal{T}_h) = \{ v \in C^0(\mathcal{T}_h) : v|_{\tau} \in \mathbb{P}_K(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is scalar valued, i.e. evaluations of the function are of
   * Rodin::Real type.
   *
   * ## Building the H1 Space
   *
   * The H1 space is built by distributing degrees of freedom (DOFs) across the
   * mesh entities based on polynomial degree K:
   * - **K ≥ 1**: Each vertex carries exactly 1 DOF (nodal value)
   * - **K ≥ 2**: Each edge carries K-1 interior DOFs (in addition to endpoint vertices)
   * - **K ≥ 3**: Each face (2D) carries interior DOFs (triangle: (K-1)(K-2)/2, quad: (K-1)²)
   * - **K ≥ 4**: Each cell (3D) carries interior DOFs
   *
   * The construction requires the mesh connectivity chain from the highest
   * dimension D down to dimension 0:
   * @f[
   *  D \to D-1 \to \ldots \to 1 \to 0
   * @f]
   *
   * This connectivity is necessary because DOF numbering proceeds top-down:
   * starting from cells (dimension D), we recursively number DOFs on faces,
   * edges, and vertices while respecting H1 conformity (shared entities have
   * shared DOFs).
   *
   * ### Nodal Ordering
   *
   * The nodal basis functions are ordered using:
   * - **Segments**: Gauss-Lobatto-Legendre (GLL) points
   * - **Triangles**: Fekete points (optimal for interpolation)
   * - **Tetrahedra**: Fekete points
   * - **Quadrilaterals**: Tensor-product GLL grid
   * - **Wedges**: Tensor-product (triangle Fekete × segment GLL)
   *
   * @tparam K Polynomial degree (1, 2, 3, ...)
   * @tparam Scalar Scalar type (Real, Complex)
   *
   * @see getClosure() for DOF numbering details
   * @see Cochain for local-to-global DOF mapping on reference elements
   */
  template <size_t K, class Scalar>
  class H1<K, Scalar, Geometry::Mesh<Context::Local>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, H1<K, Scalar, Geometry::Mesh<Context::Local>>>
  {
    public:
      static_assert(K > 0, "Polynomial degree K must be greater than 0.");

      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the H1 space
      using ContextType = Context::Local;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = H1Element<K, RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      /**
       * @brief Local cochain map for H1 DOF injection from boundary to element.
       *
       * The Cochain class encodes how degrees of freedom from boundary entities
       * (vertices, edges, faces) are injected into the DOF array of a polytope.
       * This is essential for building conforming H1 spaces where shared entities
       * have shared DOFs.
       *
       * For example, when numbering DOFs on a triangle, the Cochain::map method
       * injects vertex and edge DOFs into the triangle's local DOF array,
       * respecting the reference element's nodal ordering.
       *
       * @tparam G Polytope type (Segment, Triangle, Quadrilateral, Tetrahedron, Wedge)
       *
       * @see getClosure() which uses Cochain::map during DOF numbering
       */
      template <Geometry::Polytope::Type G>
      class Cochain;

      /**
       * @brief Pullback for the scalar/complex H1 space.
       */
      template <class Callable>
      class Pullback : public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          Pullback(const Pullback&) = default;

          decltype(auto) operator()(const Math::SpatialVector<Real>& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      /**
       * @brief Inverse Pullback for the scalar/complex H1 space.
       */
      template <class Callable>
      class Pushforward
        : public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          using CallableType = Callable;

          /**
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          template <class Function>
          Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          Pushforward(const Pushforward&) = default;

          constexpr
          decltype(auto) operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      /**
       * @brief Constructs an H1 finite element space on the given mesh.
       * @param[in] mesh Mesh on which to build the finite element space
       *
       * Creates the H1 space of degree K. The total number of DOFs depends on
       * the polynomial degree and mesh topology.
       *
       * ### Mesh Connectivity Requirements
       *
       * The mesh must have the complete connectivity chain from dimension D
       * (mesh dimension) down to dimension 0:
       * @f[
       *   D \to D-1 \to D-2 \to \ldots \to 1 \to 0
       * @f]
       *
       * This is required because DOF numbering proceeds recursively from cells
       * to faces to edges to vertices, ensuring conformity (shared entities
       * have shared DOFs).
       *
       * #### Example for 2D Triangle Mesh
       * Required connectivity:
       * - 2 → 1 (triangles to edges)
       * - 1 → 0 (edges to vertices)
       *
       * #### Example for 3D Tetrahedral Mesh
       * Required connectivity:
       * - 3 → 2 (tetrahedra to faces)
       * - 2 → 1 (faces to edges)
       * - 1 → 0 (edges to vertices)
       *
       * The constructor will automatically verify and require these connectivity
       * relationships using RODIN_GEOMETRY_REQUIRE_INCIDENCE.
       *
       * @throws Alert::Exception if required connectivity is not available
       *
       * @see getClosure() for the DOF numbering algorithm
       */
      H1(std::integral_constant<size_t, K>, const MeshType& mesh);

      /**
       * @brief Copy constructor.
       * @param[in] other H1 space to copy
       */
      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_closure(other.m_dofs),
          m_size(other.m_size)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other H1 space to move from
       */
      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_closure(std::move(other.m_dofs)),
          m_size(std::move(other.m_size))
      {}

      virtual ~H1() = default;

      /**
       * @brief Move assignment operator.
       * @param[in] other H1 space to move from
       * @return Reference to this H1 space
       */
      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_closure = std::move(other.m_dofs);
          m_size = std::move(other.m_size);
        }
        return *this;
      }

      /**
       * @brief Copy assignment operator.
       * @param[in] other H1 space to copy
       * @return Reference to this H1 space
       */
      H1& operator=(const H1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_closure = other.m_dofs;
          m_size = other.m_size;
        }
        return *this;
      }

      /**
       * @brief Recursively builds the global DOF closure for a polytope.
       *
       * This method assigns global DOF indices to a polytope (dimension d, index idx)
       * by recursively processing its boundary entities (faces, edges, vertices).
       * The process ensures H1 conformity: shared entities receive shared DOFs.
       *
       * ### Algorithm Overview
       *
       * 1. **Check if already visited**: If this polytope has been numbered, return.
       * 2. **Process boundary entities**: Recursively call getClosure on incident
       *    entities of dimension d-1 (faces for cells, edges for faces, vertices for edges).
       * 3. **Inject boundary DOFs**: Use Cochain::map to inject boundary DOFs into
       *    this polytope's local DOF array.
       * 4. **Number interior DOFs**: Assign new global DOF indices to any DOFs not
       *    on the boundary (interior edge DOFs, face DOFs, cell DOFs).
       *
       * ### Connectivity Requirements
       *
       * This method requires mesh connectivity from dimension d down to dimension 0:
       * @f[
       *   d \to (d-1) \to \ldots \to 1 \to 0
       * @f]
       *
       * For example, building H1 on a 3D tetrahedral mesh requires:
       * - 3 → 2 (cells to faces)
       * - 2 → 1 (faces to edges)  
       * - 1 → 0 (edges to vertices)
       *
       * ### Example
       *
       * For a triangle (d=2) with 3 vertices and 3 edges:
       * 1. Recursively number the 3 edges (which number their 2 vertices each)
       * 2. Inject edge DOFs into triangle using Cochain<Triangle>::map
       * 3. Number any interior triangle DOFs (for K≥3)
       *
       * The result is stored in m_closure[d][idx], which maps local DOF indices
       * to global DOF indices for this polytope.
       *
       * @param[in] d Dimension of the polytope (0=vertex, 1=edge, 2=face, D=cell)
       * @param[in] idx Global index of the polytope in dimension d
       *
       * @note This method is called recursively and uses m_visited to avoid
       *       processing the same entity multiple times.
       *
       * @see Cochain::map for the boundary injection logic
       * @see H1() constructor which initiates the process for all cells
       */
      void getClosure(size_t d, Index idx);

      /**
       * @brief Gets the finite element associated with a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] i Index of the polytope
       * @return Reference to the H1 element for this polytope type
       *
       * Returns the appropriate H1 element based on the polytope geometry.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Point);
            return s_element;
          }
          case Geometry::Polytope::Type::Segment:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Segment);
            return s_element;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Triangle);
            return s_element;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Quadrilateral);
            return s_element;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Tetrahedron);
            return s_element;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Wedge);
            return s_element;
          }
          case Geometry::Polytope::Type::Hexahedron:
          {
            static thread_local const ElementType s_element(Geometry::Polytope::Type::Hexahedron);
            return s_element;
          }
        }
        assert(false);
        static thread_local const ElementType s_null;
        return s_null;
      }

      /**
       * @brief Gets the total number of degrees of freedom.
       * @return Number of DOFs
       *
       * For H1 spaces of degree K, the number of DOFs depends on the
       * polynomial degree and mesh topology.
       */
      size_t getSize() const override
      {
        return m_size;
      }

      /**
       * @brief Gets the vector dimension of the space.
       * @return Vector dimension (1 for scalar H1)
       *
       * Returns the number of components per DOF. For scalar H1, this is 1.
       */
      size_t getVectorDimension() const override
      {
        return 1;
      }

      /**
       * @brief Gets the underlying mesh.
       * @return Reference to the mesh
       */
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * @brief Gets the global DOF indices for a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] i Index of the polytope
       * @return Array of global DOF indices
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_closure[d][i];
      }

      /**
       * @brief Converts local to global DOF index.
       * @param[in] idx Pair of (dimension, polytope index)
       * @param[in] local Local DOF index within the polytope
       * @return Global DOF index
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        return m_closure[d][i](local);
      }

      /**
       * @brief Returns the Pullback of the function from the physical element
       * to the reference element.
       * @param[in] idx Index of the element in the mesh
       * @param[in] v Function defined on an element of the mesh
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /**
       * @brief Returns the inverse Pullback of the function from the physical
       * element to the reference element.
       * @param[in] idx Index of the element in the mesh.
       * @param[in] v Callable type
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;

      size_t m_size;
      std::vector<std::vector<uint8_t>> m_visited;
      std::vector<std::vector<IndexArray>> m_closure;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for H1 from mesh - deduces to RealH1
   */
  template <size_t K, class Context>
  H1(std::integral_constant<size_t, K>, const Geometry::Mesh<Context>&)
    -> H1<K, Real, Geometry::Mesh<Context>>;

  /// Alias for a scalar real-valued H1 finite element space
  template <size_t K, class Mesh>
  using RealH1 = H1<K, Real, Mesh>;

  /// Alias for a scalar complex-valued H1 finite element space
  template <size_t K, class Mesh>
  using ComplexH1 = H1<K, Complex, Mesh>;

  /**
   * @ingroup H1Specializations
   * @brief Vector-valued continuous piecewise polynomial (degree K) Lagrange finite element space.
   *
   * Represents the finite element space composed of @f$ d @f$-dimensional
   * vector-valued, continuous, piecewise polynomial functions:
   * @f[
   *  [\mathbb{P}_K (\mathcal{T}_h)]^d = \{ \mathbf{v} \in [C^0(\mathcal{T}_h)]^d : \mathbf{v}|_{\tau} \in [\mathbb{P}_K(\tau)]^d, \ \tau \in \mathcal{T}_h \}
   * @f]
   * where @f$ d @f$ is the vector dimension (typically the spatial dimension).
   *
   * @tparam K Polynomial degree
   * @tparam Scalar Scalar type for vector components (Real or Complex)
   */
  template <size_t K, class Scalar>
  class H1<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>,
        H1<K, Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    public:
      static_assert(K > 0, "Polynomial degree K must be greater than 0.");

      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = Math::Vector<ScalarType>;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<Context::Local>;

      /// Represents the Context of the H1 space
      using ContextType = Context::Local;

      /// Type of finite element
      using ElementType = H1Element<K, RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, H1<K, RangeType, MeshType>>;

      template <class Callable>
      class Pullback :
        public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          Pullback(const Pullback&) = default;

          decltype(auto) operator()(const Math::SpatialPoint& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      template <class Callable>
      class Pushforward :
        public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          using CallableType = Callable;

          template <class Function>
          Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          Pushforward(const Pushforward&) = default;

          constexpr
          decltype(auto) operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      /**
       * @brief Constructs a vector H1 space on the given mesh.
       * @param[in] mesh Mesh on which to build the finite element space
       * @param[in] vdim Vector dimension
       */
      H1(std::integral_constant<size_t, K>, const Geometry::Mesh<ContextType>& mesh, size_t vdim);

      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim),
          m_closure(other.m_dofs),
          m_size(other.m_size)
      {}

      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(std::move(other.m_vdim)),
          m_closure(std::move(other.m_dofs)),
          m_size(std::move(other.m_size))
      {}

      virtual ~H1() = default;

      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = std::move(other.m_vdim);
          m_closure = std::move(other.m_dofs);
          m_size = std::move(other.m_size);
        }
        return *this;
      }

      H1& operator=(const H1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_vdim = other.m_vdim;
          m_closure = other.m_dofs;
          m_size = other.m_size;
        }
        return *this;
      }

      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto& g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Point, 0),
              ElementType(Geometry::Polytope::Type::Point, 1),
              ElementType(Geometry::Polytope::Type::Point, 2),
              ElementType(Geometry::Polytope::Type::Point, 3)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Segment:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Segment, 0),
              ElementType(Geometry::Polytope::Type::Segment, 1),
              ElementType(Geometry::Polytope::Type::Segment, 2),
              ElementType(Geometry::Polytope::Type::Segment, 3)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Triangle, 0),
              ElementType(Geometry::Polytope::Type::Triangle, 1),
              ElementType(Geometry::Polytope::Type::Triangle, 2),
              ElementType(Geometry::Polytope::Type::Triangle, 3)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Quadrilateral, 0),
              ElementType(Geometry::Polytope::Type::Quadrilateral, 1),
              ElementType(Geometry::Polytope::Type::Quadrilateral, 2),
              ElementType(Geometry::Polytope::Type::Quadrilateral, 3)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Tetrahedron, 0),
              ElementType(Geometry::Polytope::Type::Tetrahedron, 1),
              ElementType(Geometry::Polytope::Type::Tetrahedron, 2),
              ElementType(Geometry::Polytope::Type::Tetrahedron, 3)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Wedge, 0),
              ElementType(Geometry::Polytope::Type::Wedge, 1),
              ElementType(Geometry::Polytope::Type::Wedge, 2),
              ElementType(Geometry::Polytope::Type::Wedge, 3)
            };
            return s_elements[m_vdim];
          }
          case Geometry::Polytope::Type::Hexahedron:
          {
            static thread_local std::array<ElementType, RODIN_MAXIMAL_SPACE_DIMENSION + 1> s_elements =
            {
              ElementType(Geometry::Polytope::Type::Hexahedron, 0),
              ElementType(Geometry::Polytope::Type::Hexahedron, 1),
              ElementType(Geometry::Polytope::Type::Hexahedron, 2),
              ElementType(Geometry::Polytope::Type::Hexahedron, 3)
            };
            return s_elements[m_vdim];
          }
        }
        assert(false);
        static thread_local ElementType s_null(Geometry::Polytope::Type::Point, 0);
        return s_null;
      }

      size_t getSize() const override
      {
        return m_size;
      }

      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return m_closure[d][i];
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        return m_closure[d][i](local);
      }

      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const Geometry::Mesh<ContextType>> m_mesh;
      size_t m_vdim;

      size_t m_size;
      std::vector<std::vector<IndexArray>> m_closure;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Vector H1 from mesh and vector dimension
   */
  template <size_t K, class Context>
  H1(std::integral_constant<size_t, K>, const Geometry::Mesh<Context>&, size_t)
    -> H1<K, Math::Vector<Real>, Geometry::Mesh<Context>>;

  /// Alias for a vector-valued real H1 finite element space
  template <size_t K, class Mesh>
  using VectorH1 = H1<K, Math::Vector<Real>, Mesh>;
}

#include "H1.hpp"

#endif
