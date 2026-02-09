/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_P1_H
#define RODIN_VARIATIONAL_P1_P1_H

#include <boost/multi_array.hpp>
#include <functional>

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"

#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "P1Element.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Mesh>
  struct Traits<Variational::P1<Scalar, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = ScalarType;
    using ElementType = Variational::P1Element<RangeType>;
  };

  template <class Scalar, class Mesh>
  struct Traits<Variational::P1<Math::Vector<Scalar>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Scalar;
    using RangeType = Math::Vector<ScalarType>;
    using ElementType = Variational::P1Element<RangeType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P1Specializations P1 Template Specializations
   * @brief Template specializations of the P1 class.
   * @see P1
   */

  template <class Range, class Mesh = Geometry::Mesh<Context::Local>>
  class P1;

  /**
   * @ingroup P1Specializations
   * @brief Real valued Lagrange finite element space
   *
   * Represents the finite element space composed of scalar valued continuous,
   * piecewise linear functions:
   * @f[
   *  \mathbb{P}_1 (\mathcal{T}_h) = \{ v \in C^0(\mathcal{T}_h) : v|_{\tau} \in \mathbb{P}_1(\tau), \ \tau \in \mathcal{T}_h \} \ .
   * @f]
   *
   * This class is scalar valued, i.e. evaluations of the function are of
   * Rodin::Real type.
   */
  template <class Scalar>
  class P1<Scalar, Geometry::Mesh<Context::Local>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, P1<Scalar, Geometry::Mesh<Context::Local>>>
  {
    public:
      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the P1 space
      using ContextType = Context::Local;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

      /**
       * @brief Pullback for the scalar/complex P1 space.
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
       * @brief Inverse Pullback for the scalar/complex P1 space.
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
       * @brief Constructs a P1 finite element space on the given mesh.
       * @param[in] mesh Mesh on which to build the finite element space
       *
       * Creates the P1 space with one degree of freedom per vertex. The total
       * number of DOFs equals the number of mesh vertices.
       */
      P1(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other P1 space to copy
       */
      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other P1 space to move from
       */
      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh))
      {}

      virtual ~P1() = default;

      /**
       * @brief Move assignment operator.
       * @param[in] other P1 space to move from
       * @return Reference to this P1 space
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
       * @brief Copy assignment operator.
       * @param[in] other P1 space to copy
       * @return Reference to this P1 space
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
       * @brief Gets the finite element associated with a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] i Index of the polytope
       * @return Reference to the P1 element for this polytope type
       *
       * Returns the appropriate P1 element based on the polytope geometry
       * (point, segment, triangle, quadrilateral, tetrahedron, or wedge).
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Point);
            return s_element;
          }
          case Geometry::Polytope::Type::Segment:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Segment);
            return s_element;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Triangle);
            return s_element;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Quadrilateral);
            return s_element;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Tetrahedron);
            return s_element;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Wedge);
            return s_element;
          }
          case Geometry::Polytope::Type::Hexahedron:
          {
            static thread_local constexpr ElementType s_element(Geometry::Polytope::Type::Hexahedron);
            return s_element;
          }
        }
        assert(false);
        static thread_local constexpr ElementType s_null;
        return s_null;
      }

      /**
       * @brief Gets the total number of degrees of freedom.
       * @return Number of DOFs (equals number of vertices)
       *
       * For P1 spaces, the number of DOFs equals the number of mesh vertices
       * since each vertex has one DOF.
       */
      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount();
      }

      /**
       * @brief Gets the vector dimension of the space.
       * @return Vector dimension (1 for scalar P1)
       *
       * Returns the number of components per DOF. For scalar P1, this is 1.
       * Vector-valued P1 spaces have dimension equal to the space dimension.
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
       *
       * For P1, returns the vertex indices of the polytope, which are
       * the global DOF indices.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        return getMesh().getConnectivity().getPolytope(d, i);
      }

      /**
       * @brief Converts local to global DOF index.
       * @param[in] idx Pair of (dimension, polytope index)
       * @param[in] local Local DOF index within the polytope
       * @return Global DOF index
       *
       * Maps the local DOF index (vertex number within element) to the
       * global DOF index (global vertex number).
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        assert(p.size() >= 0);
        assert(local < static_cast<size_t>(p.size()));
        return p(local);
      }

      /**
       * @brief Returns the Pullback of the function from the physical element
       * to the reference element.
       * @param[in] idx Index of the element in the mesh
       * @param[in] v Function defined on an element of the mesh
       *
       * For all @f$ \tau \in \mathcal{T}_h @f$ in the mesh, the finite
       * element space is generated by the bijective Pullback:
       * @f[
       *  \psi_\tau : \mathbb{P}_1(\tau) \rightarrow \mathbb{P}_1(R)
       * @f]
       * taking a function @f$ v \in V(\tau) @f$ from the global element @f$
       * \tau @f$ element @f$ R @f$.
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
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for P1 from mesh - deduces to RealP1
   */
  template <class Context>
  P1(const Geometry::Mesh<Context>&) -> P1<Real, Geometry::Mesh<Context>>;

  /// Alias for a scalar real-valued P1 finite element space
  template <class Mesh>
  using RealP1 = P1<Real, Mesh>;

  /// Alias for a scalar complex-valued P1 finite element space
  template <class Mesh>
  using ComplexP1 = P1<Complex, Mesh>;

  /**
   * @ingroup P1Specializations
   * @brief Vector-valued continuous piecewise linear Lagrange finite element space.
   *
   * Represents the finite element space composed of @f$ d @f$-dimensional
   * vector-valued, continuous, piecewise linear functions:
   * @f[
   *  [\mathbb{P}_1 (\mathcal{T}_h)]^d = \{ \mathbf{v} \in [C^0(\mathcal{T}_h)]^d : \mathbf{v}|_{\tau} \in [\mathbb{P}_1(\tau)]^d, \ \tau \in \mathcal{T}_h \}
   * @f]
   * where @f$ d @f$ is the vector dimension (typically the spatial dimension).
   *
   * ## Properties
   * - **DOF count**: @f$ d \cdot n_{\text{vertices}} @f$ total
   * - **DOF structure**: Each vertex has @f$ d @f$ DOFs (one per vector component)
   * - **Continuity**: C⁰ continuous for each component
   * - **Basis functions**: @f$ \boldsymbol{\phi}_{i,j}(x) = \phi_i(x) \mathbf{e}_j @f$
   *
   * ## Use Cases
   * - Displacement fields in elasticity
   * - Velocity fields in fluid mechanics
   * - Vector-valued PDEs requiring H¹ conformity
   *
   * ## Example
   * @code{.cpp}
   * Mesh Th;
   * Th = Th.UniformGrid(Polytope::Type::Triangle, {8, 8});
   * size_t vdim = 2;  // 2D displacement
   * P1 Vh(Th, vdim);  // Vector P1 space
   * GridFunction u(Vh);
   * @endcode
   *
   * @tparam Scalar Scalar type for vector components (Real or Complex)
   *
   * @see P1Element, TrialFunction, TestFunction
   */
  template <class Scalar>
  class P1<Math::Vector<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>,
        P1<Math::Vector<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = Math::Vector<ScalarType>;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<Context::Local>;

      /// Represents the Context of the P1 space
      using ContextType = Context::Local;

      /// Type of finite element
      using ElementType = P1Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P1<RangeType, MeshType>>;

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

      P1(const Geometry::Mesh<ContextType>& mesh, size_t vdim)
        : m_mesh(mesh), m_vdim(vdim)
      {
        const size_t vn = mesh.getVertexCount();
        m_dofs.resize(mesh.getDimension() + 1);
        for (size_t d = 0; d <= mesh.getDimension(); d++)
        {
          const size_t n = mesh.getConnectivity().getCount(d);
          m_dofs[d].reserve(n);
          for (size_t i = 0; i < n; i++)
          {
            const auto& polytope = mesh.getConnectivity().getPolytope(d, i);
            const size_t count = polytope.size();
            auto& dofs = m_dofs[d].emplace_back(count * vdim);
            for (size_t local = 0; local < count * vdim; local++)
            {
              const size_t q = local / vdim;
              const size_t r = local % vdim;
              assert(q < count);
              dofs.coeffRef(local) = polytope(q) + r * vn;
            }
          }
        }
      }

      P1(const P1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim),
          m_dofs(other.m_dofs)
      {}

      P1(P1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(std::move(other.m_vdim)),
          m_dofs(std::move(other.m_dofs))
      {}

      virtual ~P1() = default;

      P1& operator=(P1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = std::move(other.m_vdim);
          m_dofs = std::move(other.m_dofs);
        }
        return *this;
      }

      P1& operator=(const P1& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
          m_vdim = other.m_vdim;
          m_dofs = other.m_dofs;
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
        return m_mesh.get().getVertexCount() * m_vdim;
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
        return m_dofs[d][i];
      }

      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto& [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        const size_t q = local / m_vdim;
        const size_t r = local % m_vdim;
        assert(q < static_cast<size_t>(p.size()));
        return p(q) + r * m_mesh.get().getVertexCount();
      }

      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>&, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const Geometry::Mesh<ContextType>> m_mesh;
      size_t m_vdim;
      std::vector<std::vector<IndexArray>> m_dofs;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Vector P1 from mesh and vector dimension
   */
  template <class Context>
  P1(const Geometry::Mesh<Context>&, size_t)
    -> P1<Math::Vector<Real>, Geometry::Mesh<Context>>;

  /// Alias for a vector-valued real P1 finite element space
  template <class Mesh>
  using VectorP1 = P1<Math::Vector<Real>, Mesh>;
}

#endif
