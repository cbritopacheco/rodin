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

#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"

#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "H1Element.h"

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
   * @tparam K Polynomial degree (1, 2, 3, ...)
   * @tparam Scalar Scalar type (Real, Complex)
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
       */
      H1(std::integral_constant<size_t, K>, const MeshType& mesh)
        : m_mesh(mesh)
      {
        const auto& conn = mesh.getConnectivity();
        const size_t D   = mesh.getDimension();

        // Ensure we have incidence from cells to all lower dimensions
        for (size_t d = 0; d <= D; ++d)
          RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, d);

        m_owned.resize(D + 1);
        for (size_t d = 0; d <= D; ++d)
          m_owned[d].resize(conn.getCount(d));

        std::vector<std::vector<uint8_t>> visited(D + 1);
        for (size_t d = 0; d <= D; ++d)
          visited[d].assign(conn.getCount(d), false);

        Index global = 0;
        for (Index c = 0; c < conn.getCount(D); ++c)
        {
          // Loop over all subentity dimensions of the cell
          for (size_t dp = 0; dp <= D; ++dp)
          {
            const auto& inc = conn.getIncidence({ D, dp }, c);
            for (Index idx : inc)
            {
              if (visited[dp][idx])
                continue;

              visited[dp][idx] = true;

              size_t count = 0;
              const auto g = mesh.getGeometry(dp, idx);
              switch (g)
              {
                case Geometry::Polytope::Type::Point:
                {
                  // Vertices: 1 DOF per vertex for K >= 1
                  if constexpr (K >= 1)
                    count = 1;
                  break;
                }

                case Geometry::Polytope::Type::Segment:
                {
                  // Edges: (K-1) interior DOFs for K >= 2
                  if constexpr (K >= 2)
                    count = K - 1;
                  break;
                }

                case Geometry::Polytope::Type::Triangle:
                {
                  // Triangular faces / cells: interior DOFs (K-1)(K-2)/2 for K >= 3
                  if constexpr (K >= 3)
                  {
                    if (dp >= 2) // 2D cells or 3D faces
                      count = (K - 1) * (K - 2) / 2;
                  }
                  break;
                }

                case Geometry::Polytope::Type::Quadrilateral:
                {
                  // Quad faces / cells: interior DOFs (K-1)^2 for K >= 2
                  if constexpr (K >= 2)
                  {
                    if (dp >= 2)
                      count = (K - 1) * (K - 1);
                  }
                  break;
                }

                case Geometry::Polytope::Type::Tetrahedron:
                {
                  // Tetrahedra: interior DOFs (K-1)(K-2)(K-3)/6 for K >= 4
                  if constexpr (K >= 4)
                  {
                    if (dp == 3)
                      count = (K - 1) * (K - 2) * (K - 3) / 6;
                  }
                  break;
                }

                case Geometry::Polytope::Type::Wedge:
                {
                  // Wedges (prisms): choose a reasonable interior DOF formula
                  if constexpr (K >= 3)
                  {
                    if (dp == 3)
                      count = (K - 1) * (K - 1) * (K - 2) / 2;
                  }
                  break;
                }
              }

              m_owned[dp][idx].resize(count);
              for (size_t k = 0; k < count; ++k)
                m_owned[dp][idx](k) = global++;
            }
          }
        }

        m_size = static_cast<size_t>(global);

        m_closure.assign(D + 1, {});
        for (size_t d = 0; d <= D; ++d)
          m_closure[d].resize(conn.getCount(d));

        // For lower-dimensional entities, closure := owned (can be enriched later)
        for (size_t d = 0; d < D; ++d)
        {
          const size_t nd = conn.getCount(d);
          for (size_t i = 0; i < nd; ++i)
            m_closure[d][i] = m_owned[d][i];
        }

        // For cells (d = D), build standard H1 closure:
        //   vertices -> edges -> faces -> cell interior
        const size_t cellCount = conn.getCount(D);
        for (size_t i = 0; i < cellCount; ++i)
        {
          std::vector<Index> buffer;
          buffer.reserve(32); // small default; will grow if needed

          if (D == 1)
          {
            // 1D: segment cell
            //  - vertex DOFs
            const auto& poly = conn.getPolytope(1, i); // vertex indices
            for (Index v : poly)
            {
              for (Index dof : m_owned[0][v])
                buffer.push_back(dof);
            }
            //  - edge interior DOFs
            for (Index dof : m_owned[1][i])
              buffer.push_back(dof);
          }
          else if (D == 2)
          {
            // 2D: triangle / quad cell
            // 1) vertex DOFs
            const auto& poly = conn.getPolytope(2, i); // vertex indices
            for (Index v : poly)
            {
              for (Index dof : m_owned[0][v])
                buffer.push_back(dof);
            }

            // 2) edge interior DOFs
            const auto& edges = conn.getIncidence({ 2, 1 }, i);
            for (Index e : edges)
            {
              for (Index dof : m_owned[1][e])
                buffer.push_back(dof);
            }

            // 3) cell interior DOFs
            for (Index dof : m_owned[2][i])
              buffer.push_back(dof);
          }
          else if (D == 3)
          {
            // 3D: tetra / wedge cell
            // 1) vertex DOFs
            const auto& poly = conn.getPolytope(3, i); // vertex indices
            for (Index v : poly)
            {
              for (Index dof : m_owned[0][v])
                buffer.push_back(dof);
            }

            // 2) edge interior DOFs
            const auto& edges = conn.getIncidence({ 3, 1 }, i);
            for (Index e : edges)
            {
              for (Index dof : m_owned[1][e])
                buffer.push_back(dof);
            }

            // 3) face interior DOFs
            const auto& faces = conn.getIncidence({ 3, 2 }, i);
            for (Index f : faces)
            {
              for (Index dof : m_owned[2][f])
                buffer.push_back(dof);
            }

            // 4) cell interior DOFs
            for (Index dof : m_owned[3][i])
              buffer.push_back(dof);
          }
          else
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Dimension " << D << " not supported." << Alert::Raise;
          }

          auto& dofs = m_closure[D][i];
          dofs.resize(buffer.size());
          std::copy(buffer.begin(), buffer.end(), dofs.data());
        }
      }

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
      std::vector<std::vector<IndexArray>> m_owned;
      std::vector<std::vector<IndexArray>> m_closure;
      size_t m_size;
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
      H1(std::integral_constant<size_t, K>, const Geometry::Mesh<ContextType>& mesh, size_t vdim)
        : m_mesh(mesh), m_vdim(vdim)
      {
        const auto& conn = mesh.getConnectivity();
        const size_t D = mesh.getDimension();

        for (size_t d = 0; d <= D; ++d)
          RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, d);

        m_owned.resize(D + 1);
        for (size_t d = 0; d <= D; ++d)
          m_owned[d].resize(conn.getCount(d));

        std::vector<std::vector<uint8_t>> visited(D + 1);
        for (size_t d = 0; d <= D; ++d)
          visited[d].assign(conn.getCount(d), false);

        Index global = 0;
        for (Index c = 0; c < conn.getCount(D); ++c)
        {
          for (size_t dp = 0; dp <= D; ++dp)
          {
            const auto& inc = conn.getIncidence({ D, dp }, c);
            for (Index idx : inc)
            {
              if (visited[dp][idx])
                continue;

              visited[dp][idx] = true;

              // scalar DOFs owned by (dp, idx)
              size_t count = 0;
              const auto g = mesh.getGeometry(dp, idx);
              switch (g)
              {
                case Geometry::Polytope::Type::Point:
                {
                  if constexpr (K >= 1)
                    count = 1;
                  break;
                }

                case Geometry::Polytope::Type::Segment:
                {
                  if constexpr (K >= 2)
                    count = K - 1;
                  break;
                }

                case Geometry::Polytope::Type::Triangle:
                {
                  if constexpr (K >= 3)
                  {
                    if (dp >= 2)
                      count = (K - 1) * (K - 2) / 2;
                  }
                  break;
                }

                case Geometry::Polytope::Type::Quadrilateral:
                {
                  if constexpr (K >= 2)
                  {
                    if (dp >= 2)
                      count = (K - 1) * (K - 1);
                  }
                  break;
                }

                case Geometry::Polytope::Type::Tetrahedron:
                {
                  if constexpr (K >= 4)
                  {
                    if (dp == 3)
                      count = (K - 1) * (K - 2) * (K - 3) / 6;
                  }
                  break;
                }

                case Geometry::Polytope::Type::Wedge:
                {
                  if constexpr (K >= 3)
                  {
                    if (dp == 3)
                      count = (K - 1) * (K - 1) * (K - 2) / 2;
                  }
                  break;
                }
              }

              count *= m_vdim;
              m_owned[dp][idx].resize(count);
              for (size_t k = 0; k < count; ++k)
                m_owned[dp][idx](k) = global++;
            }
          }
        }

        m_size = static_cast<size_t>(global);

        m_closure.assign(D + 1, {});
        for (size_t d = 0; d <= D; ++d)
          m_closure[d].resize(conn.getCount(d));

        // For lower-dimensional entities, closure := owned (for now)
        for (size_t d = 0; d < D; ++d)
        {
          const size_t nd = conn.getCount(d);
          for (size_t i = 0; i < nd; ++i)
            m_closure[d][i] = m_owned[d][i];
        }

        // For cells (d = D), build closure: vertices + edges + faces + interior
        const size_t cellCount = conn.getCount(D);
        for (size_t i = 0; i < cellCount; ++i)
        {
          std::vector<Index> buffer;
          buffer.reserve(64);

          if (D == 1)
          {
            // 1D segment cell
            const auto& poly = conn.getPolytope(1, i);
            for (Index v : poly)
            {
              for (Index dof : m_owned[0][v])
                buffer.push_back(dof);
            }
            for (Index dof : m_owned[1][i])
              buffer.push_back(dof);
          }
          else if (D == 2)
          {
            // 2D triangle / quad cell
            // 1) vertex vector DOFs
            const auto& poly = conn.getPolytope(2, i);
            for (Index v : poly)
            {
              for (Index dof : m_owned[0][v])
                buffer.push_back(dof);
            }

            // 2) edge vector DOFs
            const auto& edges = conn.getIncidence({ 2, 1 }, i);
            for (Index e : edges)
            {
              for (Index dof : m_owned[1][e])
                buffer.push_back(dof);
            }

            // 3) cell-interior vector DOFs
            for (Index dof : m_owned[2][i])
              buffer.push_back(dof);
          }
          else if (D == 3)
          {
            // 3D tetra / wedge cell
            // 1) vertex DOFs
            const auto& poly = conn.getPolytope(3, i);
            for (Index v : poly)
            {
              for (Index dof : m_owned[0][v])
                buffer.push_back(dof);
            }

            // 2) edge DOFs
            const auto& edges = conn.getIncidence({ 3, 1 }, i);
            for (Index e : edges)
            {
              for (Index dof : m_owned[1][e])
                buffer.push_back(dof);
            }

            // 3) face DOFs
            const auto& faces = conn.getIncidence({ 3, 2 }, i);
            for (Index f : faces)
            {
              for (Index dof : m_owned[2][f])
                buffer.push_back(dof);
            }

            // 4) cell-interior DOFs
            for (Index dof : m_owned[3][i])
              buffer.push_back(dof);
          }
          else
          {
            Alert::MemberFunctionException(*this, __func__)
              << "Dimension " << D << " not supported." << Alert::Raise;
          }

          m_closure[D][i].resize(buffer.size());
          std::copy(buffer.begin(), buffer.end(), m_closure[D][i].data());
        }
      }

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
      std::vector<std::vector<IndexArray>> m_owned;
      std::vector<std::vector<IndexArray>> m_closure;
      size_t m_size;
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

#endif
