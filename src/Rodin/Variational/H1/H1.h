/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_H1_H1_H
#define RODIN_VARIATIONAL_H1_H1_H

#include <functional>
#include <type_traits>

#include <boost/multi_array.hpp>

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
   * @tparam K Polynomial degree (0, 1, 2, 3, ...)
   * @tparam Scalar Scalar type (Real, Complex)
   */
  template <size_t K, class Scalar>
  class H1<K, Scalar, Geometry::Mesh<Context::Local>>
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, H1<K, Scalar, Geometry::Mesh<Context::Local>>>
  {
    public:
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
        buildDOFs();
      }

      /**
       * @brief Copy constructor.
       * @param[in] other H1 space to copy
       */
      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_dofs(other.m_dofs),
          m_size(other.m_size)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other H1 space to move from
       */
      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_dofs(std::move(other.m_dofs)),
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
          m_dofs = std::move(other.m_dofs);
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
          m_dofs = other.m_dofs;
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
        return m_dofs[d][i];
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
        return m_dofs[d][i](local);
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
      /**
       * @brief Build the DOF structure for the H1 space.
       *
       * For H1-conforming elements, DOFs are distributed across vertices,
       * edges, faces, and interior based on the polynomial degree K.
       *
       * DOF distribution for degree K:
       * - K = 0: No DOFs (constant elements, similar to P0)
       * - K >= 1: 1 DOF per vertex
       * - K >= 2: (K-1) DOFs per edge interior
       * - K >= 3: (K-1)(K-2)/2 DOFs per triangle interior (for triangles)
       *           (K-1)^2 DOFs per quadrilateral interior (for quads)
       * - K >= 4: (K-1)(K-2)(K-3)/6 DOFs per tetrahedron interior
       */
      void buildDOFs()
      {
        const auto& mesh = getMesh();
        const size_t meshDim = mesh.getDimension();

        // Initialize DOF storage for each dimension
        m_dofs.resize(meshDim + 1);

        // For H1<K>, DOFs are placed at:
        // - Vertices: 1 DOF per vertex (for K >= 1)
        // - Edges: (K-1) DOFs per edge (for K >= 2)
        // - Faces (triangles): (K-1)(K-2)/2 DOFs per face (for K >= 3)
        // - Faces (quads): (K-1)^2 DOFs per face (for K >= 2)
        // - Cells: interior DOFs

        // Count DOFs by entity type
        Index globalDOFIndex = 0;

        // Vertex DOFs (dimension 0)
        // Each vertex has 1 DOF if K >= 1, 0 DOFs if K == 0
        // For K=0, this results in a space with no DOFs (similar to P0 behavior)
        const size_t vertexCount = mesh.getVertexCount();
        const size_t dofsPerVertex = (K >= 1) ? 1 : 0;
        m_dofs[0].resize(vertexCount);
        for (size_t v = 0; v < vertexCount; ++v)
        {
          if (dofsPerVertex > 0)
          {
            m_dofs[0][v].resize(dofsPerVertex);
            m_dofs[0][v](0) = globalDOFIndex++;
          }
          else
          {
            m_dofs[0][v].resize(0);
          }
        }

        // Edge DOFs (dimension 1, if mesh dimension >= 1)
        if (meshDim >= 1)
        {
          const size_t edgeCount = mesh.getConnectivity().getCount(1);
          const size_t dofsPerEdge = (K >= 2) ? (K - 1) : 0;
          m_dofs[1].resize(edgeCount);
          for (size_t e = 0; e < edgeCount; ++e)
          {
            m_dofs[1][e].resize(dofsPerEdge);
            for (size_t local = 0; local < dofsPerEdge; ++local)
              m_dofs[1][e](local) = globalDOFIndex++;
          }
        }

        // Face DOFs (dimension 2, if mesh dimension >= 2)
        if (meshDim >= 2)
        {
          const size_t faceCount = mesh.getConnectivity().getCount(2);
          m_dofs[2].resize(faceCount);
          for (size_t f = 0; f < faceCount; ++f)
          {
            const auto g = mesh.getGeometry(2, f);
            size_t dofsPerFace = 0;
            if (g == Geometry::Polytope::Type::Triangle)
            {
              // Interior DOFs for triangle: (K-1)(K-2)/2 for K >= 3
              dofsPerFace = (K >= 3) ? ((K - 1) * (K - 2) / 2) : 0;
            }
            else if (g == Geometry::Polytope::Type::Quadrilateral)
            {
              // Interior DOFs for quadrilateral: (K-1)^2 for K >= 2
              dofsPerFace = (K >= 2) ? ((K - 1) * (K - 1)) : 0;
            }
            m_dofs[2][f].resize(dofsPerFace);
            for (size_t local = 0; local < dofsPerFace; ++local)
              m_dofs[2][f](local) = globalDOFIndex++;
          }
        }

        // Cell DOFs (dimension meshDim, if mesh dimension >= 3)
        if (meshDim >= 3)
        {
          const size_t cellCount = mesh.getCellCount();
          m_dofs[3].resize(cellCount);
          for (size_t c = 0; c < cellCount; ++c)
          {
            const auto g = mesh.getGeometry(meshDim, c);
            size_t dofsPerCell = 0;
            if (g == Geometry::Polytope::Type::Tetrahedron)
            {
              // Interior DOFs for tetrahedron: (K-1)(K-2)(K-3)/6 for K >= 4
              dofsPerCell = (K >= 4) ? ((K - 1) * (K - 2) * (K - 3) / 6) : 0;
            }
            else if (g == Geometry::Polytope::Type::Wedge)
            {
              // Interior DOFs for wedge: more complex formula
              // For simplicity, use (K-1)*(K-1)*(K-2)/2 for K >= 3
              dofsPerCell = (K >= 3) ? ((K - 1) * (K - 1) * (K - 2) / 2) : 0;
            }
            m_dofs[3][c].resize(dofsPerCell);
            for (size_t local = 0; local < dofsPerCell; ++local)
              m_dofs[3][c](local) = globalDOFIndex++;
          }
        }

        m_size = globalDOFIndex;

        // Now rebuild the DOF arrays for each cell to include all DOFs
        // (vertex + edge + face + interior)
        rebuildCellDOFs();
      }

      /**
       * @brief Rebuild the complete DOF arrays for each cell.
       *
       * Each cell's DOF array should contain all DOFs associated with that cell,
       * including vertex, edge, face, and interior DOFs.
       *
       * @note For 2D meshes, m_dofs[2] initially stores only interior DOFs for each
       * face. This function reads those interior DOFs, builds a complete DOF array
       * including vertex and edge DOFs, then overwrites m_dofs[2] with the complete
       * array. This is safe because interior DOFs are unique to each cell (not shared).
       */
      void rebuildCellDOFs()
      {
        const auto& mesh = getMesh();
        const size_t meshDim = mesh.getDimension();
        const size_t cellCount = mesh.getCellCount();

        // We need to build complete DOF arrays for cells
        // The getDOFs(d, i) should return all DOFs for the d-polytope i
        // For cells, this includes vertex, edge, face, and interior DOFs

        for (size_t cellIdx = 0; cellIdx < cellCount; ++cellIdx)
        {
          const auto& fe = getFiniteElement(meshDim, cellIdx);
          const size_t numDOFs = fe.getCount();

          // Collect all DOFs for this cell
          IndexArray allDOFs(numDOFs);
          size_t dofIdx = 0;

          // Vertex DOFs
          const auto& cellVertices = mesh.getConnectivity().getPolytope(meshDim, cellIdx);
          const size_t dofsPerVertex = (K >= 1) ? 1 : 0;
          for (size_t v = 0; v < static_cast<size_t>(cellVertices.size()); ++v)
          {
            for (size_t local = 0; local < dofsPerVertex; ++local)
            {
              allDOFs(dofIdx++) = m_dofs[0][cellVertices(v)](local);
            }
          }

          // Edge DOFs (for meshDim >= 1 and K >= 2)
          if (meshDim >= 1 && K >= 2)
          {
            const size_t dofsPerEdge = K - 1;
            // Get edges of this cell
            const auto& cellEdges = mesh.getConnectivity().getIncidence({meshDim, 1}, cellIdx);
            for (size_t e = 0; e < cellEdges.size(); ++e)
            {
              for (size_t local = 0; local < dofsPerEdge; ++local)
              {
                allDOFs(dofIdx++) = m_dofs[1][cellEdges[e]](local);
              }
            }
          }

          // Face DOFs (for meshDim >= 3 and appropriate K)
          if (meshDim >= 3)
          {
            const auto& cellFaces = mesh.getConnectivity().getIncidence({meshDim, 2}, cellIdx);
            for (size_t f = 0; f < cellFaces.size(); ++f)
            {
              const size_t faceDofsCount = m_dofs[2][cellFaces[f]].size();
              for (size_t local = 0; local < faceDofsCount; ++local)
              {
                allDOFs(dofIdx++) = m_dofs[2][cellFaces[f]](local);
              }
            }
          }
          else if (meshDim == 2)
          {
            // For 2D meshes, faces are the cells themselves
            // Interior DOFs are already handled
            const size_t faceDofsCount = m_dofs[2][cellIdx].size();
            for (size_t local = 0; local < faceDofsCount; ++local)
            {
              allDOFs(dofIdx++) = m_dofs[2][cellIdx](local);
            }
          }

          // Interior cell DOFs (for meshDim >= 3)
          if (meshDim >= 3)
          {
            const size_t interiorDofsCount = m_dofs[meshDim][cellIdx].size();
            for (size_t local = 0; local < interiorDofsCount; ++local)
            {
              allDOFs(dofIdx++) = m_dofs[meshDim][cellIdx](local);
            }
          }

          // Store the complete DOF array for this cell
          m_dofs[meshDim][cellIdx] = allDOFs;
        }
      }

      std::reference_wrapper<const MeshType> m_mesh;
      std::vector<std::vector<IndexArray>> m_dofs;
      size_t m_size;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for H1 from mesh - deduces to RealH1
   */
  template <size_t K, class Context>
  H1(std::integral_constant<size_t, K>, const Geometry::Mesh<Context>&) -> H1<K, Real, Geometry::Mesh<Context>>;

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
      H1(const Geometry::Mesh<ContextType>& mesh, size_t vdim)
        : m_mesh(mesh), m_vdim(vdim)
      {
        // First build scalar DOFs
        H1<K, ScalarType, MeshType> scalarSpace(mesh);
        const size_t scalarSize = scalarSpace.getSize();

        const size_t meshDim = mesh.getDimension();
        m_dofs.resize(meshDim + 1);

        // Build vector DOFs from scalar DOFs
        for (size_t d = 0; d <= meshDim; d++)
        {
          const size_t count = mesh.getConnectivity().getCount(d);
          m_dofs[d].reserve(count);
          for (size_t i = 0; i < count; i++)
          {
            const auto& scalarDOFs = scalarSpace.getDOFs(d, i);
            const size_t scalarCount = scalarDOFs.size();
            auto& dofs = m_dofs[d].emplace_back(scalarCount * vdim);
            for (size_t local = 0; local < scalarCount * vdim; local++)
            {
              const size_t q = local / vdim;
              const size_t r = local % vdim;
              assert(q < scalarCount);
              dofs(local) = scalarDOFs(q) + r * scalarSize;
            }
          }
        }

        m_scalarSize = scalarSize;
      }

      H1(const H1& other)
        : Parent(other),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim),
          m_dofs(other.m_dofs),
          m_scalarSize(other.m_scalarSize)
      {}

      H1(H1&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(std::move(other.m_vdim)),
          m_dofs(std::move(other.m_dofs)),
          m_scalarSize(std::move(other.m_scalarSize))
      {}

      virtual ~H1() = default;

      H1& operator=(H1&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
          m_vdim = std::move(other.m_vdim);
          m_dofs = std::move(other.m_dofs);
          m_scalarSize = std::move(other.m_scalarSize);
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
          m_dofs = other.m_dofs;
          m_scalarSize = other.m_scalarSize;
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
        return m_scalarSize * m_vdim;
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
        return m_dofs[d][i](local);
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
      std::vector<std::vector<IndexArray>> m_dofs;
      size_t m_scalarSize;
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
