/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0_H
#define RODIN_VARIATIONAL_P0_P0_H

#include <map>
#include <boost/multi_array.hpp>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "P0Element.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type traits for @c P0: exposes the mesh type, the scalar type, the range
   * type, the execution context and the finite element type.
   */
  template <class Number, class Mesh>
  struct Traits<Variational::P0<Number, Mesh>>
  {
      /// @brief Mesh type.
      using MeshType = Mesh;
      /// @brief Scalar value type.
      using ScalarType = Number;
      /// @brief Range (evaluation value) type.
      using RangeType = ScalarType;
      /// @brief Execution context type.
      using ContextType = typename FormLanguage::Traits<MeshType>::ContextType;
      /// @brief Finite element type.
      using ElementType = Variational::P0Element<RangeType>;
  };

  /**
   * @brief Type traits for @c P0: exposes the mesh type, the scalar type, the range
   * type, the execution context and the finite element type.
   */
  template <class Number, class Mesh>
  struct Traits<Variational::P0<Math::SpatialVector<Number>, Mesh>>
  {
      /// @brief Mesh type.
      using MeshType = Mesh;
      /// @brief Scalar value type.
      using ScalarType = Number;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /// @brief Execution context type.
      using ContextType = typename FormLanguage::Traits<MeshType>::ContextType;
      /// @brief Finite element type.
      using ElementType = Variational::P0Element<Math::SpatialVector<ScalarType>>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P0Specializations P0 Template Specializations
   * @brief Template specializations of the P0 class.
   * @see <a href="class_rodin_1_1_variational_1_1_p0.html">P0</a>
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref P0 "P0<SpatialMatrix<Scalar>, Mesh<Context::Local>>" | Full rectangular matrix range with interleaved row-major components. |
   * | @ref P0 "P0<SpatialMatrix<Scalar>, Mesh<Context::MPI>>" | Distributed full matrix range with expanded scalar ownership and ghost maps. |
   * | @ref P0 "P0<Range, Mesh<Context::Local>>" | Real or complex scalar/vector local-mesh discontinuous piecewise constant space. |
   * | @ref P0 "P0<Range, Mesh<Context::MPI>>" | Scalar or vector-valued distributed-mesh discontinuous piecewise constant space. |
   */

  template <class Range, class Mesh = Geometry::Mesh<Context::Local>>
  class P0;

  /**
   * @ingroup P0Specializations
   * @brief Scalar-valued piecewise constant (P0) Lagrange finite element space.
   *
   * Represents the finite element space composed of discontinuous,
   * piecewise constant functions:
   * @f[
   *  \mathbb{P}_0 (\mathcal{T}_h) = \{ v \in L^2(\Omega) : v|_{\tau} \in \mathbb{P}_0(\tau), \ \tau \in \mathcal{T}_h \}
   * @f]
   * where @f$ \mathbb{P}_0(\tau) @f$ denotes constant functions on element @f$ \tau @f$.
   *
   * ## Properties
   * - **DOF count**: One DOF per cell (total: @f$ n_{\text{cells}} @f$)
   * - **DOF location**: Element barycenter
   * - **Continuity**: Discontinuous across element boundaries (@f$ L^2 @f$-conforming)
   * - **Polynomial degree**: 0 (constant functions)
   * - **Gradient**: @f$ \nabla u|_K = 0 @f$ (zero within each element)
   *
   * ## Use Cases
   * - Discontinuous Galerkin (DG) methods
   * - Flux-based formulations
   * - Element-wise constant material properties
   * - Finite volume schemes as FEM
   * - Pressure space in mixed methods (when inf-sup stability permits)
   *
   * ## Example
   * @code{.cpp}
   * Mesh Th;
   * Th = Th.UniformGrid(Polytope::Type::Triangle, {8, 8});
   * P0 Vh(Th);  // Scalar P0 space
   * GridFunction u(Vh);
   * u = 1.0;  // Constant function
   * @endcode
   *
   * @see P0Element, GridFunction
   */
  template <class Scalar>
    requires(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>)
  class P0<Scalar, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Local>,
        P0<Scalar, Geometry::Mesh<Context::Local>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;

      /// Range type of value
      using RangeType = ScalarType;

      /// Represents the Context of the P0 space
      using ContextType = Context::Local;

      /// Type of mesh on which the finite element space is built
      using MeshType = Geometry::Mesh<ContextType>;

      /// Type of finite element
      using ElementType = P0Element<RangeType>;

      /// Parent class
      using Parent = FiniteElementSpace<MeshType, P0<RangeType, MeshType>>;

      /// @brief Pullback of a P0 function to the reference element.
      template <class Callable>
      class Pullback :
        public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          /// @brief Callable type evaluated on physical points.
          using CallableType = Callable;

          /**
           * @brief Constructs the pullback of a function on a polytope.
           * @param polytope Mesh entity used by this operation.
           * @param v Function operand.
           */
          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& v)
            : m_polytope(polytope), m_v(std::forward<Function>(v))
          {}

          /**
           * @brief Copy constructor.
           * @param other Object to copy from.
           */
          Pullback(const Pullback& other) = default;

          /**
           * @brief Evaluates at a point on the reference element.
           * @param r Reference coordinates at which to evaluate the basis.
           * @returns Value of the expression at the supplied evaluation point.
           */
          auto operator()(const Math::SpatialPoint& r) const
          {
            const Geometry::Point p(m_polytope, r);
            return m_v(p);
          }

        private:
          Geometry::Polytope m_polytope;
          CallableType m_v;
      };

      /// @brief Pushforward of a P0 function to the physical element.
      template <class Callable>
      class Pushforward :
        public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          /// @brief Callable type evaluated on physical points.
          using CallableType = Callable;

          /**
           * @param[in] v Reference to the function defined on the reference
           * space.
           */
          template <class Function>
          Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          /**
           * @brief Copy constructor.
           * @param other Object to copy from.
           */
          Pushforward(const Pushforward& other) = default;

          /**
           * @brief Evaluates at a geometric point.
           * @param p Point at which the operation is evaluated.
           * @returns Value of the expression at the supplied evaluation point.
           */
          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      /**
       * @brief Constructs a P0 finite element space on the given mesh.
       * @param[in] mesh Mesh on which to build the finite element space
       *
       * Creates the P0 space with one degree of freedom per element. The total
       * number of DOFs equals the number of mesh cells. Each DOF represents
       * a constant value over the entire element.
       */
      P0(const MeshType& mesh)
        : m_mesh(mesh)
      {
        const size_t n = mesh.getCellCount();
        m_dofs.reserve(n);
        for (size_t i = 0; i < n; ++i)
          m_dofs.push_back(IndexArray{{i}});
      }

      /**
       * @brief Copy constructor.
       * @param[in] other P0 space to copy
       */
      P0(const P0& other)
        : Parent(other),
          m_dofs(other.m_dofs),
          m_mesh(other.m_mesh)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other P0 space to move from
       */
      P0(P0&& other)
        : Parent(std::move(other)),
          m_dofs(std::move(other.m_dofs)),
          m_mesh(other.m_mesh)
      {}

      virtual ~P0() override = default;

      /**
       * @brief Move assignment operator.
       * @param[in] other P0 space to move from
       * @return Reference to this P0 space
       */
      P0& operator=(P0&& other) = default;
      /// @brief Assigns the scalar basis and matrix component dimensions.
      P0& operator=(const P0& other) = default;

      /**
       * @brief Gets the finite element associated with a polytope.
       * @param[in] d Dimension of the polytope
       * @param[in] i Index of the polytope
       * @return Reference to the P0 element for this polytope type
       *
       * Returns the appropriate P0 element based on the polytope geometry.
       * P0 elements are piecewise constant over each element.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        switch (getMesh().getGeometry(d, i))
        {
          case Geometry::Polytope::Type::Point:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Point);
            return sElement;
          }
          case Geometry::Polytope::Type::Segment:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Segment);
            return sElement;
          }
          case Geometry::Polytope::Type::Triangle:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Triangle);
            return sElement;
          }
          case Geometry::Polytope::Type::Quadrilateral:
          {
            static constexpr ElementType sElement(
              Geometry::Polytope::Type::Quadrilateral);
            return sElement;
          }
          case Geometry::Polytope::Type::Tetrahedron:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Tetrahedron);
            return sElement;
          }
          case Geometry::Polytope::Type::Pyramid:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Pyramid);
            return sElement;
          }
          case Geometry::Polytope::Type::Wedge:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Wedge);
            return sElement;
          }
          case Geometry::Polytope::Type::Hexahedron:
          {
            static constexpr ElementType sElement(Geometry::Polytope::Type::Hexahedron);
            return sElement;
          }
        }
        assert(false);
        static constexpr ElementType sNull(Geometry::Polytope::Type::Point);
        return sNull;
      }

      /**
       * @brief Gets the total number of degrees of freedom.
       * @return Number of DOFs (equals number of cells)
       *
       * For P0 spaces, the number of DOFs equals the number of mesh cells
       * since each element has one constant DOF.
       */
      size_t getSize() const override
      {
        return m_mesh.get().getCellCount();
      }

      /**
       * @brief Gets the vector dimension of the space.
       * @return Vector dimension (1 for scalar P0)
       *
       * Returns the number of components per DOF. For scalar P0, this is 1.
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
       * @param[in] d Dimension of the polytope (must equal mesh dimension)
       * @param[in] i Index of the polytope (element index)
       * @return Array containing the single DOF index for this element
       *
       * For P0, each cell has exactly one DOF. The polytope must be
       * a top-dimensional cell.
       */
      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        assert(d == getMesh().getDimension());
        (void) d;
        return m_dofs.at(i);
      }

      /**
       * @brief Converts local to global DOF index.
       * @param[in] idx Pair of (dimension, cell index)
       * @param[in] local Local DOF index (always 0 for P0)
       * @return Global DOF index (equals cell index)
       *
       * For P0, the global DOF index equals the cell index since there
       * is one DOF per cell.
       */
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        assert(d == getMesh().getDimension());
        return i;
      }

      /**
       * @brief Creates a pullback transformation for a function.
       * @tparam Callable Type of the callable function
       * @param[in] idx Pair of (dimension, polytope index)
       * @param[in] v Function to pull back
       * @return Pullback transformation object
       *
       * The pullback maps a function from physical space to reference space.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /**
       * @brief Creates a pushforward transformation for a function.
       * @tparam Callable Type of the callable function
       * @param[in] v Function to push forward
       * @return Pushforward transformation object
       *
       * The pushforward maps a function from reference space to physical space.
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>&, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      std::vector<IndexArray> m_dofs;
      std::reference_wrapper<const MeshType> m_mesh;
  };

  /**
   * @brief Cellwise constant vector space with independent component DOFs.
   *
   * Each cell owns @f$m@f$ DOFs for @f$[\mathbb P_0]^m@f$. The component
   * count is independent of the physical mesh dimension.
   */
  template <class Scalar>
    requires(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>)
  class P0<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Local>,
        P0<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    public:
      /// Scalar type of each vector component.
      using ScalarType = Scalar;
      /// Vector-valued range of the space.
      using RangeType = Math::SpatialVector<Scalar>;
      /// Local mesh context.
      using ContextType = Context::Local;
      /// Mesh supporting the cellwise constant field.
      using MeshType = Geometry::Mesh<ContextType>;
      /// Vector constant finite element.
      using ElementType = P0Element<RangeType>;
      /// Common finite element space interface.
      using Parent = FiniteElementSpace<MeshType, P0<RangeType, MeshType>>;

      /// Pulls a physical vector field back to a cell reference domain.
      template <class Callable>
      class Pullback : public FiniteElementSpacePullbackBase<Pullback<Callable>>
      {
        public:
          /**
           * Binds the physical cell and callable field.
           * @param polytope Mesh entity used by this operation.
           * @param function Function to evaluate.
           */
          template <class Function>
          Pullback(const Geometry::Polytope& polytope, Function&& function)
            : m_polytope(polytope),
              m_function(std::forward<Function>(function))
          {}

          /**
           * Evaluates the physical field at a reference coordinate.
           * @param reference Point at which the operation is evaluated.
           * @returns Reference to the entry at the supplied indices.
           */
          auto operator()(const Math::SpatialPoint& reference) const
          {
            return m_function(Geometry::Point(m_polytope, reference));
          }

        private:
          Geometry::Polytope m_polytope;
          Callable m_function;
      };

      /// Pushes a reference vector field forward to a physical cell.
      template <class Callable>
      class Pushforward : public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          /**
           * Binds a reference-domain callable.
           * @param function Function to evaluate.
           */
          template <class Function>
          explicit Pushforward(Function&& function)
            : m_function(std::forward<Function>(function))
          {}

          /**
           * Evaluates the reference field at a physical point's chart coordinate.
           * @param point Point at which the operation is evaluated.
           * @returns Reference to the entry at the supplied indices.
           */
          auto operator()(const Geometry::Point& point) const
          {
            return m_function(point.getReferenceCoordinates());
          }

        private:
          Callable m_function;
      };

      /**
       * Constructs a cellwise vector space with @p vdim components.
       * @param mesh Mesh on which the object is defined.
       * @param vdim Number of components in the value range.
       */
      explicit P0(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh),
          m_vdim(vdim)
      {
        assert(m_vdim > 0);
        m_dofs.reserve(mesh.getCellCount());
        for (size_t cell = 0; cell < mesh.getCellCount(); ++cell)
        {
          IndexArray dofs(m_vdim);
          for (size_t component = 0; component < m_vdim; ++component)
            dofs[component] = cell * m_vdim + component;
          m_dofs.push_back(std::move(dofs));
        }
      }

      /**
       * Constructs a cellwise vector space with compile-time component count.
       * @param mesh Mesh on which the object is defined.
       */
      template <size_t VDim>
      explicit P0(std::integral_constant<size_t, VDim>, const MeshType& mesh)
        : P0(mesh, VDim)
      {}

      /// Copies the space while retaining its mesh reference.
      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      P0(const P0& other) = default;
      /// Moves the space while retaining its mesh reference.
      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      P0(P0&& other) = default;
      ~P0() override = default;

      size_t getSize() const override
      {
        return m_mesh.get().getCellCount() * m_vdim;
      }
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }
      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /**
       * Returns the constant vector element for polytope @p i of dimension @p d.
       * @param i Index of the requested entry.
       * @returns The finite element.
       * @param d Topological dimension of the entity.
       */
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto geometry = getMesh().getGeometry(d, i);
        static thread_local ElementType element(Geometry::Polytope::Type::Segment, 1);
        if (element.getGeometry() != geometry || element.getCount() != m_vdim)
          element = ElementType(geometry, m_vdim);
        return element;
      }

      const IndexArray& getDOFs(size_t d, Index i) const override
      {
        assert(d == getMesh().getDimension());
        return m_dofs.at(i);
      }

      Index getGlobalIndex(
        const std::pair<size_t, Index>& idx, Index local) const override
      {
        assert(idx.first == getMesh().getDimension());
        assert(static_cast<size_t>(local) < m_vdim);
        return idx.second * m_vdim + local;
      }

      /**
       * Creates the physical-to-reference field pullback on a cell.
       * @param idx Index of the requested entry.
       * @param function Function to evaluate.
       * @returns The pullback.
       */
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& function) const
      {
        return Pullback<Callable>(*getMesh().getPolytope(idx.first, idx.second),
          std::forward<Callable>(function));
      }

      /**
       * Creates the reference-to-physical field pushforward on a cell.
       * @param function Function to evaluate.
       * @returns The pushforward.
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>&, Callable&& function) const
      {
        return Pushforward<Callable>(std::forward<Callable>(function));
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_vdim;
      std::vector<IndexArray> m_dofs;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for P0 from mesh - deduces to RealP0
   */
  template <class Context>
  P0(const Geometry::Mesh<Context>&) -> P0<Real, Geometry::Mesh<Context>>;

  /// Alias for a scalar real-valued P0 finite element space
  template <class Mesh>
  using RealP0 = P0<Real, Mesh>;

  /// Alias for a scalar complex-valued P0 finite element space
  template <class Mesh>
  using ComplexP0 = P0<Complex, Mesh>;
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Scalar, class Mesh>
  struct Traits<Variational::P0<Math::SpatialMatrix<Scalar>, Mesh>>
  {
      /// @brief Mesh type supplying topology and physical transformations.
      using MeshType = Mesh;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = Scalar;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<Scalar>;
      /// @brief Local or distributed execution context.
      using ContextType = typename Traits<Mesh>::ContextType;
      /// @brief Reference finite element type.
      using ElementType = Variational::P0Element<Math::SpatialMatrix<Scalar>>;
  };
}

namespace Rodin::Variational
{
  /// @brief Matrix-valued finite element space with component-wise scalar basis replication.
  template <class Scalar>
  class P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Local>,
        P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    public:
      /// @brief Mesh and execution context of this specialization.
      using MeshType = Geometry::Mesh<Context::Local>;
      /// @brief Scalar space supplying this family's DOF numbering.
      using ScalarSpace = P0<Scalar, Geometry::Mesh<Context::Local>>;
      /// @brief Matrix reference element for this family.
      using ElementType = P0Element<Math::SpatialMatrix<Scalar>>;
      /// @brief Existing finite-element-space interface.
      using Parent = FiniteElementSpace<MeshType,
        P0<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>>>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename ScalarSpace::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Local or distributed execution context.
      using ContextType = typename ScalarSpace::ContextType;
      /// @brief Expands scalar DOF maps into interleaved row-major matrix components.
      P0(const MeshType& mesh, size_t rows, size_t cols)
        : m_scalar(ScalarSpace(mesh)),
          m_rows(rows),
          m_cols(cols)
      {
        if (rows == 0 || cols == 0 || rows > RODIN_MAXIMAL_SPACE_DIMENSION ||
          cols > RODIN_MAXIMAL_SPACE_DIMENSION)
          Alert::Exception() << "SpatialMatrix ranges require 1 to 3 rows and columns."
                             << Alert::Raise;
        // Distributed scalar spaces may compute their size collectively.
        // Resolve it at construction, never from a local evaluation loop.
        m_size = m_scalar.getSize() * rows * cols;
        m_dofs.resize(mesh.getDimension() + 1);
        for (size_t d = mesh.getDimension(); d <= mesh.getDimension(); ++d)
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
      }

      /// @brief Copies the space and its DOF maps.
      P0(const P0&) = default;
      /// @brief Moves the space and its DOF maps.
      P0(P0&&) = default;
      /// @brief Copies the space and its DOF maps.
      P0& operator=(const P0&) = default;
      /// @brief Moves the space and its DOF maps.
      P0& operator=(P0&&) = default;

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

    private:
      ScalarSpace m_scalar;
      size_t m_rows, m_cols;
      size_t m_size = 0;
      std::vector<std::vector<IndexArray>> m_dofs;
      std::map<Geometry::Polytope::Type, ElementType> m_elements;
  };

  /// @brief Deduces a matrix range from explicit rows and columns.
  template <class Context>
  P0(const Geometry::Mesh<Context>&, size_t,
    size_t) -> P0<Math::SpatialMatrix<Real>, Geometry::Mesh<Context>>;

  /// @brief Matrix-valued discontinuous constant finite element space.
  template <class Mesh>
  using MatrixP0 = P0<Math::SpatialMatrix<Real>, Mesh>;
}

#endif
