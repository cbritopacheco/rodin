/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0_P0_H
#define RODIN_VARIATIONAL_P0_P0_H

#include <boost/multi_array.hpp>

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"
#include "P0Element.h"

namespace Rodin::FormLanguage
{
  template <class Number, class Mesh>
  struct Traits<Variational::P0<Number, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Number;
    using RangeType = ScalarType;
    using ContextType = typename MeshType::Context;
    using ElementType = Variational::P0Element<RangeType>;
  };

  template <class Number, class Mesh>
  struct Traits<Variational::P0<Math::Vector<Number>, Mesh>>
  {
    using MeshType = Mesh;
    using ScalarType = Number;
    using RangeType = Math::Vector<ScalarType>;
    using ContextType = typename MeshType::Context;
    using ElementType = Variational::P0Element<RangeType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup P0Specializations P0 Template Specializations
   * @brief Template specializations of the P0 class.
   * @see P0
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
  template <>
  class P0<Real, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>, P0<Real, Geometry::Mesh<Context::Local>>>
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using ScalarType = Real;

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
          m_mesh(other.m_mesh)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other P0 space to move from
       */
      P0(P0&& other)
        : Parent(std::move(other)),
          m_mesh(other.m_mesh)
      {}

      virtual ~P0() override = default;

      /**
       * @brief Move assignment operator.
       * @param[in] other P0 space to move from
       * @return Reference to this P0 space
       */
      P0& operator=(P0&& other) = default;

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
        static thread_local constexpr ElementType s_null(Geometry::Polytope::Type::Point);
        return s_null;
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
       * @param[in] idx Pair of (dimension, polytope index)
       * @param[in] v Function to push forward
       * @return Pushforward transformation object
       *
       * The pushforward maps a function from reference space to physical space.
       */
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      std::vector<IndexArray> m_dofs;
      std::reference_wrapper<const MeshType> m_mesh;
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

#endif

