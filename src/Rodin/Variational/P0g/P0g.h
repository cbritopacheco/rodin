/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_P0G_H
#define RODIN_VARIATIONAL_P0G_P0G_H

#include <map>
#include <cassert>
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include "Rodin/Types.h"

#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"

#include "Rodin/Variational/FiniteElementSpace.h"

#include "ForwardDecls.h"

#include "P0gElement.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type traits for @c P0g: exposes the mesh type, the scalar type, the range
   * type, the execution context and the finite element type.
   */
  template <class Number, class Mesh>
  struct Traits<Variational::P0g<Number, Mesh>>
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
      using ElementType = Variational::P0gElement<RangeType>;
  };

  /**
   * @brief Type traits for @c P0g: exposes the mesh type, the scalar type, the range
   * type, the execution context and the finite element type.
   */
  template <class Number, class Mesh>
  struct Traits<Variational::P0g<Math::SpatialVector<Number>, Mesh>>
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
      using ElementType = Variational::P0gElement<Math::SpatialVector<ScalarType>>;
  };
}

namespace Rodin::Variational
{
  /**
   * @brief Global constant (P0g) finite element space.
   *
   * - Scalar-valued: 1 DOF total (span{1})
   * - Vector-valued: vdim DOFs total (span{e_0,...,e_{vdim-1}})
   *
   * This is the global analogue of P0 (which is elementwise constant).
   */
  /**
   * @brief Supported P0g range specializations.
   * | Range | Local context | MPI context |
   * |-------|---------------|-------------|
   * | Scalar | Supported | Supported |
   * | SpatialVector | Supported | Supported |
   * | SpatialMatrix | Explicit rows and columns | Component-wise ownership and ghosts |
   */
  template <class Range, class Mesh>
  class P0g;

  // --------------------------------------------------------------------------
  // Scalar P0g<Real, Mesh<Local>>
  // --------------------------------------------------------------------------
  /// @brief Cellwise-constant scalar finite element space with a global basis.
  template <class Scalar>
    requires(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>)
  class P0g<Scalar, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Local>,
        P0g<Scalar, Geometry::Mesh<Context::Local>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Range (evaluation value) type.
      using RangeType   = ScalarType;
      /// @brief Execution context type.
      using ContextType = Context::Local;
      /// @brief Mesh type.
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<RangeType>;
      /// @brief Parent class type.
      using Parent      = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

      /// @brief Pullback of a P0g function to the reference element.
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

      /// @brief Pushforward of a P0g function to the physical element.
      template <class Callable>
      class Pushforward :
        public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          /// @brief Callable type evaluated on physical points.
          using CallableType = Callable;

          /**
           * @brief Constructs the pushforward of a function.
           * @param v Function operand.
           */
          template <class Function>
          explicit Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

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

      /// @brief Constructs the P0g from the given arguments.
      explicit P0g(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      /// @brief Copy constructor.
      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      /// @brief Move constructor.
      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh))
      {}

      ~P0g() override = default;

      /// @brief Copy assignment.
      P0g& operator=(const P0g& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
        }
        return *this;
      }

      /// @brief Move assignment.
      P0g& operator=(P0g&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_mesh = std::move(other.m_mesh);
        }
        return *this;
      }

      size_t getSize() const override
      {
        return 1;
      }

      size_t getVectorDimension() const override
      {
        return 1;
      }

      const MeshType& getMesh() const override
      {
        return m_mesh.get();
      }

      /// @brief Gets the finite element attached to a polytope.
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          case Geometry::Polytope::Type::Segment:
          case Geometry::Polytope::Type::Triangle:
          case Geometry::Polytope::Type::Quadrilateral:
          case Geometry::Polytope::Type::Tetrahedron:
          case Geometry::Polytope::Type::Pyramid:
          case Geometry::Polytope::Type::Wedge:
          case Geometry::Polytope::Type::Hexahedron:
          {
            // Unlike the fixed P0/P1/H1 elements, this cached element is
            // reassigned when the requested geometry changes. It must remain
            // thread-local while the FES API returns it by reference.
            static thread_local ElementType e(Geometry::Polytope::Type::Point);
            if (e.getGeometry() != g)
              e = ElementType(g);
            return e;
          }
        }

        assert(false);
        static thread_local ElementType nullElem(Geometry::Polytope::Type::Point);
        return nullElem;
      }

      const IndexArray& getDOFs(size_t, Index) const override
      {
        static const IndexArray s_dofs{{0}};
        return s_dofs;
      }

      Index getGlobalIndex(const std::pair<size_t, Index>&, Index) const override
      {
        return 0;
      }

      /// @brief Gets the pullback of a callable on a polytope.
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /// @brief Gets the pushforward of a callable on a polytope.
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>&, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      std::reference_wrapper<const MeshType> m_mesh;
  };

  // --------------------------------------------------------------------------
  // Vector P0g<Math::SpatialVector<Real>, Mesh<Local>>
  // --------------------------------------------------------------------------
  /// @brief Cellwise-constant vector-valued finite element space with a global basis.
  template <class Scalar>
    requires(std::is_same_v<Scalar, Real> || std::is_same_v<Scalar, Complex>)
  class P0g<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Local>,
        P0g<Math::SpatialVector<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<Scalar>;
      /// @brief Execution context type.
      using ContextType = Context::Local;
      /// @brief Mesh type.
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<Math::SpatialVector<ScalarType>>;
      /// @brief Parent class type.
      using Parent = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

      /// @brief Pullback of a vector-valued P0g function to the reference element.
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
           * @brief Evaluates at a point on the reference element.
           * @param r Reference coordinates at which to evaluate the basis.
           * @returns Reference to the entry at the supplied indices.
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

      /// @brief Pushforward of a vector-valued P0g function to the physical element.
      template <class Callable>
      class Pushforward :
        public FiniteElementSpacePushforwardBase<Pushforward<Callable>>
      {
        public:
          /// @brief Callable type evaluated on physical points.
          using CallableType = Callable;

          /**
           * @brief Constructs the pushforward of a function.
           * @param v Function operand.
           */
          template <class Function>
          explicit Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          /**
           * @brief Evaluates at a geometric point.
           * @param p Point at which the operation is evaluated.
           * @returns Reference to the entry at the supplied indices.
           */
          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      /// @brief Constructs the P0g from the given arguments.
      explicit P0g(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh), m_vdim(vdim)
      {
        assert(m_vdim > 0);
        m_dofs.resize(m_vdim);
        for (size_t k = 0; k < m_vdim; ++k)
          m_dofs[k] = static_cast<Index>(k);
      }

      /// @brief Constructs the P0g from the given arguments.
      template <size_t VDim>
      explicit P0g(std::integral_constant<size_t, VDim>, const MeshType& mesh)
        : P0g(mesh, VDim)
      {}

      /// @brief Copy constructor.
      P0g(const P0g& other)
        : Parent(other),
          m_dofs(other.m_dofs),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim)
      {}

      /// @brief Move constructor.
      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_dofs(std::move(other.m_dofs)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(std::move(other.m_vdim))
      {}

      ~P0g() override = default;

      /// @brief Copy assignment.
      P0g& operator=(const P0g& other)
      {
        Parent::operator=(other);
        if (this != &other)
        {
          m_dofs = other.m_dofs;
          m_vdim = other.m_vdim;
          m_mesh = other.m_mesh;
        }
        return *this;
      }

      /// @brief Move assignment.
      P0g& operator=(P0g&& other)
      {
        Parent::operator=(std::move(other));
        if (this != &other)
        {
          m_dofs = std::move(other.m_dofs);
          m_vdim = std::move(other.m_vdim);
          m_mesh = std::move(other.m_mesh);
        }
        return *this;
      }

      size_t getSize() const override { return m_vdim; }

      size_t getVectorDimension() const override { return m_vdim; }

      const MeshType& getMesh() const override { return m_mesh.get(); }

      /// @brief Gets the finite element attached to a polytope.
      const ElementType& getFiniteElement(size_t d, Index i) const
      {
        const auto g = getMesh().getGeometry(d, i);
        switch (g)
        {
          case Geometry::Polytope::Type::Point:
          case Geometry::Polytope::Type::Segment:
          case Geometry::Polytope::Type::Triangle:
          case Geometry::Polytope::Type::Quadrilateral:
          case Geometry::Polytope::Type::Tetrahedron:
          case Geometry::Polytope::Type::Pyramid:
          case Geometry::Polytope::Type::Wedge:
          case Geometry::Polytope::Type::Hexahedron:
          {
            // Geometry and vector dimension both select the returned element;
            // sharing this mutable cache would race between evaluations.
            static thread_local ElementType e(Geometry::Polytope::Type::Point, 1);
            if (e.getGeometry() != g || e.getCount() != m_vdim)
              e = ElementType(g, m_vdim);
            return e;
          }
        }

        assert(false);
        static thread_local ElementType nullElem(Geometry::Polytope::Type::Point, 1);
        if (nullElem.getCount() != m_vdim)
          nullElem = ElementType(Geometry::Polytope::Type::Point, m_vdim);
        return nullElem;
      }

      const IndexArray& getDOFs(size_t, Index) const override
      {
        return m_dofs;
      }

      Index getGlobalIndex(const std::pair<size_t, Index>&, Index local) const override
      {
        assert(static_cast<size_t>(local) < m_vdim);
        return local;
      }

      /// @brief Gets the pullback of a callable on a polytope.
      template <class Callable>
      auto getPullback(const std::pair<size_t, Index>& idx, Callable&& v) const
      {
        const auto& [d, i] = idx;
        const auto& mesh = getMesh();
        return Pullback<Callable>(*mesh.getPolytope(d, i), std::forward<Callable>(v));
      }

      /// @brief Gets the pushforward of a callable on a polytope.
      template <class Callable>
      auto getPushforward(const std::pair<size_t, Index>&, Callable&& v) const
      {
        return Pushforward<Callable>(std::forward<Callable>(v));
      }

    private:
      IndexArray m_dofs;
      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_vdim;
  };

  // CTAD (scalar)
  /// @brief Deduction guide for @c P0g.
  template <class Context>
  P0g(const Geometry::Mesh<Context>&) -> P0g<Real, Geometry::Mesh<Context>>;

  // Aliases (scalar spaces)
  /// @brief Cellwise-constant real space with a global basis.
  template <class Mesh>
  using RealP0g = P0g<Real, Mesh>;

  /// @brief Cellwise-constant complex space with a global basis.
  template <class Mesh>
  using ComplexP0g = P0g<Complex, Mesh>;

  // Aliases (vector spaces)
  /// @brief Cellwise-constant vector-valued space with a global basis.
  template <class Mesh>
  using VectorP0g = P0g<Math::SpatialVector<Real>, Mesh>;
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Scalar, class Mesh>
  struct Traits<Variational::P0g<Math::SpatialMatrix<Scalar>, Mesh>>
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
      using ElementType = Variational::P0gElement<Math::SpatialMatrix<Scalar>>;
  };
}

namespace Rodin::Variational
{
  /// @brief Matrix-range finite element or expression specialization.
  template <class Scalar>
  class P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<Geometry::Mesh<Context::Local>,
        P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>>>
  {
    public:
      /// @brief Mesh and execution context of this specialization.
      using MeshType = Geometry::Mesh<Context::Local>;
      /// @brief Scalar space supplying this family's DOF numbering.
      using ScalarSpace = P0g<Scalar, Geometry::Mesh<Context::Local>>;
      /// @brief Matrix reference element for this family.
      using ElementType = P0gElement<Math::SpatialMatrix<Scalar>>;
      /// @brief Existing finite-element-space interface.
      using Parent = FiniteElementSpace<MeshType,
        P0g<Math::SpatialMatrix<Scalar>, Geometry::Mesh<Context::Local>>>;
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

    private:
      ScalarSpace m_scalar;
      size_t m_rows, m_cols;
      size_t m_size = 0;
      std::vector<std::vector<IndexArray>> m_dofs;
      std::map<Geometry::Polytope::Type, ElementType> m_elements;
  };

  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class Context>
  P0g(const Geometry::Mesh<Context>&, size_t,
    size_t) -> P0g<Math::SpatialMatrix<Real>, Geometry::Mesh<Context>>;

  /// @brief Matrix-valued globally constant finite element space.
  template <class Mesh>
  using MatrixP0g = P0g<Math::SpatialMatrix<Real>, Mesh>;
}

#endif
