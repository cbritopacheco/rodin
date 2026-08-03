/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_P0G_H
#define RODIN_VARIATIONAL_P0G_P0G_H

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

/// @cond RODIN_DOXYGEN_INTERNAL
namespace Rodin::FormLanguage
{
  template <class Number, class Mesh>
  struct Traits<Variational::P0g<Number, Mesh>>
  {
    using MeshType    = Mesh;
    using ScalarType  = Number;
    using RangeType   = ScalarType;
    /// @brief Execution context type.
    using ContextType = typename FormLanguage::Traits<MeshType>::ContextType;
    /// @brief Finite element type.
    using ElementType = Variational::P0gElement<RangeType>;
  };

  template <class Number, class Mesh>
  struct Traits<Variational::P0g<Math::SpatialVector<Number>, Mesh>>
  {
    using MeshType    = Mesh;
    using ScalarType  = Number;
    using RangeType   = Math::SpatialVector<ScalarType>;
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
  template <class Range, class Mesh>
  class P0g;

  // --------------------------------------------------------------------------
  // Scalar P0g<Real, Mesh<Local>>
  // --------------------------------------------------------------------------
  template <>
  class P0g<Real, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>,
        P0g<Real, Geometry::Mesh<Context::Local>>>
  {
    public:
      using ScalarType  = Real;
      using RangeType   = ScalarType;
      /// @brief Execution context type.
      using ContextType = Context::Local;
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<RangeType>;
      using Parent      = FiniteElementSpace<MeshType, P0g<RangeType, MeshType>>;

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

          auto operator()(const Math::SpatialPoint& r) const
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
          explicit Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      explicit P0g(const MeshType& mesh)
        : m_mesh(mesh)
      {}

      P0g(const P0g& other)
        : Parent(other),
          m_mesh(other.m_mesh)
      {}

      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_mesh(std::move(other.m_mesh))
      {}

      ~P0g() override = default;

      P0g& operator=(const P0g& other)
      {
        if (this != &other)
        {
          Parent::operator=(other);
          m_mesh = other.m_mesh;
        }
        return *this;
      }

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
      std::reference_wrapper<const MeshType> m_mesh;
  };

  // --------------------------------------------------------------------------
  // Vector P0g<Math::SpatialVector<Real>, Mesh<Local>>
  // --------------------------------------------------------------------------
  template <>
  class P0g<Math::SpatialVector<Real>, Geometry::Mesh<Context::Local>> final
    : public FiniteElementSpace<
        Geometry::Mesh<Context::Local>,
        P0g<Math::SpatialVector<Real>, Geometry::Mesh<Context::Local>>>
  {
    public:
      using ScalarType  = Real;
      using RangeType   = Math::SpatialVector<Real>;
      /// @brief Execution context type.
      using ContextType = Context::Local;
      using MeshType    = Geometry::Mesh<ContextType>;
      /// @brief Finite element type.
      using ElementType = P0gElement<Math::SpatialVector<ScalarType>>;
      using Parent      = FiniteElementSpace<MeshType, P0g<Math::SpatialVector<Real>, MeshType>>;

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

          auto operator()(const Math::SpatialPoint& r) const
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
          explicit Pushforward(Function&& v)
            : m_v(std::forward<Function>(v))
          {}

          constexpr
          auto operator()(const Geometry::Point& p) const
          {
            return m_v(p.getReferenceCoordinates());
          }

        private:
          CallableType m_v;
      };

      explicit P0g(const MeshType& mesh, size_t vdim)
        : m_mesh(mesh), m_vdim(vdim)
      {
        assert(m_vdim > 0);
        m_dofs.resize(m_vdim);
        for (size_t k = 0; k < m_vdim; ++k)
          m_dofs[k] = static_cast<Index>(k);
      }

      template <size_t VDim>
      explicit P0g(std::integral_constant<size_t, VDim>, const MeshType& mesh)
        : P0g(mesh, VDim)
      {}

      P0g(const P0g& other)
        : Parent(other),
          m_dofs(other.m_dofs),
          m_mesh(other.m_mesh),
          m_vdim(other.m_vdim)
      {}

      P0g(P0g&& other)
        : Parent(std::move(other)),
          m_dofs(std::move(other.m_dofs)),
          m_mesh(std::move(other.m_mesh)),
          m_vdim(std::move(other.m_vdim))
      {}

      ~P0g() override = default;

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
      IndexArray m_dofs;
      std::reference_wrapper<const MeshType> m_mesh;
      size_t m_vdim;
  };

  // CTAD (scalar)
  template <class Context>
  P0g(const Geometry::Mesh<Context>&) -> P0g<Real, Geometry::Mesh<Context>>;

  // Aliases (scalar spaces)
  template <class Mesh>
  using RealP0g = P0g<Real, Mesh>;

  template <class Mesh>
  using ComplexP0g = P0g<Complex, Mesh>;

  // Aliases (vector spaces)
  template <class Mesh>
  using VectorP0g = P0g<Math::SpatialVector<Real>, Mesh>;
}

/// @endcond
#endif
