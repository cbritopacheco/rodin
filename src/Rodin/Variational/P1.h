/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1_H
#define RODIN_VARIATIONAL_P1_H

#include <boost/multi_array.hpp>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

#include "ForwardDecls.h"
#include "FiniteElement.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  /**
   * @brief First order Lagrange finite element
   *
   * @see @m_defelement{Lagrange,https://defelement.com/elements/lagrange.html}
   */
  class P1Element final : public FiniteElementBase<P1Element>
  {
    public:
      using Parent = FiniteElementBase<P1Element>;
      using G = Geometry::Polytope::Geometry;

      constexpr
      P1Element(Geometry::Polytope::Geometry geometry)
        : Parent(geometry)
      {}

      constexpr
      P1Element(const P1Element& other)
        : Parent(other)
      {}

      constexpr
      P1Element(P1Element&& other)
        : Parent(std::move(other))
      {}

      inline
      size_t getCount() const
      {
        return 3;
      }

      inline
      const Math::Matrix& getDOFs() const
      {
        const size_t g = static_cast<size_t>(getGeometry());
        assert(g > 0);
        assert(g < s_dofs.size());
        return s_dofs[g];
      }

      inline
      Math::Vector getBasis(const Math::Vector& r) const
      {
        const auto g = getGeometry();
        assert(static_cast<size_t>(r.size()) == Geometry::Polytope::getGeometryDimension(g));
        switch (g)
        {
          case G::Point:
            return Math::Vector{{1}};
          case G::Segment:
            return Math::Vector{{1 - r.x(), r.x()}};
          case G::Triangle:
            return Math::Vector{{-r.x() - r.y() + 1, r.x(), r.y()}};
          case G::Quadrilateral:
            return Math::Vector{{r.x() * r.y() - r.x() - r.y() + 1, r.x() * (1 - r.y()), r.y() * (1 - r.x()), r.x() * r.y()}};
          case G::Tetrahedron:
            return Math::Vector{{-r.x() - r.y() - r.z() + 1, r.x(), r.y(), r.z()}};
        }
      }

      inline
      Math::Matrix getGradient(const Math::Vector& r) const
      {
        const auto g = getGeometry();
        const size_t dim = Geometry::Polytope::getGeometryDimension(g);
        assert(static_cast<size_t>(r.size()) == dim);
        switch (g)
        {
          case G::Point:
            return Math::Matrix{{0}};
          case G::Segment:
            return Math::Matrix{{-1}, {1}};
          case G::Triangle:
            return Math::Matrix{{-1, -1}, {1, 0}, {0, 1}};
          case G::Quadrilateral:
            return Math::Matrix{{r.y() - 1, r.x() - 1}, {1 - r.y(), -r.x()}, {1 - r.x(), -r.y()}, {r.y(), r.x()}};
          case G::Tetrahedron:
            return Math::Matrix{{-1, -1, -1}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        }
      }

    private:
      static const std::array<Math::Matrix, 3> s_dofs;
  };

  template <>
  class P1<Scalar, Context::Serial, Geometry::Mesh<Context::Serial>> final
    : public FiniteElementSpaceBase
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using RangeType = Scalar;

      P1(const Geometry::Mesh<Context::Serial>& mesh);

      P1(const P1& other)
        : FiniteElementSpaceBase(other),
          m_mesh(other.m_mesh),
          m_elements(other.m_elements)
      {}

      P1(P1&& other)
        : FiniteElementSpaceBase(std::move(other)),
          m_mesh(other.m_mesh),
          m_elements(other.m_elements)
      {}

      P1& operator=(P1&& other) = default;

      inline
      const P1Element& getFiniteElement(size_t d, Index i) const
      {
        return m_elements[d][i];
      }

      inline
      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount();
      }

      inline
      size_t getVectorDimension() const override
      {
        return 1;
      }

      inline
      const Geometry::Mesh<Context::Serial>& getMesh() const override
      {
        return m_mesh.get();
      }

      inline
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        assert(i < static_cast<size_t>(p.size()));
        return p(local);
      }

    private:
      std::reference_wrapper<const Geometry::Mesh<Context::Serial>> m_mesh;
      std::vector<std::vector<P1Element>> m_elements;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&) -> P1<Scalar, Context, Geometry::Mesh<Context>>;

  template <>
  class P1<Math::Vector, Context::Serial, Geometry::Mesh<Context::Serial>> final
    : public FiniteElementSpaceBase
  {
    using KeyLeft = std::tuple<size_t, Index, Index>;
    using KeyRight = Index;
    using IndexMap = FlatMap<Index, Index>;

    public:
      using RangeType = Math::Vector;

      P1(const Geometry::Mesh<Context::Serial>& mesh, size_t vdim);

      P1(const P1& other) = default;

      P1(P1&& other) = default;

      P1& operator=(P1&& other) = default;

      inline
      const P1Element& getFiniteElement(size_t d, Index i) const
      {
        return m_elements[d][i];
      }

      inline
      size_t getSize() const override
      {
        return m_mesh.get().getVertexCount();
      }

      inline
      size_t getVectorDimension() const override
      {
        return m_vdim;
      }

      inline
      const Geometry::Mesh<Context::Serial>& getMesh() const override
      {
        return m_mesh.get();
      }

      inline
      Index getGlobalIndex(const std::pair<size_t, Index>& idx, Index local) const override
      {
        const auto [d, i] = idx;
        const auto& p = getMesh().getConnectivity().getPolytope(d, i);
        assert(i < static_cast<size_t>(p.size()));
        return p(local);
      }

    private:
      std::reference_wrapper<const Geometry::Mesh<Context::Serial>> m_mesh;
      size_t m_vdim;
      std::vector<std::vector<P1Element>> m_elements;
  };

  template <class Context>
  P1(const Geometry::Mesh<Context>&, size_t) -> P1<Math::Vector, Context, Geometry::Mesh<Context>>;
}

#endif
