/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P1ELEMENT_H
#define RODIN_VARIATIONAL_P1ELEMENT_H

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"
#include "Rodin/Geometry/GeometryIndexed.h"

#include "ForwardDecls.h"
#include "FiniteElement.h"
#include "FiniteElementSpace.h"

namespace Rodin::Variational
{
  /**
   * @brief Degree 1 scalar Lagrange element
   *
   * @see @m_defelement{Lagrange,https://defelement.com/elements/lagrange.html}
   */
  template <>
  class P1Element<Scalar> final : public FiniteElementBase<P1Element<Scalar>>
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
      constexpr
      FiniteElementMapping getMapping() const
      {
        return FiniteElementMapping::Identity;
      }

      inline
      size_t getCount() const
      {
        return Geometry::Polytope::getVertexCount(getGeometry());
      }

      inline
      const Math::Matrix& getDOFs() const
      {
        return s_dofs[getGeometry()];
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
        assert(false);
        return {};
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
        assert(false);
        return {};
      }

    private:
      static const Geometry::GeometryIndexed<Math::Matrix> s_dofs;
  };

  /**
   * @brief Degree 1 vector Lagrange element
   *
   * @see @m_defelement{Vector Lagrange,https://defelement.com/elements/vector-lagrange.html}
   */
  template <>
  class P1Element<Math::Vector> final : public FiniteElementBase<P1Element<Math::Vector>>
  {
    using G = Geometry::Polytope::Geometry;

    public:
      using Parent = FiniteElementBase<P1Element>;

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
      constexpr
      FiniteElementMapping getMapping() const
      {
        return FiniteElementMapping::Identity;
      }

      inline
      constexpr
      size_t getCount() const
      {
        return Geometry::Polytope::getVertexCount(getGeometry()
            ) * Geometry::Polytope::getGeometryDimension(getGeometry());
      }

      inline
      constexpr
      const Math::Matrix& getDOFs() const
      {
        return s_dofs[getGeometry()];
      }

      inline
      Math::Matrix getBasis(const Math::Vector& r) const
      {
        const auto g = getGeometry();
        switch (g)
        {
          case G::Point:
          {
            return Math::Matrix{{1}};
          }
          case G::Segment:
          {
            return Math::Matrix{{1 - r.x(), r.x()}};
          }
          case G::Triangle:
          {
            return Math::Matrix{
              {-r.x() - r.y() + 1, 0, r.x(), 0, r.y(), 0},
              {0, -r.x() - r.y() + 1, 0, r.x(), 0, r.y()}};
          }
          case G::Quadrilateral:
          {
            return Math::Matrix{
              {r.x() * r.y() - r.x() - r.y() + 1, 0, r.x() * (1 - r.y()), 0, r.y() * (1 - r.x()), 0, r.x() * r.y(), 0},
              {0, r.x() * r.y() - r.x() - r.y() + 1, 0, r.x() * (1 - r.y()), 0, r.y() * (1 - r.x()), 0, r.x() * r.y()}};
          }
          case G::Tetrahedron:
          {
            return Math::Matrix{
              {-r.x() - r.y() - r.z() + 1, 0, 0, r.x(), 0, 0, r.y(), 0, 0, r.z(), 0, 0},
              {0, -r.x() - r.y() - r.z() + 1, 0, 0, r.x(), 0, 0, r.y(), 0, 0, r.z(), 0},
              {0, 0, -r.x() - r.y() - r.z() + 1, 0, 0, r.x(), 0, 0, r.y(), 0, 0, r.z()}
            };
          }
        }
        assert(false);
        return {};
      }

      inline
      Math::Tensor<3> getJacobian(const Math::Vector& r) const
      {
        assert(false);
        return {};
        // const auto g = getGeometry();
        // switch (g)
        // {
        //   case G::Point:
        //   {
        //     return Math::Tensor<3>{{{0}}};
        //   }
        //   case G::Segment:
        //   {
        //     Math::Tensor<3> res(1, 1, getCount());
        //     res.setValues({{{-1}}, {{1}}});
        //     return res;
        //   }
        //   case G::Triangle:
        //   {
        //     Math::Tensor<3> res(2, 2, getCount());
        //     res.setValues({{{-1}}, {{1}}});
        //     return res;
        //   }
        // }
        // assert(false);
        // return {};
      }

    private:
      static const Geometry::GeometryIndexed<Math::Matrix> s_dofs;
  };

  /**
   * @brief Alias for P1Element<Scalar>
   */
  using ScalarP1Element = P1Element<Scalar>;

  /**
   * @brief Alias for P1Element<Math::Vector>
   */
  using VectorP1Element = P1Element<Math::Vector>;
}

#endif
