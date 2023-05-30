#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Connectivity.h"

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
      auto getBasis(size_t local) const
      {
        const size_t count = getCount();
        return Math::Matrix::Identity(count, count).row(local);
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

    private:
  };

  using ScalarP1Element = P1Element<Scalar>;

  using VectorP1Element = P1Element<Math::Vector>;
}
