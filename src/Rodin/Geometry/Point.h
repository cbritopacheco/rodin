/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POINT_H
#define RODIN_GEOMETRY_POINT_H

#include <set>
#include <iostream>
#include <array>
#include <optional>

#include "Rodin/Configure.h"

#include "Rodin/Array.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Threads/Mutable.h"

#include "Polytope.h"
#include "ForwardDecls.h"

#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class for spatial points on a discrete mesh.
   *
   * This class represents the tuple @f$ (x, r, p) @f$
   * such that:
   * @f[
   *  p = x(r)
   * @f]
   * for a polytope @f$ \tau \in \mathcal{T}_h @f$ belonging to the mesh @f$
   * \mathcal{T}_h @f$. Here, @f$ p \in \tau @f$ denotes the physical
   * coordinates of the point, while @f$ x : K \rightarrow \tau @f$ represents
   * the transformation taking reference coordinates @f$ r \in K @f$, for a
   * reference geometry @f$ K @f$.
   *
   * @section rodin-geometry-point-thread_safety Thread safety
   * This class is not thread safe.
   *
   * @see PolytopeTransformation
   */
  class PointBase
  {
    enum class PolytopeStorage
    {
      Value,
      Reference
    };

    public:
      /// Denotes the type of coordinates.
      enum class Coordinates
      {
        Reference, ///< Reference coordinates
        Physical ///< Physical coordinates
      };

      explicit
      PointBase(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          const Math::SpatialVector<Scalar>& pc);

      explicit
      PointBase(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const PolytopeTransformation> trans);

      explicit
      PointBase(
          Polytope&& polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          const Math::SpatialVector<Scalar>& pc);

      explicit
      PointBase(
          Polytope&& polytope,
          std::reference_wrapper<const PolytopeTransformation> trans);

      /**
       * @brief Copy constructor.
       */
      PointBase(const PointBase& other);

      /**
       * @brief Move constructor.
       */
      PointBase(PointBase&& other);

      /**
       * @brief Gets the space dimension of the physical coordinates.
       * @returns Dimension of the physical coordinates.
       */
      size_t getDimension(Coordinates coords = Coordinates::Physical) const;

      /**
       * @brief Gets the i-th physical coordinate.
       * @returns Physical i-th coordinate.
       */
      inline
      Scalar operator()(size_t i, Coordinates coords = Coordinates::Physical) const
      {
        switch (coords)
        {
          case Coordinates::Physical:
          {
            return getPhysicalCoordinates()(i);
          }
          case Coordinates::Reference:
          {
            return getReferenceCoordinates()(i);
          }
        }
        assert(false);
        return NAN;
      }

      inline
      const auto& asVector() const
      {
        return getPhysicalCoordinates();
      }

      /**
       * @brief Gets the @f$ x @f$ coordinate.
       * @returns @f$ x @f$ coordinate of the point.
       */
      inline
      Scalar x(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(0, coords);
      }

      /**
       * @brief Gets the @f$ y @f$ coordinate.
       * @returns @f$ y @f$ coordinate of the point.
       */
      inline
      Scalar y(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(1, coords);
      }

      /**
       * @brief Gets the @f$ z @f$ coordinate.
       * @returns @f$ z @f$ coordinate of the point.
       */
      inline
      Scalar z(Coordinates coords = Coordinates::Physical) const
      {
        return operator()(2, coords);
      }

      inline
      Scalar norm() const
      {
        return asVector().norm();
      }

      inline
      Scalar stableNorm() const
      {
        return asVector().stableNorm();
      }

      inline
      Scalar blueNorm() const
      {
        return asVector().blueNorm();
      }

      template <size_t p>
      inline
      Scalar lpNorm() const
      {
        return asVector().lpNorm<p>();
      }

      inline
      Scalar squaredNorm() const
      {
        return asVector().squaredNorm();
      }

      /**
       * @brief Lexicographical comparison.
       */
      bool operator<(const PointBase& p) const;

      const Polytope& getPolytope() const;

      inline
      const PolytopeTransformation& getTransformation() const
      {
        return m_trans.get();
      }

      const Math::SpatialVector<Scalar>& getPhysicalCoordinates() const;

      const Math::SpatialVector<Scalar>& getCoordinates(Coordinates coords = Coordinates::Physical) const;

      /**
       * @brief Computes the Jacobian matrix of the transformation at the
       * point.
       */
      virtual const Math::SpatialMatrix<Scalar>& getJacobian() const;

      Scalar getJacobianDeterminant() const;

      /**
       * @brief Computes the inverse of the Jacobian matrix of the
       * transformation at the point.
       */
      const Math::SpatialMatrix<Scalar>& getJacobianInverse() const;

      /**
       * @brief Computes the distortion of space of the transformation at the
       * point.
       */
      Scalar getDistortion() const;

      virtual const Math::SpatialVector<Scalar>& getReferenceCoordinates() const = 0;

    private:
      PolytopeStorage m_polytopeStorage;
      std::variant<const Polytope, std::reference_wrapper<const Polytope>> m_polytope;
      std::reference_wrapper<const PolytopeTransformation> m_trans;

      mutable Threads::Mutable<std::optional<const Math::SpatialVector<Scalar>>> m_pc;
      mutable Threads::Mutable<std::optional<const Math::SpatialMatrix<Scalar>>> m_jacobian;
      mutable Threads::Mutable<std::optional<const Math::SpatialMatrix<Scalar>>> m_jacobianInverse;
      mutable Threads::Mutable<std::optional<const Scalar>>              m_jacobianDeterminant;
      mutable Threads::Mutable<std::optional<const Scalar>>              m_distortion;
  };

  /**
   * @brief Represents a spatial point on a discrete mesh.
   *
   * This class represents the tuple @f$ (x, r, p) @f$
   * such that:
   * @f[
   *  p = x(r)
   * @f]
   * for a polytope @f$ \tau \in \mathcal{T}_h @f$ belonging to the mesh @f$
   * \mathcal{T}_h @f$. Here, @f$ p \in \tau @f$ denotes the physical
   * coordinates of the point, while @f$ x : K \rightarrow \tau @f$ represents
   * the transformation taking reference coordinates @f$ r \in K @f$, for a
   * reference geometry @f$ K @f$.
   *
   * @section rodin-geometry-point-thread_safety Thread safety
   * This class is not thread safe.
   *
   * @see PolytopeTransformation
   */
  class Point final : public PointBase
  {
    enum class RCStorage
    {
      Value,
      Reference
    };

    public:
      using Parent = PointBase;

      explicit
      Point(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const Math::SpatialVector<Scalar>> rc,
          const Math::SpatialVector<Scalar>& pc)
        : Point(polytope, std::cref(polytope.get().getTransformation()), rc, pc)
      {}

      explicit
      Point(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          std::reference_wrapper<const Math::SpatialVector<Scalar>> rc,
          const Math::SpatialVector<Scalar>& pc);

      explicit
      Point(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          Math::SpatialVector<Scalar>&& rc,
          const Math::SpatialVector<Scalar>& pc);

      explicit
      Point(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          std::reference_wrapper<const Math::SpatialVector<Scalar>> rc);

      explicit
      Point(
          std::reference_wrapper<const Polytope> polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          Math::SpatialVector<Scalar>&& rc);

      explicit
      Point(
          Polytope&& polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          std::reference_wrapper<const Math::SpatialVector<Scalar>> rc,
          const Math::SpatialVector<Scalar>& pc);

      explicit
      Point(
          Polytope&& polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          Math::SpatialVector<Scalar>&& rc,
          const Math::SpatialVector<Scalar>& pc);

      explicit
      Point(
          Polytope&& polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          std::reference_wrapper<const Math::SpatialVector<Scalar>> rc);

      explicit
      Point(
          Polytope&& polytope,
          std::reference_wrapper<const PolytopeTransformation> trans,
          Math::SpatialVector<Scalar>&& rc);

      Point(const Point& other);

      Point(Point&& other);

      const Math::SpatialVector<Scalar>& getReferenceCoordinates() const override;

    private:
      const RCStorage m_rcStorage;
      std::variant<const Math::SpatialVector<Scalar>, std::reference_wrapper<const Math::SpatialVector<Scalar>>> m_rc;
  };

  inline
  std::ostream& operator<<(std::ostream& os, Polytope::Type g)
  {
    switch (g)
    {
      case Polytope::Type::Point:
      {
        os << "Point";
        break;
      }
      case Polytope::Type::Segment:
      {
        os << "Segment";
        break;
      }
      case Polytope::Type::Triangle:
      {
        os << "Triangle";
        break;
      }
      case Polytope::Type::Quadrilateral:
      {
        os << "Quadrilateral";
        break;
      }
      case Polytope::Type::Tetrahedron:
      {
        os << "Tetrahedron";
        break;
      }
    }
    return os;
  }

  inline
  auto
  operator+(const Geometry::Point& p, const Geometry::Point& q)
  {
    return p.asVector() + q.asVector();
  }

  inline
  auto
  operator-(const Geometry::Point& p, const Geometry::Point& q)
  {
    return p.asVector() - q.asVector();
  }

  template <class EigenDerived>
  inline
  auto
  operator+(const Geometry::Point& p, const Eigen::MatrixBase<EigenDerived>& q)
  {
    return p.asVector() + q;
  }

  template <class EigenDerived>
  inline
  auto
  operator+(const Eigen::MatrixBase<EigenDerived>& p, const Geometry::Point& q)
  {
    return p + q.asVector();
  }

  template <class EigenDerived>
  inline
  auto
  operator-(const Geometry::Point& p, const Eigen::MatrixBase<EigenDerived>& q)
  {
    return p.asVector() - q;
  }

  template <class EigenDerived>
  inline
  auto
  operator-(const Eigen::MatrixBase<EigenDerived>& p, const Geometry::Point& q)
  {
    return p - q.asVector();
  }

  inline
  auto
  operator*(Scalar s, const Geometry::Point& p)
  {
    return s * p.asVector();
  }

  inline
  auto
  operator*(const Geometry::Point& p, Scalar s)
  {
    return p.asVector() * s;
  }
}
#endif
