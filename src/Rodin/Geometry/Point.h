/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POINT_H
#define RODIN_GEOMETRY_POINT_H

/**
 * @file
 * @brief Spatial points on mesh polytopes with reference and physical coordinates.
 */

#include <set>
#include <iostream>
#include <array>
#include <optional>

#include "Rodin/Configure.h"

#include "Rodin/Array.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "Polytope.h"
#include "ForwardDecls.h"

#include "Types.h"

namespace Rodin::Geometry
{
  /**
   * @brief Base class for spatial points on a discrete mesh.
   *
   * This class represents a point in the computational domain, maintaining
   * both reference and physical coordinates along with the associated polytope.
   *
   * # Mathematical Foundation
   *
   * A Point represents the tuple @f$ (x, r, p, \tau) @f$ where:
   * - @f$ \tau @f$ is a polytope in the mesh @f$ \mathcal{T}_h @f$
   * - @f$ x : K \rightarrow \tau @f$ is the geometric transformation
   * - @f$ r \in K @f$ are the reference coordinates
   * - @f$ p = x(r) \in \tau @f$ are the physical coordinates
   *
   * The transformation satisfies:
   * @f[
   *  p = x(r)
   * @f]
   *
   * # Coordinate Systems
   *
   * Points maintain two coordinate systems:
   * - **Reference coordinates** @f$ r @f$: coordinates in the reference element @f$ K @f$
   * - **Physical coordinates** @f$ p @f$: coordinates in the physical element @f$ \tau @f$
   *
   * # Jacobian Information
   *
   * The class provides access to:
   * - Jacobian matrix @f$ \mathbf{J}_x(r) = \frac{\partial x}{\partial r} @f$
   * - Jacobian determinant @f$ |\mathbf{J}_x(r)| @f$
   * - Inverse Jacobian @f$ \mathbf{J}_x^{-1}(r) @f$
   * - Distortion measure
   *
   * # Thread Safety
   * This class is **not** thread-safe. Each thread should use its own Point instances.
   *
   * @see PolytopeTransformation, Polytope, Point
   */
  class PointBase
  {
    public:
      /**
       * @brief Enumeration of coordinate types.
       */
      enum class Coordinates
      {
        Reference, ///< Reference coordinates @f$ r \in K @f$
        Physical   ///< Physical coordinates @f$ p \in \tau @f$
      };

      /**
       * @brief Constructs a point on a polytope.
       * @param[in] polytope The polytope containing this point
       */
      explicit
      PointBase(const Polytope& polytope);

      /**
       * @brief Constructs a point with physical coordinates.
       * @param[in] polytope The polytope containing this point
       * @param[in] pc Physical coordinates
       */
      explicit
      PointBase(const Polytope& polytope, const Math::SpatialPoint& pc);

      /**
       * @brief Constructs a point with physical coordinates (move).
       * @param[in] polytope The polytope containing this point
       * @param[in] pc Physical coordinates
       */
      explicit
      PointBase(const Polytope& polytope, Math::SpatialPoint&& pc);

      /**
       * @brief Copy constructor.
       */
      PointBase(const PointBase& other);

      /**
       * @brief Move constructor.
       */
      PointBase(PointBase&& other);

      /**
       * @brief Gets the dimension of coordinates.
       * @param[in] coords Type of coordinates (Reference or Physical)
       * @returns Dimension of the specified coordinate system
       */
      size_t getDimension(Coordinates coords = Coordinates::Physical) const;

      /**
       * @brief Accesses the i-th physical coordinate.
       * @param[in] i Coordinate index
       * @returns i-th physical coordinate value
       */
      Real operator()(size_t i) const
      {
        return getPhysicalCoordinates()(i);
      }

      /**
       * @brief Gets physical coordinates as a vector.
       * @returns Reference to physical coordinate vector
       */
      const auto& vector() const
      {
        return getPhysicalCoordinates();
      }

      /**
       * @brief Gets the x-coordinate (first component).
       * @returns x-coordinate value
       */
      Real x() const
      {
        return operator()(0);
      }

      /**
       * @brief Gets the y-coordinate (second component).
       * @returns y-coordinate value
       */
      Real y() const
      {
        return operator()(1);
      }

      /**
       * @brief Gets the z-coordinate (third component).
       * @returns z-coordinate value
       */
      Real z() const
      {
        return operator()(2);
      }

      /**
       * @brief Computes the Euclidean norm.
       * @returns Euclidean norm @f$ \|p\|_2 @f$
       */
      Real norm() const
      {
        return vector().norm();
      }

      /**
       * @brief Computes the stable norm (avoids overflow).
       * @returns Stable norm
       */
      Real stableNorm() const
      {
        return vector().stableNorm();
      }

      /**
       * @brief Computes the Blue norm.
       * @returns Blue norm
       */
      Real blueNorm() const
      {
        return vector().blueNorm();
      }

      /**
       * @brief Computes the Lp norm.
       * @tparam p Norm order
       * @returns Lp norm @f$ \|p\|_p @f$
       */
      template <size_t p>
      Real lpNorm() const
      {
        return vector().lpNorm<p>();
      }

      /**
       * @brief Computes the squared Euclidean norm.
       * @returns Squared norm @f$ \|p\|_2^2 @f$
       */
      Real squaredNorm() const
      {
        return vector().squaredNorm();
      }

      /**
       * @brief Lexicographical comparison operator.
       * @param[in] p Point to compare with
       * @returns True if this point is lexicographically less than p
       */
      bool operator<(const PointBase& p) const;

      /**
       * @brief Gets the associated polytope.
       * @returns Reference to the polytope containing this point
       */
      const Polytope& getPolytope() const;

      /**
       * @brief Gets the physical coordinates.
       * @returns Reference to physical coordinate vector
       */
      const Math::SpatialVector<Real>& getPhysicalCoordinates() const;

      /**
       * @brief Gets coordinates of specified type.
       * @param[in] coords Type of coordinates to retrieve
       * @returns Reference to coordinate vector
       */
      const Math::SpatialVector<Real>& getCoordinates(Coordinates coords = Coordinates::Physical) const;

      /**
       * @brief Computes the Jacobian matrix at this point.
       * @returns Jacobian matrix @f$ \mathbf{J}_x(r) @f$
       *
       * The Jacobian is computed and cached on first access.
       */
      virtual const Math::SpatialMatrix<Real>& getJacobian() const;

      /**
       * @brief Gets the Jacobian determinant at this point.
       * @returns Determinant @f$ |\mathbf{J}_x(r)| @f$
       *
       * The determinant is cached on first access.
       */
      Real getJacobianDeterminant() const;

      /**
       * @brief Computes the inverse Jacobian matrix.
       * @returns Inverse Jacobian @f$ \mathbf{J}_x^{-1}(r) @f$
       *
       * The inverse is computed and cached on first access.
       */
      const Math::SpatialMatrix<Real>& getJacobianInverse() const;

      /**
       * @brief Computes the distortion measure at this point.
       * @returns Distortion value
       *
       * Measures how much the transformation distorts space at this point.
       */
      Real getDistortion() const;

      /**
       * @brief Sets the polytope for this point.
       * @param[in] polytope New polytope
       * @returns Reference to this point for method chaining
       */
      PointBase& setPolytope(const Polytope& polytope);

      /**
       * @brief Gets the reference coordinates (pure virtual).
       * @returns Reference coordinate vector
       *
       * Must be implemented by derived classes.
       */
      virtual const Math::SpatialVector<Real>& getReferenceCoordinates() const = 0;

    private:
      Polytope m_polytope;

      mutable Optional<Math::SpatialVector<Real>> m_pc;
      mutable Optional<Math::SpatialMatrix<Real>> m_jacobian;
      mutable Optional<Math::SpatialMatrix<Real>> m_jacobianInverse;
      mutable Optional<Real>                      m_jacobianDeterminant;
      mutable Optional<Real>                      m_distortion;
  };

  /**
   * @brief Concrete implementation of a spatial point on a mesh.
   *
   * This class provides a complete implementation of PointBase with storage
   * for reference coordinates. Points can be constructed with various
   * combinations of reference and physical coordinates.
   *
   * # Usage Examples
   *
   * @code{.cpp}
   * // Create a point at the centroid of a cell in reference coordinates
   * Math::SpatialPoint rc(3);
   * rc << 0.25, 0.25, 0.25;  // Barycentric coords for tetrahedron
   * Point p(cell, rc);
   *
   * // Access physical coordinates
   * std::cout << "x = " << p.x() << ", y = " << p.y() << ", z = " << p.z() << std::endl;
   *
   * // Compute Jacobian information
   * const auto& J = p.getJacobian();
   * Real detJ = p.getJacobianDeterminant();
   * @endcode
   *
   * # Thread Safety
   * Not thread-safe. Each thread should use separate Point instances.
   *
   * @see PointBase, Polytope, PolytopeTransformation
   */
  class Point final : public PointBase
  {
    public:
      /**
       * @brief Parent class type.
       */
      using Parent = PointBase;

      /**
       * @brief Constructs a point with reference coordinates.
       * @param[in] polytope Polytope containing the point
       * @param[in] rc Reference coordinates
       */
      explicit
      Point(const Polytope& polytope, const Math::SpatialPoint& rc);

      /**
       * @brief Constructs a point with reference coordinates (move).
       * @param[in] polytope Polytope containing the point
       * @param[in] rc Reference coordinates
       */
      explicit
      Point(const Polytope& polytope, Math::SpatialPoint&& rc);

      /**
       * @brief Constructs with both reference and physical coordinates.
       * @param[in] polytope Polytope containing the point
       * @param[in] rc Reference coordinates
       * @param[in] pc Physical coordinates
       */
      explicit
      Point(const Polytope& polytope, const Math::SpatialPoint& rc, const Math::SpatialPoint& pc);

      /**
       * @brief Constructs with reference and physical coordinates (mixed).
       * @param[in] polytope Polytope containing the point
       * @param[in] rc Reference coordinates
       * @param[in] pc Physical coordinates (move)
       */
      explicit
      Point(const Polytope& polytope, const Math::SpatialPoint& rc, Math::SpatialPoint&& pc);

      /**
       * @brief Constructs with reference and physical coordinates (mixed).
       * @param[in] polytope Polytope containing the point
       * @param[in] rc Reference coordinates (move)
       * @param[in] pc Physical coordinates
       */
      explicit
      Point(const Polytope& polytope, Math::SpatialPoint&& rc, const Math::SpatialPoint& pc);

      /**
       * @brief Constructs with reference and physical coordinates (move both).
       * @param[in] polytope Polytope containing the point
       * @param[in] rc Reference coordinates (move)
       * @param[in] pc Physical coordinates (move)
       */
      explicit
      Point(const Polytope& polytope, Math::SpatialPoint&& rc, Math::SpatialPoint&& pc);

      /**
       * @brief Copy constructor.
       */
      Point(const Point& other);

      /**
       * @brief Move constructor.
       */
      Point(Point&& other);

      /**
       * @brief Gets the reference coordinates.
       * @returns Reference coordinate vector
       */
      const Math::SpatialPoint& getReferenceCoordinates() const override;

    private:
      std::variant<const Math::SpatialPoint, std::reference_wrapper<Math::SpatialPoint>> m_rc;
  };

  /**
   * @brief Addition operator for Point and Eigen vector.
   * @tparam EigenDerived Eigen vector type
   * @param[in] p Point
   * @param[in] q Vector
   * @returns Result vector
   */
  template <class EigenDerived>
  auto
  operator+(const Geometry::Point& p, const Eigen::MatrixBase<EigenDerived>& q)
  {
    return p.vector() + q;
  }

  /**
   * @brief Addition operator for Eigen vector and Point.
   * @tparam EigenDerived Eigen vector type
   * @param[in] p Vector
   * @param[in] q Point
   * @returns Result vector
   */
  template <class EigenDerived>
  auto
  operator+(const Eigen::MatrixBase<EigenDerived>& p, const Geometry::Point& q)
  {
    return p + q.vector();
  }

  /**
   * @brief Subtraction operator for Point and Eigen vector.
   * @tparam EigenDerived Eigen vector type
   * @param[in] p Point
   * @param[in] q Vector
   * @returns Result vector
   */
  template <class EigenDerived>
  auto
  operator-(const Geometry::Point& p, const Eigen::MatrixBase<EigenDerived>& q)
  {
    return p.vector() - q;
  }

  /**
   * @brief Subtraction operator for Eigen vector and Point.
   * @tparam EigenDerived Eigen vector type
   * @param[in] p Vector
   * @param[in] q Point
   * @returns Result vector
   */
  template <class EigenDerived>
  auto
  operator-(const Eigen::MatrixBase<EigenDerived>& p, const Geometry::Point& q)
  {
    return p - q.vector();
  }

  /**
   * @brief Addition operator for two Points.
   * @param[in] p First point
   * @param[in] q Second point
   * @returns Result vector (sum of physical coordinates)
   */
  inline
  auto
  operator+(const Geometry::Point& p, const Geometry::Point& q)
  {
    return p.vector() + q.vector();
  }

  /**
   * @brief Subtraction operator for two Points.
   * @param[in] p First point
   * @param[in] q Second point
   * @returns Result vector (difference of physical coordinates)
   */
  inline
  auto
  operator-(const Geometry::Point& p, const Geometry::Point& q)
  {
    return p.vector() - q.vector();
  }

  /**
   * @brief Scalar multiplication (scalar * Point).
   * @param[in] s Scalar value
   * @param[in] p Point
   * @returns Scaled vector
   */
  inline
  auto
  operator*(Real s, const Geometry::Point& p)
  {
    return s * p.vector();
  }

  /**
   * @brief Scalar multiplication (Point * scalar).
   * @param[in] p Point
   * @param[in] s Scalar value
   * @returns Scaled vector
   */
  inline
  auto
  operator*(const Geometry::Point& p, Real s)
  {
    return p.vector() * s;
  }

  template <class Scalar>
  auto operator+(
    const Math::SpatialVector<Scalar>& v, const Geometry::Point& p)
  {
    return v + p.vector();
  }

  template <class Scalar>
  auto operator+(
    const Geometry::Point& p, const Math::SpatialVector<Scalar>& v)
  {
    return p.vector() + v;
  }

  template <class Scalar>
  auto operator-(
    const Math::SpatialVector<Scalar>& v, const Geometry::Point& p)
  {
    return v - p.vector();
  }

  template <class Scalar>
  auto operator-(
    const Geometry::Point& p, const Math::SpatialVector<Scalar>& v)
  {
    return p.vector() - v;
  }
}
#endif
