/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_TRANSFORMATION_H
#define RODIN_GEOMETRY_TRANSFORMATION_H

/**
 * @file
 * @brief Base class for polytope geometric transformations.
 */

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>

#include "Rodin/Copyable.h"
#include "Rodin/Math.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Geometry/Polytope.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Abstract base class for polytope geometric transformations.
   *
   * This class represents the transformation function that maps reference
   * coordinates to physical coordinates for a polytope in a mesh.
   *
   * # Mathematical Foundation
   *
   * Let @f$ \tau @f$ denote a polytope in a triangulation @f$ \mathcal{T}_h @f$
   * with an associated reference element @f$ K @f$. This class represents the
   * transformation:
   * @f[
   *    x : K \subset \mathbb{R}^k \rightarrow \tau \subset \mathbb{R}^s
   * @f]
   * that maps a reference point @f$ r \in K @f$ to a physical point @f$ p \in \tau @f$:
   * @f[
   *    p = x(r)
   * @f]
   * where:
   * - @f$ k @f$ is the reference dimension (topological dimension of @f$ K @f$)
   * - @f$ s @f$ is the physical dimension (spatial dimension of @f$ \tau @f$)
   * - @f$ K @f$ is the reference element
   * - @f$ \tau @f$ is the physical element
   *
   * # Implementations
   *
   * Concrete implementations include:
   * - IdentityTransformation: @f$ x(r) = r @f$ (identity map)
   * - IsoparametricTransformation: Uses finite element basis functions
   *
   * # Thread Safety
   *
   * Transformation objects are immutable after construction and are thread-safe
   * for concurrent read access. The transform(), jacobian(), and inverse()
   * methods can be called concurrently from multiple threads.
   *
   * @see IdentityTransformation, IsoparametricTransformation, Point
   */
  class PolytopeTransformation : public Copyable
  {
    friend class boost::serialization::access;

    public:
      /**
       * @brief Constructs a transformation with given dimensions.
       * @param[in] rdim Reference dimension @f$ k @f$
       * @param[in] pdim Physical dimension @f$ s @f$
       */
      constexpr
      PolytopeTransformation(size_t rdim, size_t pdim)
        : m_rdim(rdim), m_pdim(pdim)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      PolytopeTransformation(const PolytopeTransformation&) = default;

      /**
       * @brief Move constructor.
       */
      constexpr
      PolytopeTransformation(PolytopeTransformation&&) = default;

      /**
       * @brief Move assignment operator.
       */
      PolytopeTransformation& operator=(PolytopeTransformation&&) = default;

      /**
       * @brief Virtual destructor.
       */
      virtual ~PolytopeTransformation() = default;

      /**
       * @brief Gets the reference dimension @f$ k @f$.
       * @returns Reference dimension (topological dimension of reference element)
       */
      constexpr
      size_t getReferenceDimension() const
      {
        return m_rdim;
      }

      /**
       * @brief Gets the physical dimension @f$ s @f$.
       * @returns Physical dimension (spatial dimension of physical element)
       */
      constexpr
      size_t getPhysicalDimension() const
      {
        return m_pdim;
      }

      /**
       * @brief Gets the polynomial order of the transformation.
       * @returns Polynomial order
       *
       * For linear transformations, returns 1. For higher-order geometric
       * approximations (e.g., quadratic elements), returns the appropriate order.
       */
      virtual size_t getOrder() const = 0;

      /**
       * @brief Gets the polynomial order of the Jacobian.
       * @returns Polynomial order of the Jacobian matrix entries
       *
       * For linear transformations, returns 0 (constant Jacobian). For higher-order
       * transformations, returns the order of the derivative.
       */
      virtual size_t getJacobianOrder() const = 0;

      /**
       * @brief Computes the physical coordinates from reference coordinates.
       * @param[out] pc Physical coordinates @f$ p @f$ (resized automatically)
       * @param[in] rc Reference coordinates @f$ r \in K @f$
       *
       * Computes:
       * @f[
       *    pc = x(rc)
       * @f]
       * The output vector is resized to the physical dimension @f$ s @f$.
       */
      virtual void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const = 0;

      /**
       * @brief Computes the Jacobian matrix of the transformation.
       * @param[out] jacobian Jacobian matrix @f$ \mathbf{J}_x @f$ (resized automatically)
       * @param[in] rc Reference coordinates @f$ r \in K @f$
       *
       * Computes the Jacobian matrix:
       * @f[
       *  \mathbf{J}_x (r) = \begin{bmatrix}
       * \dfrac{\partial x_1}{\partial r_1} & \ldots & \dfrac{\partial x_1}{\partial r_k}\\
       * \vdots & \ddots & \vdots\\
       * \dfrac{\partial x_s}{\partial r_1} & \ldots & \dfrac{\partial x_s}{\partial r_k}
       * \end{bmatrix}
       * @f]
       * for the transformation @f$ x : K \rightarrow \tau @f$.
       *
       * The output matrix is resized to @f$ s \times k @f$ where @f$ k @f$ is
       * the reference dimension and @f$ s @f$ is the physical dimension.
       */
      virtual void jacobian(Math::SpatialMatrix<Real>& jacobian, const Math::SpatialPoint& rc) const = 0;

      /**
       * @brief Computes the reference coordinates from physical coordinates.
       * @param[out] rc Reference coordinates @f$ r @f$ (resized automatically)
       * @param[in] pc Physical coordinates @f$ p \in \tau @f$
       *
       * Computes the inverse transformation:
       * @f[
       *    rc = x^{-1}(pc)
       * @f]
       *
       * @note Default implementation assumes the reference element has origin at 0.
       * Override for more complex inverse transformations.
       */
      virtual void inverse(Math::SpatialPoint& rc, const Math::SpatialPoint& pc) const;

      /**
       * @brief Serialization method for Boost.Serialization.
       * @param[in,out] ar Archive object
       */
      template<class Archive>
      void serialize(Archive & ar, const unsigned int)
      {
        ar & m_rdim;
        ar & m_pdim;
      }

      /**
       * @brief Creates a polymorphic copy of this transformation.
       * @returns Pointer to a new transformation object
       *
       * Derived classes must implement this to return a copy of their specific type.
       */
      virtual PolytopeTransformation* copy() const noexcept override = 0;

    private:
      size_t m_rdim; ///< Reference dimension @f$ k @f$
      size_t m_pdim; ///< Physical dimension @f$ s @f$
  };
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Rodin::Geometry::PolytopeTransformation)

#endif
