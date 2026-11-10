/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_POLYTOPETRANSFORMATION_H
#define RODIN_GEOMETRY_POLYTOPETRANSFORMATION_H

#include <vector>
#include <initializer_list>
#include <type_traits>

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
#include "Rodin/QF/ForwardDecls.h"

#include "Bounds.h"
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
   * transformation.
   * @f[
   *    x : K \subset \mathbb{R}^k \rightarrow \tau \subset \mathbb{R}^s
   * @f]
   * that maps a reference point @f$ r \in K @f$ to a physical point @f$ p \in \tau @f$.
   * @f[
   *    p = x(r)
   * @f]
   * Here,
   * - @f$ k @f$ is the reference dimension (topological dimension of @f$ K @f$)
   * - @f$ s @f$ is the physical dimension (spatial dimension of @f$ \tau @f$)
   * - @f$ K @f$ is the reference element
   * - @f$ \tau @f$ is the physical element
   *
   * # Implementations
   *
   * Concrete implementations include:
   * - IdentityTransformation: @f$ x(r) = r @f$ (identity map)
   * - ParametricTransformation: Uses finite element basis functions
   *
   * # Thread Safety
   *
   * Concurrent read access is supported. Mutating a derived map requires
   * exclusive access and fresh point/quadrature and cell bindings before
   * evaluation resumes. Geometry must remain fixed during each binding.
   *
   * @see <a href="class_rodin_1_1_geometry_1_1_identity_transformation.html">IdentityTransformation</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_parametric_transformation.html">ParametricTransformation</a>
   * @see <a href="class_rodin_1_1_geometry_1_1_point.html">Point</a>
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
      constexpr PolytopeTransformation(size_t rdim, size_t pdim)
        : m_rdim(rdim),
          m_pdim(pdim)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      constexpr PolytopeTransformation(const PolytopeTransformation& other) noexcept
        : Copyable(other),
          m_rdim(other.m_rdim),
          m_pdim(other.m_pdim)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      constexpr PolytopeTransformation(PolytopeTransformation&& other) noexcept
        : Copyable(std::move(other)),
          m_rdim(other.m_rdim),
          m_pdim(other.m_pdim)
      {}

      /**
       * @brief Move assignment operator.
       * @returns Reference to this object after the operation.
       * @param other Object to move from.
       */
      PolytopeTransformation& operator=(PolytopeTransformation&& other) noexcept
      {
        if (this != &other)
        {
          m_rdim = other.m_rdim;
          m_pdim = other.m_pdim;
        }
        return *this;
      }

      /// @brief Virtual destructor.
      virtual ~PolytopeTransformation() = default;

      /**
       * @brief Gets the reference dimension @f$ k @f$.
       * @returns Reference dimension (topological dimension of reference element)
       */
      constexpr size_t getReferenceDimension() const
      {
        return m_rdim;
      }

      /**
       * @brief Gets the physical dimension @f$ s @f$.
       * @returns Physical dimension (spatial dimension of physical element)
       */
      constexpr size_t getPhysicalDimension() const
      {
        return m_pdim;
      }

      /**
       * @brief Returns a total polynomial degree bound for coordinate components.
       * @returns A certified upper bound, or no established polynomial bound.
       */
      virtual Optional<size_t> getDegree() const
      {
        return std::nullopt;
      }

      /**
       * @brief Bounds coordinate-map degree in selected reference coordinates.
       *
       * The degree is total in the selected coordinates, with the other
       * coordinates treated as coefficients. The default total-degree bound
       * is conservative for every subset. Indices are zero-based, duplicates
       * are ignored, and the empty set has degree zero.
       * @param directions Reference-coordinate indices.
       * @returns A conservative bound, or no bound when unavailable or invalid.
       */
      virtual Optional<size_t> getDegree(std::initializer_list<size_t> directions) const
      {
        for (const size_t direction : directions)
        {
          if (direction >= getReferenceDimension())
            return std::nullopt;
        }
        if (directions.size() == 0)
          return 0;
        return getDegree();
      }

      /**
       * @brief Certifies membership in the native vertex-defined geometry space.
       * @returns True when certified; false also includes unsupported maps.
       */
      virtual bool isFirstOrder() const
      {
        const auto degree = getDegree();
        return degree && *degree <= 1;
      }

      /**
       * @brief Constructs a conservative enclosure of the reference-cell image.
       * @returns An owned enclosure, or no certified enclosure for this representation.
       */
      virtual Optional<Bounds> getBounds() const
      {
        return std::nullopt;
      }

      /**
       * @brief Computes the physical coordinates from reference coordinates.
       * @param[out] pc Physical coordinates @f$ p @f$ (resized automatically)
       * @param[in] rc Reference coordinates @f$ r \in K @f$
       *
       * Computes the following quantity.
       * @f[
       *    pc = x(rc)
       * @f]
       * The output vector is resized to the physical dimension @f$ s @f$.
       */
      virtual void transform(
        Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const = 0;

      /**
       * @brief Computes the Jacobian matrix of the transformation.
       * @param[out] jacobian Jacobian matrix @f$ \mathbf{J}_x @f$ (resized automatically)
       * @param[in] rc Reference coordinates @f$ r \in K @f$
       *
       * Computes the Jacobian matrix.
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
      virtual void jacobian(
        Math::SpatialMatrix<Real>& jacobian, const Math::SpatialPoint& rc) const = 0;

      /**
       * @brief Computes the Jacobians at a reference quadrature's points.
       *
       * The output contains @f$J_x(\hat x_q)@f$ in formula order. The default
       * implementation evaluates the pointwise operation. Parametric elements
       * may reuse their reference tabulation while retaining cell-dependent
       * control points. This bulk operation is local and stores no cache.
       *
       * @param[out] jacobians Jacobian matrices, resized to the formula size.
       * @param[in] qf Reference quadrature formula defining the points.
       */
      virtual void jacobian(std::vector<Math::SpatialMatrix<Real>>& jacobians,
        const QF::QuadratureFormula& qf) const;

      /**
       * @brief Computes the reference coordinates from physical coordinates.
       * @param[out] rc Reference coordinates @f$ r @f$ (resized automatically)
       * @param[in] pc Physical coordinates @f$ p \in \tau @f$
       *
       * Computes the inverse transformation.
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
       * @param version Boost.Serialization class version; unused by this implementation.
       */
      template <class Archive>
      void serialize(Archive& ar, [[maybe_unused]] const unsigned int version)
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

      /**
       * @brief Maps an owned reference sample using available reference data.
       *
       * The caller must establish that @p rc is the sample at @p qp in the
       * persistent @p formula. Cache misses retain direct pointwise evaluation.
       * No reference-table storage is retained by the point.
       *
       * @param[out] pc Physical coordinates.
       * @param rc Owned reference coordinates of the sample.
       * @param formula Persistent immutable reference formula.
       * @param qp Logical reference-sample index.
       */
      virtual void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc,
        const QF::QuadratureFormula& formula, size_t qp) const;

    private:
      size_t m_rdim; ///< Reference dimension @f$ k @f$
      size_t m_pdim; ///< Physical dimension @f$ s @f$
  };
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Rodin::Geometry::PolytopeTransformation)

#endif
