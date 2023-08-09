/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_TRANSFORMATION_H
#define RODIN_GEOMETRY_TRANSFORMATION_H

#include "Rodin/Math.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Geometry/Simplex.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents the transformation function of a simplex, taking
   * reference coordinates to physical coordinates.
   *
   * Let @f$ \tau @f$ denote a polytope in a triangulation @f$ \mathcal{T}_h
   * @f$ which has an associated reference element @f$ K @f$, This class
   * represents the transformation @f$ x : K \subset \mathbb{R}^k \rightarrow
   * \tau \subset \mathbb{R}^s @f$ of a reference point @f$ r @f$ into a
   * physical point @f$ p @f$:
   * @f[
   *    p = x(r)
   * @f]
   * Here, @f$ k @f$ and @f$ s @f$ represent the reference and physical
   * dimensions, @f$ p \in \tau @f$ denotes the physical coordinates of the
   * point, while @f$ x : K \rightarrow \tau @f$ represents the transformation
   * taking reference coordinates @f$ r \in K @f$, for a reference geometry @f$
   * K @f$.
   *
   * @see @ref Geometry::Point "Point"
   */
  class PolytopeTransformation
  {
    public:
      constexpr
      PolytopeTransformation(size_t rdim, size_t pdim)
        : m_rdim(rdim), m_pdim(pdim)
      {}

      constexpr
      PolytopeTransformation(const PolytopeTransformation&) = default;

      constexpr
      PolytopeTransformation(PolytopeTransformation&&) = default;

      virtual ~PolytopeTransformation() = default;

      inline
      constexpr
      size_t getReferenceDimension() const
      {
        return m_rdim;
      }

      inline
      constexpr
      size_t getPhysicalDimension() const
      {
        return m_pdim;
      }

      /**
       * @brief Computes the physical coordinates of the given reference point.
       *
       * Given @f$ r \in K @f$, computes the point:
       * @f[
       *    p = x(r)
       * @f]
       * in physical coordinates.
       *
       * @param[in] rc Reference coordinates of the point.
       * @returns A vector of size @f$ s @f$ where @f$ s @f$ represents the
       * physical dimension.
       */
      Math::SpatialVector transform(const Math::SpatialVector& rc) const
      {
        Math::SpatialVector res;
        transform(rc, res);
        return res;
      }

      virtual void transform(const Math::SpatialVector& rc, Math::SpatialVector&) const = 0;

      /**
       * @brief Computes the Jacobian matrix of the transformation.
       *
       * Given @f$ r \in K @f$, computes the Jacobian matrix:
       * @f[
       *  \mathbf{J}_x (r) = \begin{bmatrix}
       * \dfrac{\partial x_1}{\partial r_1} & \ldots & \dfrac{\partial x_s}{\partial r_k}\\
       * \vdots & \ddots & \vdots\\
       * \dfrac{\partial x_s}{\partial r_1} & \ldots & \dfrac{\partial x_s}{\partial r_k}
       * \end{bmatrix} ,
       * @f]
       * for the given transformation @f$ x : K \rightarrow \tau @f$.
       *
       * @returns A matrix of dimensions @f$ s \times k @f$ where @f$ k @f$
       * represents the reference dimension and @f$ s @f$ represents the
       * physical dimension.
       */
      Math::SpatialMatrix jacobian(const Math::SpatialVector& rc) const
      {
        Math::SpatialMatrix res;
        jacobian(rc, res);
        return res;
      }

      virtual void jacobian(const Math::SpatialVector& rc, Math::SpatialMatrix& jacobian) const = 0;

      /**
       * @brief Computes the reference coordinates of the given physical point.
       *
       * Given @f$ p \in \tau @f$, computes the point:
       * @f[
       *    r = x^{-1}(p)
       * @f]
       * in reference coordinates.
       *
       * @param[in] pc Physical coordinates of the point.
       */
      Math::SpatialVector inverse(const Math::SpatialVector& pc) const
      {
        Math::SpatialVector res;
        inverse(pc, res);
        return res;
      }

      virtual void inverse(const Math::SpatialVector& pc, Math::SpatialVector& rc) const;

    private:
      const size_t m_rdim;
      const size_t m_pdim;
  };
}

#endif
