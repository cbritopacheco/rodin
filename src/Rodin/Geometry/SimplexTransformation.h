/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_TRANSFORMATION_H
#define RODIN_GEOMETRY_TRANSFORMATION_H

#include <mfem.hpp>

#include "Rodin/Math/Vector.h"
#include "Rodin/Math/Matrix.h"

#include "ForwardDecls.h"

namespace Rodin::Geometry
{
  /**
   * @brief Represents the transformation function of a simplex, taking
   * reference coordinates to physical coordinates.
   *
   * This class represents the transformation @f$ x : K \rightarrow \tau @f$ of
   * a point:
   * @f[
   *    p = x(r)
   * @f]
   * on some simplex @f$ \tau \in \mathcal{T}_h @f$ belonging to
   * some discrete mesh @f$ \mathcal{T}_h @f$. Here @f$ p \in \tau @f$ denotes
   * the physical coordinates of the point, while @f$ x : K \rightarrow \tau
   * @f$ represents the transformation taking reference coordinates @f$ r \in K
   * @f$, for a reference geometry @f$ K @f$.
   *
   * @see Geometry::Point
   */
  class Transformation
  {
    public:
      /**
       * @brief Performs the transformation, taking reference coordinates into
       * physical coordinates.
       *
       * Given @f$ r \in K @f$, computes the point:
       * @f[
       *    p = x(r)
       * @f]
       * in physical coordinates.
       *
       * @param[in] rc Reference coordinates of the point.
       * @returns Physical coordinates
       */
      virtual Math::Vector transform(const Math::Vector& rc) const = 0;

      /**
       * @brief Performs the inverse transformation, taking physical
       * coordinates into reference coordinates.
       *
       * Given @f$ p \in \tau @f$, computes the point:
       * @f[
       *    r = x^{-1}(p)
       * @f]
       * in reference coordinates.
       *
       * @param[in] pc Physical coordinates of the point.
       */
      virtual Math::Vector inverse(const Math::Vector& pc) const = 0;

      virtual Math::Matrix jacobian(const Math::Vector& rc) const = 0;

      virtual const mfem::ElementTransformation& getHandle() const = 0;
  };

  /**
   * @brief Represents a transformation between the reference space to the
   * physical space of a simplex.
   */
  class IsoparametricTransformation final : public Transformation
  {
    public:
      IsoparametricTransformation(const Simplex& simplex)
      {

      }

      inline
      Math::Vector transform(const Math::Vector& rc) const final override
      {
        assert(false);
      }

      inline
      Math::Vector inverse(const Math::Vector& pc) const final override
      {
        assert(false);
      }

      inline
      Math::Matrix jacobian(const Math::Vector& rc) const final override
      {
        assert(false);
      }

      inline
      const mfem::IsoparametricTransformation& getHandle() const final override
      {
        return *m_handle;
      }

    private:
      const mfem::IsoparametricTransformation* m_handle;
  };
}

#endif
