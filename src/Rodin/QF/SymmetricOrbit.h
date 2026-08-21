/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_SYMMETRICORBIT_H
#define RODIN_QF_SYMMETRICORBIT_H

/**
 * @file
 * @brief Defines SymmetricOrbit, the barycentric symmetry orbit of a
 * fully symmetric simplex quadrature rule.
 */

#include <array>
#include <vector>
#include <algorithm>

#include "Rodin/Math.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief One symmetry orbit of a fully symmetric simplex quadrature rule.
   *
   * A rule on a simplex is *fully symmetric* when its point set is invariant
   * under every permutation of the barycentric coordinates. Such a rule is
   * therefore a union of orbits, each orbit being the set of distinct
   * permutations of a single barycentric tuple, and every point of an orbit
   * carrying the same weight.
   *
   * This is the representation used by the rule families of Xiao and Gimbutas
   * @cite xiao2010numerical and of Witherden and Vincent
   * @cite witherden2015identification: an orbit is stored once, as a tuple and
   * a weight, and expanded on construction. It is compact, it makes the
   * symmetry structural rather than incidental, and it reduces the unknowns of
   * the moment equations from @f$ d \cdot n + n @f$ to a handful of orbit
   * parameters.
   *
   * The orbit classes conventionally named @f$ S_3 @f$, @f$ S_{21}(a) @f$ and
   * @f$ S_{111}(a,b) @f$ on the triangle, and @f$ S_4 @f$, @f$ S_{31}(a) @f$,
   * @f$ S_{22}(a) @f$, @f$ S_{211}(a,b) @f$, @f$ S_{1111}(a,b,c) @f$ on the
   * tetrahedron, are not enumerated here. They are exactly the multiplicity
   * patterns of the stored tuple, so expanding the distinct permutations
   * recovers each class and its cardinality without a separate taxonomy.
   *
   * @see expand
   */
  class SymmetricOrbit
  {
    public:
      /// @brief Barycentric coordinates; @f$ d + 1 @f$ entries for a
      /// @f$ d @f$-simplex.
      using Barycentric = std::vector<Real>;

      /**
       * @brief Constructs an orbit from a barycentric tuple and its weight.
       * @param barycentric Representative tuple, summing to one.
       * @param weight Weight carried by *each* point of the orbit, not by the
       * orbit as a whole.
       */
      SymmetricOrbit(Barycentric barycentric, Real weight)
        : m_barycentric(std::move(barycentric)),
          m_weight(weight)
      {}

      /// @brief The representative barycentric tuple.
      const Barycentric& getBarycentric() const
      {
        return m_barycentric;
      }

      /// @brief The weight carried by each point of the orbit.
      Real getWeight() const
      {
        return m_weight;
      }

      /**
       * @brief The distinct permutations of the representative tuple.
       *
       * The cardinality is the multinomial coefficient of the multiplicity
       * pattern: 1, 3 or 6 on the triangle, and 1, 4, 6, 12 or 24 on the
       * tetrahedron, recovering the conventional orbit classes.
       */
      std::vector<Barycentric> expand() const
      {
        Barycentric sorted = m_barycentric;
        std::sort(sorted.begin(), sorted.end());
        std::vector<Barycentric> out;
        do
        {
          out.push_back(sorted);
        }
        while (std::next_permutation(sorted.begin(), sorted.end()));
        return out;
      }

      /// @brief Number of points contributed by this orbit.
      size_t getSize() const
      {
        return expand().size();
      }

      /**
       * @brief Maps a barycentric tuple to reference coordinates of @p g.
       *
       * The image is @f$ \sum_i \lambda_i v_i @f$ over the reference vertices
       * published by Geometry::Polytope::Traits, so the mapping follows the
       * library's reference element rather than restating it.
       */
      static Math::SpatialVector<Real> toReference(
        Geometry::Polytope::Type g, const Barycentric& barycentric)
      {
        const Geometry::Polytope::Traits traits(g);
        const size_t d = traits.getDimension();
        assert(barycentric.size() == traits.getVertexCount());
        Math::SpatialVector<Real> x;
        x.resize(static_cast<Eigen::Index>(d));
        x.setZero();
        for (size_t i = 0; i < barycentric.size(); ++i)
        {
          const auto& v = traits.getVertex(i);
          for (size_t k = 0; k < d; ++k)
          {
            x[static_cast<Eigen::Index>(k)] +=
              barycentric[i] * v[static_cast<Eigen::Index>(k)];
          }
        }
        return x;
      }

    private:
      Barycentric m_barycentric;
      Real m_weight;
  };
}

#endif
