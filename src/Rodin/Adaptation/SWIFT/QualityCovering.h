/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_QUALITYCOVERING_H
#define RODIN_ADAPTATION_SWIFT_QUALITYCOVERING_H

#include <cassert>
#include <cmath>
#include <memory>
#include <mutex>
#include <vector>

#include "QualityLattice.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Cached closed-form reference coverings with equal positive base weights.
   *
   * The integer indices of the subdivision lattice are shifted and contracted.
   * Except on tetrahedra, coordinates are @f$(I+\frac12\boldsymbol1)/(m+1)@f$.
   * Tetrahedra use @f$(I+a\boldsymbol1)/(m+3a)@f$, with @f$a=1/3@f$ for
   * @f$m=1@f$ and @f$a=1/(2\sqrt2)@f$ otherwise. Point counts are unchanged.
   * These structured coverings have known reference radii; unrestricted
   * optimality is not claimed. Vertices are not supplemental witnesses.
   *
   * @par Architecture
   * Generation visits the reference geometry's integer index set once and
   * applies its closed-form contraction. Immutable formulas are cached by
   * geometry and subdivision. Equal base weights sum to reference volume.
   * The quadrature interface supplies mapped-point caching, not polynomial
   * exactness or a continuous quality certificate.
   */
  class QualityCovering final : public QF::QuadratureFormulaBase
  {
    public:
      /**
       * @brief Gets immutable witnesses retained until program termination.
       * @param geometry Reference-cell geometry.
       * @param subdivisions Positive number of subdivisions per reference edge.
       * @returns The immutable cached covering.
       */
      static const QualityCovering& get(Geometry::Polytope::Type geometry,
        size_t subdivisions)
      {
        assert(subdivisions > 0);
        using Key = std::pair<Geometry::Polytope::Type, size_t>;
        static std::mutex mutex;
        static FlatMap<Key, std::unique_ptr<QualityCovering>> cache;
        const std::lock_guard lock(mutex);
        const Key key{geometry, subdivisions};
        auto& formula = cache[key];
        if (!formula)
          formula.reset(new QualityCovering(geometry, subdivisions));
        return *formula;
      }

      /**
       * @brief Gets the number of distinct witnesses.
       * @returns The covering cardinality.
       */
      size_t getSize() const override { return m_points.size(); }

      /**
       * @brief Gets the reference coordinates of a witness.
       * @param i Witness index.
       * @returns Reference coordinates retained by this covering.
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      { return m_points[i]; }

      /**
       * @brief Gets the equal reference penalty weight.
       * @param i Witness index; unused because base weights are equal.
       * @returns Reference-cell volume divided by covering cardinality.
       */
      Real getWeight([[maybe_unused]] size_t i) const override { return m_weight; }

      /**
       * @brief Gets the Euclidean covering radius in reference coordinates.
       * @returns The exact radius of this structured construction.
       */
      Real getCoveringRadius() const { return m_radius; }

      /**
       * @brief Copies the witness set with a fresh cache identity.
       * @returns An owning pointer to the copy.
       */
      QualityCovering* copy() const noexcept override
      { return new QualityCovering(*this); }

    private:
      /**
       * @brief Constructs the shifted and contracted generating lattice.
       * @param geometry Reference-cell geometry.
       * @param m Positive generating-lattice subdivision count.
       */
      QualityCovering(Geometry::Polytope::Type geometry, size_t m)
      {
        using G = Geometry::Polytope::Type;
        const size_t dimension = Geometry::Polytope::Traits(geometry).getDimension();
        const bool tetrahedron = geometry == G::Tetrahedron;
        const Real shift = tetrahedron
          ? m == 1 ? Real(1) / 3 : Real(1) / (2 * std::sqrt(Real(2)))
          : Real(0.5);
        const Real scale = tetrahedron ? m + 3 * shift : m + Real(1);
        m_radius = tetrahedron && m == 1 ? Real(1) / std::sqrt(Real(6))
          : std::sqrt(Real(dimension)) / (2 * scale);
        const auto& lattice = QualityLattice::get(geometry, m);
        m_points.reserve(lattice.getSize());
        for (size_t q = 0; q < lattice.getSize(); ++q)
        {
          Math::SpatialVector<Real> point(lattice.getPoint(q));
          for (size_t d = 0; d < dimension; ++d)
            point[d] = (m * point[d] + shift) / scale;
          m_points.push_back(std::move(point));
        }
        m_weight = lattice.getWeight(0);
      }

      std::vector<Math::SpatialVector<Real>> m_points;
      Real m_weight;
      Real m_radius;
  };
}

#endif
