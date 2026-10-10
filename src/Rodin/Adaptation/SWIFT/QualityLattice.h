/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_QUALITYLATTICE_H
#define RODIN_ADAPTATION_SWIFT_QUALITYLATTICE_H

#include <cassert>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Cached generating lattice for reference coverings and offline comparisons.
   *
   * Simplex points have barycentric coordinates in integer multiples of
   * @f$1/m@f$. Tensor cells use Cartesian grids, wedges use triangular grids
   * times segments, and pyramids use shrinking square layers. All vertices
   * are included without supplemental points. Weights sum to reference volume.
   * SWIFT samples @ref QualityCovering instead; this lattice supplies its
   * geometry-specific indices and supports offline covering comparisons.
   * The quadrature interface supplies mapped-point caching, not polynomial
   * exactness or a continuous quality certificate.
   */
  class QualityLattice final : public QF::QuadratureFormulaBase
  {
    public:
      /**
       * @brief Gets immutable witnesses retained until program termination.
       * @param geometry Reference-cell geometry.
       * @param subdivisions Positive number of subdivisions per reference edge.
       */
      static const QualityLattice& get(Geometry::Polytope::Type geometry,
        size_t subdivisions)
      {
        assert(subdivisions > 0);
        using Key = std::pair<Geometry::Polytope::Type, size_t>;
        static std::mutex mutex;
        static std::map<Key, std::unique_ptr<QualityLattice>> cache;
        const std::lock_guard lock(mutex);
        const Key key{geometry, subdivisions};
        auto& formula = cache[key];
        if (!formula)
          formula.reset(new QualityLattice(geometry, subdivisions));
        return *formula;
      }

      /// @brief Gets the number of distinct witnesses.
      size_t getSize() const override { return m_points.size(); }

      /// @brief Gets the reference coordinates of a witness.
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      { return m_points[i]; }

      /// @brief Gets the equal reference penalty weight.
      Real getWeight(size_t) const override { return m_weight; }

      /// @brief Copies the witness set with a fresh cache identity.
      QualityLattice* copy() const noexcept override
      { return new QualityLattice(*this); }

    private:
      QualityLattice(Geometry::Polytope::Type geometry, size_t m)
      {
        using G = Geometry::Polytope::Type;
        const size_t dimension = Geometry::Polytope::Traits(geometry).getDimension();
        Real volume = 1;
        switch (geometry)
        {
          case G::Triangle: volume = Real(1) / 2; break;
          case G::Tetrahedron: volume = Real(1) / 6; break;
          case G::Wedge: volume = Real(1) / 2; break;
          case G::Pyramid: volume = Real(1) / 3; break;
          default: break;
        }
        const size_t zMax = dimension == 3 ? m : 0;
        for (size_t k = 0; k <= zMax; ++k)
        {
          const size_t yMax = dimension < 2 ? 0
            : geometry == G::Tetrahedron || geometry == G::Pyramid ? m - k : m;
          for (size_t j = 0; j <= yMax; ++j)
          {
            const size_t xMax = dimension == 0 ? 0
              : geometry == G::Triangle || geometry == G::Wedge ? m - j
              : geometry == G::Tetrahedron ? m - j - k
              : geometry == G::Pyramid ? m - k : m;
            for (size_t i = 0; i <= xMax; ++i)
            {
              Math::SpatialVector<Real> point(dimension);
              if (dimension > 0)
                point[0] = Real(i) / m;
              if (dimension > 1)
                point[1] = Real(j) / m;
              if (dimension > 2)
                point[2] = Real(k) / m;
              m_points.push_back(std::move(point));
            }
          }
        }
        m_weight = volume / m_points.size();
      }

      std::vector<Math::SpatialVector<Real>> m_points;
      Real m_weight;
  };
}

#endif
