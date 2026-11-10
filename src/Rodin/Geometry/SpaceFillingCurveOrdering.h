/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SPACEFILLINGCURVEORDERING_H
#define RODIN_GEOMETRY_SPACEFILLINGCURVEORDERING_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Mesh.h"

namespace Rodin::Geometry::Internal
{
  /**
   * @brief Bounding-box quantization and sorting for spatial ordering policies.
   * @tparam Derived Policy transforming integer coordinates before bit interleaving.
   *
   * The stages are coordinate validation, bounding-box normalization, integer
   * quantization, policy transformation and stable sorting of interleaved keys.
   * No mesh incidence is computed or modified.
   */
  template <class Derived>
  class SpaceFillingCurveOrdering
  {
    public:
      SpaceFillingCurveOrdering()
        : m_bits(std::nullopt)
      {}

      /**
       * @brief Sets the quantization precision.
       * @param[in] bits Bits per coordinate, or std::nullopt for automatic precision.
       * @returns Reference to the concrete ordering policy.
       *
       * Explicit precision must be positive and at most 52. Application to a
       * mesh additionally requires at most 63 bits across all coordinates.
       * Automatic precision is the smaller of 52 and 63 divided by the space
       * dimension, rounded down, preserving Mesh::renumber()'s Morton precision.
       */
      Derived& setBitsPerCoordinate(Optional<size_t> bits)
      {
        constexpr size_t maxBits = std::numeric_limits<Real>::digits - 1;
        if (bits && (*bits == 0 || *bits > maxBits))
        {
          Alert::MemberFunctionException(static_cast<const Derived&>(*this), __func__)
            << "Bits per coordinate must be between 1 and " << maxBits << "."
            << Alert::Raise;
        }
        m_bits = bits;
        return static_cast<Derived&>(*this);
      }

      /**
       * @brief Gets the requested precision.
       * @returns Explicit bit depth, or std::nullopt for automatic precision.
       */
      Optional<size_t> getBitsPerCoordinate() const
      {
        return m_bits;
      }

      /**
       * @brief Computes a spatial vertex ordering.
       * @param[in] mesh Local mesh with finite vertex coordinates.
       * @returns Permutation with newToOld[newIndex] equal to the old vertex index.
       *
       * Each coordinate is normalized independently in the mesh bounding box.
       * Constant coordinate directions quantize to zero. Equal keys retain
       * the original vertex order, including coincident vertices.
       */
      IndexVector operator()(const Mesh<Context::Local>& mesh) const
      {
        constexpr size_t keyBits = std::numeric_limits<std::uint64_t>::digits - 1;
        constexpr size_t coordinateBits = std::numeric_limits<Real>::digits - 1;
        const size_t dim = mesh.getSpaceDimension();
        const size_t count = mesh.getVertexCount();
        const size_t bits =
          m_bits.value_or(dim > 0 ? std::min(keyBits / dim, coordinateBits) : 0);
        if (dim > 0 && bits > keyBits / dim)
        {
          Alert::MemberFunctionException(static_cast<const Derived&>(*this), __func__)
            << "The requested bit depth exceeds the 63-bit spatial key capacity."
            << Alert::Raise;
        }
        std::vector<Real> lower(dim, std::numeric_limits<Real>::max());
        std::vector<Real> upper(dim, std::numeric_limits<Real>::lowest());
        for (Index v = 0; v < count; ++v)
        {
          const auto point = mesh.getVertexCoordinates(v);
          for (size_t i = 0; i < dim; ++i)
          {
            if (!std::isfinite(point(i)))
            {
              Alert::MemberFunctionException(static_cast<const Derived&>(*this), __func__)
                << "Spatial ordering requires finite vertex coordinates." << Alert::Raise;
            }
            lower[i] = std::min(lower[i], point(i));
            upper[i] = std::max(upper[i], point(i));
          }
        }
        const Real cells = std::ldexp(Real(1), static_cast<int>(bits)) - 1;
        std::vector<std::uint64_t> codes(count, 0), coordinates(dim);
        for (Index v = 0; v < count; ++v)
        {
          const auto point = mesh.getVertexCoordinates(v);
          for (size_t i = 0; i < dim; ++i)
          {
            const Real extent = upper[i] - lower[i];
            Real scaled = extent > 0 ? (point(i) - lower[i]) / extent : 0;
            // Halving avoids overflow for finite endpoints of opposite sign.
            if (!std::isfinite(extent))
              scaled = (point(i) / 2 - lower[i] / 2) / (upper[i] / 2 - lower[i] / 2);
            coordinates[i] = static_cast<std::uint64_t>(
              std::floor(std::clamp(scaled, Real(0), Real(1)) * cells));
          }
          if (dim > 0 && bits > 0)
            Derived::transform(coordinates, bits);
          for (size_t b = bits; b-- > 0;)
          {
            for (size_t i = 0; i < dim; ++i)
              codes[v] = (codes[v] << 1) | ((coordinates[i] >> b) & 1);
          }
        }
        IndexVector order(count);
        std::iota(order.begin(), order.end(), Index(0));
        std::stable_sort(order.begin(), order.end(),
          [&](Index a, Index b) { return codes[a] < codes[b]; });
        return order;
      }

    private:
      Optional<size_t> m_bits;
  };
}

#endif
