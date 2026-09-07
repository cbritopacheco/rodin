/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_LOCATION_AABB_H
#define RODIN_LOCATION_AABB_H

#include <cmath>
#include <array>
#include <atomic>
#include <mutex>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <vector>
#include <map>
#include <utility>
#include <functional>

#include "Rodin/Types.h"
#include "Rodin/Geometry.h"

#include "ForwardDecls.h"

namespace Rodin::Location
{
  /**
   * @brief Bounding-volume-hierarchy point locator over axis-aligned boxes.
   *
   * @section AABBArchitecture Architecture
   *
   * The locator answers the following geometric query. Given a physical point
   * @f$ x @f$, find a polytope @f$ K @f$ and a reference coordinate
   * @f$ \widehat x \in \widehat K @f$ such that
   * @f[
   *   x = T_K(\widehat x),
   * @f]
   * where @f$ T_K : \widehat K \to K @f$ is the polytope transformation. The
   * returned value is a Geometry::Point carrying the polytope, the recovered
   * reference coordinate and the physical coordinate.
   *
   * The implementation separates the query into two stages:
   *
   * - Broad phase. For each requested polytope dimension, a packed
   *   median-split bounding-volume hierarchy is built lazily over
   *   per-polytope axis-aligned boxes @f$ B_K @f$ satisfying
   *   @f[
   *     T_K(\widehat K) \subset B_K .
   *   @f]
   *   Tree nodes reject groups of polytopes, and leaf entries reject
   *   individual polytopes, using only componentwise box containment.
   *
   * - Narrow phase. For each surviving candidate, the actual transformation
   *   is inverted by Newton iteration. The candidate is accepted only when
   *   the physical residual is below tolerance and the recovered reference
   *   coordinate lies in the reference polytope. Consequently, the AABB never
   *   certifies membership; it only decides which candidates are worth testing.
   *
   * Boxes bound the whole curved image, not just sampled points on it. The
   * transformation is resampled on a unisolvent reference lattice and
   * converted into control points of a non-negative partition of unity
   * (Bernstein on simplices, tensor Bernstein on quadrilaterals, hexahedra
   * and wedges, the collapsed-coordinate basis on pyramids), whose extrema
   * bound the image exactly for every polytope type and order. Degree-one
   * factors need no conversion: the mapped vertices are already the control
   * points, so affine and multilinear cells keep the cheap vertex box. Should
   * a geometry ever fall outside this construction the box degrades to a
   * heuristic inflation, and an exhaustive narrow-phase fallback can be
   * enabled with setExhaustiveFallback(). The degree is taken from
   * PolytopeTransformation::getOrder(), so a transformation that understates
   * its own order understates its box.
   *
   * Tolerances are relative to the mesh bounding-box diagonal. Queries are
   * thread-safe: the lazy per-dimension build is guarded by a mutex and
   * transformations are immutable after construction.
   */
  template <class MeshType>
  class AABB
  {
    public:
      /// @brief Builds a locator bound to a fixed mesh.
      explicit AABB(const MeshType& mesh)
        : m_mesh(mesh),
          m_tolerance(Real(1e-10)),
          m_referenceTolerance(Real(1e-10)),
          m_maxNewtonIterations(16),
          m_exhaustiveFallback(false),
          m_index(mesh.getDimension() + 1)
      {
        computeScale();
      }

      /// Relative physical tolerance (scaled by the mesh diagonal).
      Real getTolerance() const
      {
        return m_tolerance;
      }

      /// Sets the relative physical tolerance and invalidates the index.
      AABB& setTolerance(Real tolerance)
      {
        m_tolerance = tolerance;
        invalidate();
        return *this;
      }

      /// Reference-space containment slack (reference coordinates are O(1)).
      Real getReferenceTolerance() const
      {
        return m_referenceTolerance;
      }

      /// @brief Sets the reference-space containment slack.
      AABB& setReferenceTolerance(Real tolerance)
      {
        m_referenceTolerance = tolerance;
        return *this;
      }

      /**
       * @brief Enables the exhaustive narrow-phase fallback on broad-phase
       * miss.
       *
       * Only useful for transformations of order three or higher, whose true
       * extent may exceed the sampled, inflated boxes. Costs one full
       * narrow-phase sweep per miss.
       */
      AABB& setExhaustiveFallback(bool fallback)
      {
        m_exhaustiveFallback = fallback;
        return *this;
      }

      /**
       * @brief Locates a physical point on a polytope of the requested dimension.
       *
       * Returns an empty optional when the point is outside every candidate
       * polytope, when the coordinate dimension is incompatible with the mesh,
       * or when the inverse transformation does not pass the residual checks.
       */
      Optional<Geometry::Point> locate(
        size_t dimension, const Math::SpatialPoint& x) const
      {
        if (dimension >= m_index.size())
          return {};
        if (static_cast<size_t>(x.size()) != m_mesh.get().getSpaceDimension())
          return {};
        if (!isFinite(x))
          return {};

        const DimensionIndex& index = ensureBuilt(dimension);
        if (auto p = traverse(index, dimension, x))
          return p;
        if (m_exhaustiveFallback)
          return exhaustive(index, dimension, x);
        return {};
      }

      /// @brief Locates a physical point on a cell of the mesh dimension.
      Optional<Geometry::Point> locate(const Math::SpatialPoint& x) const
      {
        return locate(m_mesh.get().getDimension(), x);
      }

    private:
      static constexpr size_t MaxSpaceDimension = 3;
      static constexpr size_t LeafSize = 8;
      static constexpr int32_t StackDepth = 64;

      static bool isFinite(const Math::SpatialPoint& v)
      {
        for (Eigen::Index i = 0; i < v.size(); ++i)
        {
          if (!std::isfinite(v[i]))
            return false;
        }
        return true;
      }

      static bool isFinite(const Math::SpatialMatrix<Real>& m)
      {
        for (Eigen::Index i = 0; i < m.rows(); ++i)
        {
          for (Eigen::Index j = 0; j < m.cols(); ++j)
          {
            if (!std::isfinite(m(i, j)))
              return false;
          }
        }
        return true;
      }

      using Bound = std::array<Real, MaxSpaceDimension>;

      struct Node
      {
          Bound lo;
          Bound hi;
          int32_t left;    ///< Left child, or -1 for a leaf
          int32_t right;   ///< Right child (valid when left >= 0)
          uint32_t begin;  ///< First entry (valid for leaves)
          uint32_t end;    ///< One past last entry (valid for leaves)
      };

      struct DimensionIndex
      {
          std::vector<Node> nodes;
          std::vector<Index> entries;      ///< Polytope indices, leaf-ordered
          std::vector<Bound> entryLo;      ///< Per-entry box, leaf-ordered
          std::vector<Bound> entryHi;
          std::atomic<bool> built{false};
          mutable std::mutex mutex;
      };

      void invalidate()
      {
        for (auto& index : m_index)
        {
          std::lock_guard lock(index.mutex);
          index.nodes.clear();
          index.entries.clear();
          index.entryLo.clear();
          index.entryHi.clear();
          index.built.store(false, std::memory_order_release);
        }
        computeScale();
      }

      void computeScale()
      {
        const auto& mesh = m_mesh.get();
        const size_t sdim = mesh.getSpaceDimension();
        Bound lo, hi;
        lo.fill(std::numeric_limits<Real>::infinity());
        hi.fill(-std::numeric_limits<Real>::infinity());
        for (Index v = 0; v < mesh.getVertexCount(); ++v)
        {
          const auto& x = mesh.getVertexCoordinates(v);
          for (size_t i = 0; i < sdim; ++i)
          {
            lo[i] = std::min(lo[i], x[static_cast<Eigen::Index>(i)]);
            hi[i] = std::max(hi[i], x[static_cast<Eigen::Index>(i)]);
          }
        }
        Real diag2 = 0;
        for (size_t i = 0; i < sdim; ++i)
        {
          const Real e = hi[i] - lo[i];
          if (std::isfinite(e))
            diag2 += e * e;
        }
        const Real diag = std::sqrt(diag2);
        m_scale = diag > Real(0) ? diag : Real(1);
      }

      /// Effective physical tolerance in mesh units.
      Real physicalTolerance() const
      {
        return m_tolerance * m_scale;
      }

      const DimensionIndex& ensureBuilt(size_t dimension) const
      {
        DimensionIndex& index = m_index[dimension];
        if (!index.built.load(std::memory_order_acquire))
        {
          std::lock_guard lock(index.mutex);
          if (!index.built.load(std::memory_order_relaxed))
          {
            build(index, dimension);
            index.built.store(true, std::memory_order_release);
          }
        }
        return index;
      }

      void build(DimensionIndex& index, size_t dimension) const
      {
        const auto& mesh = m_mesh.get();
        const size_t sdim = mesh.getSpaceDimension();
        const size_t count = mesh.getPolytopeCount(dimension);

        index.entries.clear();
        index.nodes.clear();
        index.entryLo.clear();
        index.entryHi.clear();
        if (count == 0)
          return;

        // Per-entry boxes and centroids in flat storage.
        std::vector<Bound> lo(count), hi(count);
        std::vector<Bound> mid(count);
        index.entries.reserve(count);
        size_t n = 0;
        for (auto it = mesh.getPolytope(dimension); it; ++it, ++n)
        {
          makeBox(*it, lo[n], hi[n]);
          for (size_t i = 0; i < sdim; ++i)
            mid[n][i] = Real(0.5) * (lo[n][i] + hi[n][i]);
          index.entries.push_back(it->getIndex());
        }
        assert(n == count);

        std::vector<uint32_t> order(count);
        for (uint32_t i = 0; i < count; ++i)
          order[i] = i;

        index.nodes.reserve(2 * count / LeafSize + 2);
        buildNode(index, order, 0, static_cast<uint32_t>(count), lo, hi, mid, sdim);

        // Reorder entries to leaf order so leaves are contiguous. The boxes
        // follow, so that a leaf can reject a candidate before inverting it.
        std::vector<Index> reordered(count);
        index.entryLo.resize(count);
        index.entryHi.resize(count);
        for (size_t i = 0; i < count; ++i)
        {
          reordered[i] = index.entries[order[i]];
          index.entryLo[i] = lo[order[i]];
          index.entryHi[i] = hi[order[i]];
        }
        index.entries = std::move(reordered);
      }

      int32_t buildNode(DimensionIndex& index, std::vector<uint32_t>& order,
        uint32_t begin, uint32_t end, const std::vector<Bound>& lo,
        const std::vector<Bound>& hi, const std::vector<Bound>& mid, size_t sdim) const
      {
        const int32_t self = static_cast<int32_t>(index.nodes.size());
        index.nodes.emplace_back();
        {
          Node& node = index.nodes.back();
          node.lo.fill(std::numeric_limits<Real>::infinity());
          node.hi.fill(-std::numeric_limits<Real>::infinity());
          for (uint32_t k = begin; k < end; ++k)
          {
            for (size_t i = 0; i < sdim; ++i)
            {
              node.lo[i] = std::min(node.lo[i], lo[order[k]][i]);
              node.hi[i] = std::max(node.hi[i], hi[order[k]][i]);
            }
          }
          node.left = -1;
          node.right = -1;
          node.begin = begin;
          node.end = end;
        }

        if (end - begin <= LeafSize)
          return self;

        // Split on the widest centroid axis; fall back to a leaf when the
        // centroids are degenerate (all identical).
        Bound clo, chi;
        clo.fill(std::numeric_limits<Real>::infinity());
        chi.fill(-std::numeric_limits<Real>::infinity());
        for (uint32_t k = begin; k < end; ++k)
        {
          for (size_t i = 0; i < sdim; ++i)
          {
            clo[i] = std::min(clo[i], mid[order[k]][i]);
            chi[i] = std::max(chi[i], mid[order[k]][i]);
          }
        }
        size_t axis = 0;
        Real extent = Real(0);
        for (size_t i = 0; i < sdim; ++i)
        {
          const Real e = chi[i] - clo[i];
          if (e > extent)
          {
            extent = e;
            axis = i;
          }
        }
        if (!(extent > Real(0)))
          return self;

        const uint32_t midIdx = (begin + end) / 2;
        std::nth_element(order.begin() + begin, order.begin() + midIdx,
          order.begin() + end,
          [&](uint32_t a, uint32_t b) { return mid[a][axis] < mid[b][axis]; });

        const int32_t left = buildNode(index, order, begin, midIdx, lo, hi, mid, sdim);
        const int32_t right = buildNode(index, order, midIdx, end, lo, hi, mid, sdim);
        index.nodes[self].left = left;
        index.nodes[self].right = right;
        return self;
      }

      /**
       * @brief Per-factor polynomial degree of a geometry transformation.
       *
       * PolytopeTransformation::getOrder() reports the total degree, which on
       * a tensor-product geometry is the sum over its factors: @f$ 2k @f$ on
       * quadrilaterals, wedges and pyramids, @f$ 3k @f$ on hexahedra. The
       * control-point basis is built per factor, so the total is divided by
       * the number of factors. Rounding up keeps the recovered degree
       * conservative if a transformation ever reports a total that is not an
       * exact multiple.
       */
      static size_t factorDegree(Geometry::Polytope::Type g, size_t order)
      {
        using G = Geometry::Polytope::Type;
        switch (g)
        {
          case G::Point:
            return 0;
          case G::Segment:
          case G::Triangle:
          case G::Tetrahedron:
            return order;
          case G::Quadrilateral:
          case G::Wedge:
          case G::Pyramid:
            return (order + 1) / 2;
          case G::Hexahedron:
            return (order + 2) / 3;
        }
        assert(false);
        return order;
      }

      static Real binomial(size_t n, size_t i)
      {
        if (i > n)
          return 0;
        if (i > n - i)
          i = n - i;
        Real out = 1;
        for (size_t p = 1; p <= i; ++p)
        {
          out *= static_cast<Real>(n + 1 - p);
          out /= static_cast<Real>(p);
        }
        return out;
      }

      /// @brief Univariate Bernstein polynomial of degree @f$ n @f$ on [0, 1].
      static Real bernstein(size_t n, size_t i, Real x)
      {
        if (i > n)
          return 0;
        Real out = binomial(n, i);
        for (size_t p = 0; p < i; ++p)
          out *= x;
        for (size_t p = 0; p < n - i; ++p)
          out *= Real(1) - x;
        return out;
      }

      /// @brief Multinomial coefficient @f$ k! / \prod_i \alpha_i! @f$.
      static Real multinomial(size_t k, const std::array<size_t, 4>& alpha, size_t count)
      {
        Real out = 1;
        for (size_t p = 2; p <= k; ++p)
          out *= static_cast<Real>(p);
        for (size_t i = 0; i < count; ++i)
        {
          for (size_t p = 2; p <= alpha[i]; ++p)
            out /= static_cast<Real>(p);
        }
        return out;
      }

      /**
       * @brief Bernstein polynomial of degree @f$ k @f$ on the reference
       * simplex of dimension @p rdim, indexed by a barycentric multi-index.
       */
      static Real simplexBernstein(size_t k, const std::array<size_t, 4>& alpha,
        size_t rdim, const Math::SpatialPoint& r)
      {
        Real lambda = 1;
        for (size_t i = 0; i < rdim; ++i)
          lambda -= r[static_cast<Eigen::Index>(i)];
        Real out = multinomial(k, alpha, rdim + 1);
        for (size_t p = 0; p < alpha[0]; ++p)
          out *= lambda;
        for (size_t i = 0; i < rdim; ++i)
        {
          for (size_t p = 0; p < alpha[i + 1]; ++p)
            out *= r[static_cast<Eigen::Index>(i)];
        }
        return out;
      }

      /// @brief Evaluates one control-point basis function of degree @p k.
      static Real controlBasisValue(Geometry::Polytope::Type g, size_t k,
        const std::array<size_t, 4>& mode, const Math::SpatialPoint& r)
      {
        using G = Geometry::Polytope::Type;
        switch (g)
        {
          case G::Segment:
            return bernstein(k, mode[0], r.x());
          case G::Triangle:
            return simplexBernstein(k, mode, 2, r);
          case G::Tetrahedron:
            return simplexBernstein(k, mode, 3, r);
          case G::Quadrilateral:
            return bernstein(k, mode[0], r.x()) * bernstein(k, mode[1], r.y());
          case G::Hexahedron:
            return bernstein(k, mode[0], r.x()) * bernstein(k, mode[1], r.y())
              * bernstein(k, mode[2], r.z());
          case G::Wedge:
            return simplexBernstein(k, mode, 2, r) * bernstein(k, mode[3], r.z());
          case G::Pyramid:
          {
            // Collapsed coordinates a = x / (1 - z), b = y / (1 - z), as in
            // the pyramid element itself. The apex layer carries the whole
            // partition of unity at z = 1.
            const size_t m = mode[2];
            const size_t n = k - m;
            const Real z = r.z();
            const Real q = Real(1) - z;
            const Real bz = bernstein(k, m, z);
            if (n == 0)
              return bz;
            if (!(q > std::numeric_limits<Real>::epsilon()))
              return 0;
            return bernstein(n, mode[0], r.x() / q)
              * bernstein(n, mode[1], r.y() / q) * bz;
          }
          case G::Point:
            return 1;
        }
        assert(false);
        return 0;
      }

      /**
       * @brief Reference lattice and control-point basis indices of degree
       * @p k on @p g.
       *
       * The lattice is the equispaced (collapsed, on pyramids) node set of
       * the corresponding finite element, so it is unisolvent for the space
       * the geometry transformation lives in.
       */
      static void makeControlLattice(Geometry::Polytope::Type g, size_t k,
        std::vector<std::array<size_t, 4>>& modes,
        std::vector<Math::SpatialPoint>& samples)
      {
        using G = Geometry::Polytope::Type;
        const Real kr = static_cast<Real>(k);
        const size_t rdim = Geometry::Polytope::Traits(g).getDimension();
        auto emplace = [&](const std::array<size_t, 4>& mode,
                          std::initializer_list<Real> coords)
        {
          modes.push_back(mode);
          Math::SpatialPoint r;
          r.resize(static_cast<Eigen::Index>(rdim));
          Eigen::Index i = 0;
          for (Real c : coords)
            r[i++] = c;
          assert(i == static_cast<Eigen::Index>(rdim));
          samples.push_back(std::move(r));
        };

        switch (g)
        {
          case G::Segment:
          {
            for (size_t i = 0; i <= k; ++i)
              emplace({ i, 0, 0, 0 }, { static_cast<Real>(i) / kr });
            return;
          }
          case G::Triangle:
          {
            for (size_t a2 = 0; a2 <= k; ++a2)
            {
              for (size_t a1 = 0; a1 + a2 <= k; ++a1)
              {
                emplace({ k - a1 - a2, a1, a2, 0 },
                  { static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr });
              }
            }
            return;
          }
          case G::Tetrahedron:
          {
            for (size_t a3 = 0; a3 <= k; ++a3)
            {
              for (size_t a2 = 0; a2 + a3 <= k; ++a2)
              {
                for (size_t a1 = 0; a1 + a2 + a3 <= k; ++a1)
                {
                  emplace({ k - a1 - a2 - a3, a1, a2, a3 },
                    { static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr,
                      static_cast<Real>(a3) / kr });
                }
              }
            }
            return;
          }
          case G::Quadrilateral:
          {
            for (size_t j = 0; j <= k; ++j)
            {
              for (size_t i = 0; i <= k; ++i)
              {
                emplace({ i, j, 0, 0 },
                  { static_cast<Real>(i) / kr, static_cast<Real>(j) / kr });
              }
            }
            return;
          }
          case G::Hexahedron:
          {
            for (size_t l = 0; l <= k; ++l)
            {
              for (size_t j = 0; j <= k; ++j)
              {
                for (size_t i = 0; i <= k; ++i)
                {
                  emplace({ i, j, l, 0 },
                    { static_cast<Real>(i) / kr, static_cast<Real>(j) / kr,
                      static_cast<Real>(l) / kr });
                }
              }
            }
            return;
          }
          case G::Wedge:
          {
            for (size_t m = 0; m <= k; ++m)
            {
              for (size_t a2 = 0; a2 <= k; ++a2)
              {
                for (size_t a1 = 0; a1 + a2 <= k; ++a1)
                {
                  emplace({ k - a1 - a2, a1, a2, m },
                    { static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr,
                      static_cast<Real>(m) / kr });
                }
              }
            }
            return;
          }
          case G::Pyramid:
          {
            for (size_t m = 0; m <= k; ++m)
            {
              const size_t n = k - m;
              const Real z = static_cast<Real>(m) / kr;
              const Real q = Real(1) - z;
              if (n == 0)
              {
                emplace({ 0, 0, m, 0 }, { Real(0), Real(0), z });
                continue;
              }
              const Real nr = static_cast<Real>(n);
              for (size_t j = 0; j <= n; ++j)
              {
                for (size_t i = 0; i <= n; ++i)
                {
                  emplace({ i, j, m, 0 },
                    { q * static_cast<Real>(i) / nr, q * static_cast<Real>(j) / nr, z });
                }
              }
            }
            return;
          }
          case G::Point:
            return;
        }
        assert(false);
      }

      /**
       * @brief Reference samples and the map taking sampled coordinates to
       * control points of a geometry transformation.
       *
       * The transformation of a polytope of per-factor degree @f$ k @f$ lies
       * in the span of a basis @f$ \psi_n @f$ that is non-negative and sums
       * to one on the reference domain: Bernstein polynomials on simplices,
       * their tensor products on quadrilaterals, hexahedra and wedges, and
       * the collapsed-coordinate product basis on pyramids. Writing
       * @f[
       *   T(r) = \sum_n C_n \, \psi_n(r),
       * @f]
       * every mapped point is a convex combination of the control points
       * @f$ C_n @f$, so their componentwise extrema bound the entire curved
       * image -- exactly, not heuristically, and for every order.
       *
       * The control points are recovered from values sampled on a unisolvent
       * lattice: @f$ F_m = T(r_m) = \sum_n C_n \psi_n(r_m) @f$ reads
       * @f$ F = C V^\top @f$ with @f$ V_{mn} = \psi_n(r_m) @f$, so
       * @f$ C = F (V^{-1})^\top @f$. Only the conversion depends on the
       * geometry and the degree, so it is built once and shared by every
       * polytope of that kind.
       */
      struct ControlBasis
      {
        std::vector<Math::SpatialPoint> samples;
        Math::Matrix<Real> conversion;
        bool valid = false;
      };

      /// @brief Returns the cached control-point basis of degree @p k on @p g.
      static const ControlBasis& getControlBasis(Geometry::Polytope::Type g, size_t k)
      {
        using Key = std::pair<int, size_t>;
        static std::map<Key, ControlBasis> s_cache;
        static std::mutex s_mutex;

        const Key key{ static_cast<int>(g), k };
        std::lock_guard<std::mutex> lock(s_mutex);
        const auto it = s_cache.find(key);
        if (it != s_cache.end())
          return it->second;

        ControlBasis basis;
        std::vector<std::array<size_t, 4>> modes;
        makeControlLattice(g, k, modes, basis.samples);
        const size_t n = modes.size();
        if (n > 0)
        {
          Math::Matrix<Real> vandermonde(n, n);
          for (size_t m = 0; m < n; ++m)
          {
            for (size_t c = 0; c < n; ++c)
            {
              vandermonde(m, c) = controlBasisValue(g, k, modes[c], basis.samples[m]);
            }
          }
          const Eigen::FullPivLU<Math::Matrix<Real>> lu(vandermonde);
          if (lu.isInvertible())
          {
            basis.conversion = lu.inverse().transpose();
            basis.valid = true;
          }
        }
        return s_cache.emplace(key, std::move(basis)).first->second;
      }

      void padBox(Bound& lo, Bound& hi, Real pad, size_t sdim) const
      {
        for (size_t i = 0; i < sdim; ++i)
        {
          lo[i] -= pad;
          hi[i] += pad;
        }
      }

      void makeBox(const Geometry::Polytope& polytope, Bound& lo, Bound& hi) const
      {
        const auto& mesh = m_mesh.get();
        const size_t sdim = mesh.getSpaceDimension();
        lo.fill(std::numeric_limits<Real>::infinity());
        hi.fill(-std::numeric_limits<Real>::infinity());

        const Real physTol = physicalTolerance();
        if (polytope.getDimension() == 0)
        {
          const auto& x = mesh.getVertexCoordinates(polytope.getIndex());
          for (size_t i = 0; i < sdim; ++i)
          {
            lo[i] = x[static_cast<Eigen::Index>(i)] - physTol;
            hi[i] = x[static_cast<Eigen::Index>(i)] + physTol;
          }
          return;
        }

        const auto& transformation = polytope.getTransformation();
        const auto g = polytope.getGeometry();
        const Geometry::Polytope::Traits traits(g);
        const size_t nv = traits.getVertexCount();

        Math::SpatialPoint x;
        auto add = [&](const Math::SpatialPoint& p) {
          for (size_t i = 0; i < sdim; ++i)
          {
            lo[i] = std::min(lo[i], p[static_cast<Eigen::Index>(i)]);
            hi[i] = std::max(hi[i], p[static_cast<Eigen::Index>(i)]);
          }
        };

        // Degree one per factor: the control-point basis is the vertex
        // partition of unity itself -- barycentric on a simplex, multilinear
        // on a tensor geometry -- so the mapped vertices are the control
        // points and their box already bounds the image. Affine cells stay
        // on this path and pay nothing for the curved machinery.
        const size_t degree = factorDegree(g, transformation.getOrder());
        if (degree <= 1)
        {
          for (size_t i = 0; i < nv; ++i)
          {
            transformation.transform(x, traits.getVertex(i));
            add(x);
          }
          padBox(lo, hi, physTol, sdim);
          return;
        }

        const ControlBasis& basis = getControlBasis(g, degree);
        if (basis.valid)
        {
          const size_t n = basis.samples.size();
          Math::Matrix<Real> sampled(sdim, n);
          for (size_t m = 0; m < n; ++m)
          {
            transformation.transform(x, basis.samples[m]);
            for (size_t i = 0; i < sdim; ++i)
              sampled(i, m) = x[static_cast<Eigen::Index>(i)];
          }
          const Math::Matrix<Real> control = sampled * basis.conversion;
          for (size_t i = 0; i < sdim; ++i)
          {
            lo[i] = control.row(static_cast<Eigen::Index>(i)).minCoeff();
            hi[i] = control.row(static_cast<Eigen::Index>(i)).maxCoeff();
          }
          padBox(lo, hi, physTol, sdim);
          return;
        }

        // No conservative basis available for this geometry and degree: fall
        // back to sampled vertices, centroid and edge midpoints inflated by
        // the largest chord deviation. This bounds quadratic geometries and
        // is heuristic beyond them, which is what setExhaustiveFallback()
        // exists for.
        std::vector<Math::SpatialPoint> vertex(nv);
        for (size_t i = 0; i < nv; ++i)
        {
          transformation.transform(vertex[i], traits.getVertex(i));
          add(vertex[i]);
        }
        transformation.transform(x, traits.getCentroid());
        add(x);

        Real deviation = 0;
        for (size_t i = 0; i < nv; ++i)
        {
          for (size_t j = i + 1; j < nv; ++j)
          {
            transformation.transform(
              x, Real(0.5) * (traits.getVertex(i) + traits.getVertex(j)));
            add(x);
            deviation =
              std::max(deviation, (x - Real(0.5) * (vertex[i] + vertex[j])).norm());
          }
        }

        padBox(lo, hi, physTol + deviation, sdim);
      }

      bool boxContains(
        const Bound& lo, const Bound& hi, const Math::SpatialPoint& x, size_t sdim) const
      {
        for (size_t i = 0; i < sdim; ++i)
        {
          const Real xi = x[static_cast<Eigen::Index>(i)];
          if (xi < lo[i] || xi > hi[i])
            return false;
        }
        return true;
      }

      bool boxContains(const Node& node, const Math::SpatialPoint& x, size_t sdim) const
      {
        return boxContains(node.lo, node.hi, x, sdim);
      }

      bool containsReference(
        const Geometry::Polytope::Traits& traits, const Math::SpatialPoint& rc) const
      {
        if (traits.getDimension() == 0)
          return rc.size() == 0;
        if (!isFinite(rc))
          return false;

        const auto& hs = traits.getHalfSpace();
        for (Eigen::Index i = 0; i < hs.vector.size(); ++i)
        {
          const Real margin = hs.vector[i] - rc.dot(hs.matrix.row(i).transpose());
          if (!(margin >= -m_referenceTolerance))
            return false;
        }
        return true;
      }

      /**
       * Newton inversion of the polytope transformation. Exact after one
       * iteration for affine maps; iterative for bilinear and curved maps.
       * Returns true only when the physical residual is below tolerance, so
       * off-manifold points (facet queries) and diverged iterations are
       * rejected.
       */
      bool invert(const Geometry::PolytopeTransformation& transformation,
        const Geometry::Polytope::Traits& traits, const Math::SpatialPoint& x,
        Math::SpatialPoint& rc) const
      {
        const Real physTol = physicalTolerance();
        const size_t rdim = traits.getDimension();
        const size_t pdim = static_cast<size_t>(x.size());

        rc = traits.getCentroid();
        Math::SpatialPoint mapped;
        Math::SpatialMatrix<Real> jac;
        Math::SpatialPoint residual;
        Math::SpatialPoint step;

        for (size_t iteration = 0; iteration < m_maxNewtonIterations; ++iteration)
        {
          transformation.transform(mapped, rc);
          residual = x - mapped;
          if (residual.norm() <= physTol)
            return true;

          transformation.jacobian(jac, rc);
          if (!isFinite(jac))
            return false;

          if (rdim == pdim)
            step = jac.solve(residual);
          else
          {
            const Math::SpatialMatrix<Real> normal = jac.transpose() * jac;
            const Math::SpatialPoint rhs = jac.transpose() * residual;
            step = normal.solve(rhs);
          }
          if (!isFinite(step))
            return false;

          rc += step;
          // Diverging iterates cannot represent a contained point: reference
          // coordinates of interest live in an O(1) neighborhood.
          if (rc.norm() > Real(1e3))
            return false;
        }

        transformation.transform(mapped, rc);
        return (x - mapped).norm() <= physTol;
      }

      Optional<Geometry::Point> narrowPhase(
        size_t dimension, Index polytopeIndex, const Math::SpatialPoint& x) const
      {
        const auto& mesh = m_mesh.get();
        const Geometry::Polytope polytope = *mesh.getPolytope(dimension, polytopeIndex);

        if (dimension == 0)
        {
          if ((mesh.getVertexCoordinates(polytopeIndex) - x).norm() <=
            physicalTolerance())
            return Geometry::Point(polytope, Math::SpatialPoint(0), x);
          return {};
        }

        const Geometry::Polytope::Traits traits(polytope.getGeometry());
        Math::SpatialPoint rc;
        if (!invert(polytope.getTransformation(), traits, x, rc))
          return {};
        if (!containsReference(traits, rc))
          return {};
        return Geometry::Point(polytope, rc, x);
      }

      Optional<Geometry::Point> traverse(
        const DimensionIndex& index, size_t dimension, const Math::SpatialPoint& x) const
      {
        if (index.nodes.empty())
          return {};

        const size_t sdim = m_mesh.get().getSpaceDimension();
        std::array<int32_t, StackDepth> stack;
        int32_t top = 0;
        stack[top++] = 0;
        while (top > 0)
        {
          const Node& node = index.nodes[stack[--top]];
          if (!boxContains(node, x, sdim))
            continue;
          if (node.left < 0)
          {
            for (uint32_t k = node.begin; k < node.end; ++k)
            {
              // The entry box bounds the polytope, curvature included, so a
              // miss here rules the candidate out without a Newton inversion.
              if (!boxContains(index.entryLo[k], index.entryHi[k], x, sdim))
                continue;
              if (auto p = narrowPhase(dimension, index.entries[k], x))
                return p;
            }
            continue;
          }
          assert(top + 2 <= StackDepth);
          stack[top++] = node.left;
          stack[top++] = node.right;
        }
        return {};
      }

      Optional<Geometry::Point> exhaustive(
        const DimensionIndex& index, size_t dimension, const Math::SpatialPoint& x) const
      {
        for (const Index polytopeIndex : index.entries)
        {
          if (auto p = narrowPhase(dimension, polytopeIndex, x))
            return p;
        }
        return {};
      }

      std::reference_wrapper<const MeshType> m_mesh;
      Real m_tolerance;
      Real m_referenceTolerance;
      size_t m_maxNewtonIterations;
      bool m_exhaustiveFallback;
      Real m_scale;
      mutable std::vector<DimensionIndex> m_index;
  };
}

#endif
