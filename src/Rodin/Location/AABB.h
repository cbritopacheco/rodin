/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_LOCATION_AABB_H
#define RODIN_LOCATION_AABB_H

#include <cmath>
#include <stdexcept>
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
#include <Eigen/QR>

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
   *   individual polytopes, using componentwise box containment. Optional
   *   control-hull projections reject further candidates before inversion.
   *
   * - Narrow phase. For each surviving candidate, the actual transformation
   *   is inverted by controlled Newton iteration, retrying from other
   *   reference points if the centroid attempt fails. A reference coordinate
   *   slightly outside the polytope is clipped back to it, then accepted only
   *   if the clipped point maps within physical tolerance of the query. The
   *   AABB never certifies membership; it only selects candidates to test.
   *
   * Boxes bound the whole curved image, not just sampled points on it. The
   * transformation is resampled on a unisolvent reference lattice and
   * converted into control points of a non-negative partition of unity
   * (Bernstein on simplices, tensor Bernstein on quadrilaterals, hexahedra
   * and wedges, the collapsed-coordinate basis on pyramids), whose extrema
   * bound the image in exact arithmetic. Numerical boxes include a conversion
   * roundoff allowance. Degree-one
   * factors need no conversion: the mapped vertices are already the control
   * points, so affine and multilinear cells keep the cheap vertex box. If a
   * box cannot be computed reliably, that entry remains unpruned. An
   * exhaustive narrow-phase fallback can also be enabled with
   * setExhaustiveFallback(). The degree is obtained from
   * PolytopeTransformation::getFactorOrder(); its default total-degree bound
   * is conservative for custom polynomial maps.
   *
   * Physical tolerance is relative to the mesh bounding-box diagonal. The
   * reference tolerance only permits a small inversion overshoot before
   * clipping; physical residual is checked again after clipping. Queries are
   * thread-safe for concurrent queries on a fixed mesh and configuration:
   * setters and mesh mutation must not overlap queries. The lazy build is
   * guarded by a mutex and
   * transformations are immutable after construction.
   */
  template <class MeshType>
  class AABB
  {
    public:
      /// @brief Builds a locator bound to a fixed mesh.
      /// @param mesh Mesh on which the object is defined.
      explicit AABB(const MeshType& mesh)
        : m_mesh(mesh),
          m_tolerance(DefaultPhysicalTolerance),
          m_referenceTolerance(DefaultReferenceTolerance),
          m_maxNewtonIterations(DefaultMaxNewtonIterations),
          m_exhaustiveFallback(false),
          m_projectionPruning(false),
          m_index(mesh.getDimension() + 1)
      {
        computeScale();
      }

      /// Relative physical tolerance (scaled by the mesh diagonal).
      /// @returns The tolerance.
      Real getTolerance() const
      {
        return m_tolerance;
      }

      /// Sets the relative physical tolerance and invalidates the index.
      /// @param tolerance Tolerance used by the operation.
      /// @returns Reference to this object after the operation.
      AABB& setTolerance(Real tolerance)
      {
        if (!std::isfinite(tolerance) || tolerance < Real(0))
          throw std::invalid_argument(
            "AABB physical tolerance must be finite and nonnegative.");
        m_tolerance = tolerance;
        invalidate();
        return *this;
      }

      /// Maximum reference-space overshoot before clipping and residual check.
      /// @returns The reference tolerance.
      Real getReferenceTolerance() const
      {
        return m_referenceTolerance;
      }

      /// @brief Sets the maximum reference-space overshoot before clipping.
      /// @param tolerance Tolerance used by the operation.
      /// @returns Reference to this object after the operation.
      AABB& setReferenceTolerance(Real tolerance)
      {
        if (!std::isfinite(tolerance) || tolerance < Real(0))
          throw std::invalid_argument(
            "AABB reference tolerance must be finite and nonnegative.");
        m_referenceTolerance = tolerance;
        return *this;
      }

      /**
       * @brief Enables the exhaustive narrow-phase fallback on broad-phase
       * miss.
       *
       * Useful as a diagnostic when a transformation's degree metadata does
       * not describe its image. Unreliable control-point conversions already
       * leave their entries unpruned. Costs one full narrow-phase sweep per miss.
       * @returns Reference to this object after the operation.
       * @param fallback Whether to enable the exhaustive fallback search.
       */
      AABB& setExhaustiveFallback(bool fallback)
      {
        m_exhaustiveFallback = fallback;
        return *this;
      }

      /**
       * @brief Enables control-hull projections in addition to axis-aligned boxes.
       *
       * Disabled by default because construction and storage costs depend on
       * geometry and query reuse. Enable to reduce Newton retries on overlapping
       * boxes during repeated point location. Membership checks and
       * Newton seed retries are unchanged. Invalidates the existing index.
       * @returns Reference to this object after the operation.
       * @param enabled Whether to enable projection-based pruning.
       */
      AABB& setProjectionPruning(bool enabled)
      {
        m_projectionPruning = enabled;
        invalidate(false);
        return *this;
      }

      /**
       * @brief Locates a physical point on a polytope of the requested dimension.
       *
       * Returns an empty optional when the point is outside every candidate
       * polytope, when the coordinate dimension is incompatible with the mesh,
       * or when the inverse transformation does not pass the residual checks.
       * @param x Point at which the operation is evaluated.
       * @param dimension Topological dimension of the entities to search.
       * @returns Located mesh point, or an empty optional when the query cannot be certified.
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
      /// @param x Point at which the operation is evaluated.
      /// @returns Located mesh point, or an empty optional when the query cannot be certified.
      Optional<Geometry::Point> locate(const Math::SpatialPoint& x) const
      {
        return locate(m_mesh.get().getDimension(), x);
      }

    private:
      /// Ambient dimension supported by the spatial coordinate storage.
      static constexpr size_t MaxSpaceDimension = 3;
      /// Leaf capacity: a policy balance between tree traversal and candidate scans.
      static constexpr size_t LeafSize = 8;
      /// Median splits give logarithmic depth; this capacity accommodates the
      /// supported 32-bit entry range with spare traversal slots.
      static constexpr int32_t StackDepth = 64;

      /// Default physical residual tolerance, relative to the mesh diagonal.
      static constexpr Real DefaultPhysicalTolerance = Real(1e-10);
      /// Default permitted reference-space overshoot before clipping.
      static constexpr Real DefaultReferenceTolerance = Real(1e-10);
      /// Work limit per Newton seed; it does not guarantee convergence.
      static constexpr size_t DefaultMaxNewtonIterations = 16;
      /// Roundoff allowance in reference coordinates, measured in machine epsilons.
      /// This is a numerical policy margin, not an error bound for the inverse map.
      static constexpr Real ReferenceRoundoffFactor = Real(16);
      /// Divergence guard: reference cells have coordinates of order one.
      /// The generous radius permits intermediate overshoot without accepting it.
      static constexpr Real MaxReferenceNorm = Real(1e3);
      /// Squared divergence radius, derived for squared-length comparisons.
      static constexpr Real MaxReferenceNormSquared = MaxReferenceNorm * MaxReferenceNorm;
      /// Try a vertex and its midpoint with the centroid after the centroid seed.
      static constexpr size_t SeedsPerVertex = 2;
      /// Equal weights place the second seed exactly at that midpoint.
      static constexpr Real SeedCentroidWeight = Real(0.5);
      /// Halve a rejected step to search toward the current iterate.
      static constexpr Real BacktrackingContraction = Real(0.5);
      /// Work limit including the full-step trial; the last scale is 2^-23.
      /// Exhausting this budget rejects the seed, rather than certifying a miss.
      static constexpr size_t MaxBacktrackingTrials = 24;
      /// Require the estimated correction to fit within one quarter of the
      /// reference accuracy before skipping another Jacobian evaluation.
      /// This conservative policy margin is heuristic, not a convergence proof.
      static constexpr Real CorrectionEstimateMargin = Real(0.25);

      /// Conversion roundoff allowance in machine epsilons per basis entry.
      /// This heuristic accounts for conditioning; it is not an interval proof.
      static constexpr Real ControlRoundoffFactor = Real(32);
      /// Cubic and higher tensor factors use line conversions. Small dense
      /// products avoid the line-loop overhead on quadratic factors.
      static constexpr size_t MinSeparableTensorDegree = 3;
      /// Euclidean norm with the cheap squared-norm path at ordinary scales.
      /// Scaling avoids overflow and avoids classifying nonzero tiny vectors as zero.
      static Real stableNorm(const Math::SpatialPoint& v)
      {
        const Real squared = v.squaredNorm();
        if (std::isnormal(squared))
          return std::sqrt(squared);
        Real norm = 0;
        for (Eigen::Index i = 0; i < v.size(); ++i)
          norm = std::hypot(norm, v[i]);
        return norm;
      }

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

      struct ProjectionBound
      {
          Bound normal;
          Real upper;
      };

      struct ProjectionRange
      {
          size_t begin = 0;
          size_t end = 0;
      };

      struct DimensionIndex
      {
          std::vector<Node> nodes;
          std::vector<Index> entries;      ///< Polytope indices, leaf-ordered
          std::vector<Bound> entryLo;      ///< Per-entry box, leaf-ordered
          std::vector<Bound> entryHi;
          std::vector<ProjectionBound> projections;
          std::vector<ProjectionRange> projectionRanges;
          std::atomic<bool> built{false};
          mutable std::mutex mutex;
      };

      void invalidate(bool updateScale = true)
      {
        for (auto& index : m_index)
        {
          std::lock_guard lock(index.mutex);
          index.nodes.clear();
          index.entries.clear();
          index.entryLo.clear();
          index.entryHi.clear();
          index.projections.clear();
          index.projectionRanges.clear();
          index.built.store(false, std::memory_order_release);
        }
        if (updateScale)
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
        Real diag = 0;
        for (size_t i = 0; i < sdim; ++i)
        {
          const Real e = hi[i] - lo[i];
          if (std::isfinite(e))
            diag = std::hypot(diag, e);
          else if (lo[i] <= hi[i])
            diag = std::numeric_limits<Real>::max();
        }
        m_scale =
          diag > Real(0) ? std::min(diag, std::numeric_limits<Real>::max()) : Real(1);
      }

      /// Effective physical tolerance in mesh units.
      Real physicalTolerance() const
      {
        return std::min(m_tolerance * m_scale, std::numeric_limits<Real>::max());
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
        const bool projectionsEnabled =
          m_projectionPruning && dimension == sdim && dimension > 1;

        index.entries.clear();
        index.nodes.clear();
        index.entryLo.clear();
        index.entryHi.clear();
        index.projections.clear();
        index.projectionRanges.clear();
        if (count == 0)
          return;

        // Per-entry boxes and centroids in flat storage.
        std::vector<Bound> lo(count), hi(count);
        std::vector<Bound> mid(count);
        index.entries.reserve(count);
        std::vector<ProjectionRange> ranges(projectionsEnabled ? count : 0);
        size_t n = 0;
        for (auto it = mesh.getPolytope(dimension); it; ++it, ++n)
        {
          if (projectionsEnabled)
            ranges[n].begin = index.projections.size();
          makeBox(*it, lo[n], hi[n], index.projections);
          if (projectionsEnabled)
            ranges[n].end = index.projections.size();
          for (size_t i = 0; i < sdim; ++i)
            mid[n][i] = std::isfinite(lo[n][i]) && std::isfinite(hi[n][i])
              ? Real(0.5) * lo[n][i] + Real(0.5) * hi[n][i]
              : Real(0);
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
        index.projectionRanges.resize(projectionsEnabled ? count : 0);
        for (size_t i = 0; i < count; ++i)
        {
          reordered[i] = index.entries[order[i]];
          index.entryLo[i] = lo[order[i]];
          index.entryHi[i] = hi[order[i]];
          if (projectionsEnabled)
            index.projectionRanges[i] = ranges[order[i]];
        }
        index.entries = std::move(reordered);
        if (index.projections.empty())
          std::vector<ProjectionRange>().swap(index.projectionRanges);
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
            return bernstein(k, mode[0], r.x()) * bernstein(k, mode[1], r.y()) *
              bernstein(k, mode[2], r.z());
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
            return bernstein(n, mode[0], r.x() / q) * bernstein(n, mode[1], r.y() / q) *
              bz;
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
                         std::initializer_list<Real> coords) {
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
              emplace({i, 0, 0, 0}, {static_cast<Real>(i) / kr});
            return;
          }
          case G::Triangle:
          {
            for (size_t a2 = 0; a2 <= k; ++a2)
            {
              for (size_t a1 = 0; a1 + a2 <= k; ++a1)
              {
                emplace({k - a1 - a2, a1, a2, 0},
                  {static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr});
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
                  emplace({k - a1 - a2 - a3, a1, a2, a3},
                    {static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr,
                      static_cast<Real>(a3) / kr});
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
                emplace(
                  {i, j, 0, 0}, {static_cast<Real>(i) / kr, static_cast<Real>(j) / kr});
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
                  emplace({i, j, l, 0},
                    {static_cast<Real>(i) / kr, static_cast<Real>(j) / kr,
                      static_cast<Real>(l) / kr});
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
                  emplace({k - a1 - a2, a1, a2, m},
                    {static_cast<Real>(a1) / kr, static_cast<Real>(a2) / kr,
                      static_cast<Real>(m) / kr});
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
                emplace({0, 0, m, 0}, {Real(0), Real(0), z});
                continue;
              }
              const Real nr = static_cast<Real>(n);
              for (size_t j = 0; j <= n; ++j)
              {
                for (size_t i = 0; i <= n; ++i)
                {
                  emplace({i, j, m, 0},
                    {q * static_cast<Real>(i) / nr, q * static_cast<Real>(j) / nr, z});
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
       * image in exact arithmetic. Floating-point conversion requires a
       * roundoff allowance; ill-conditioned conversions cannot certify a box.
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
          Real roundoffAmplification = 0;
          size_t tensorFactors = 0;
          bool valid = false;
      };

      /// @brief Returns the cached control-point basis of degree @p k on @p g.
      static const ControlBasis& getControlBasis(Geometry::Polytope::Type g, size_t k)
      {
        using Key = std::pair<int, size_t>;
        static std::map<Key, ControlBasis> s_cache;
        static std::mutex s_mutex;

        const Key key{static_cast<int>(g), k};
        std::lock_guard<std::mutex> lock(s_mutex);
        const auto it = s_cache.find(key);
        if (it != s_cache.end())
          return it->second;

        ControlBasis basis;
        std::vector<std::array<size_t, 4>> modes;
        makeControlLattice(g, k, modes, basis.samples);
        const size_t n = modes.size();
        using G = Geometry::Polytope::Type;
        basis.tensorFactors = k >= MinSeparableTensorDegree
          ? (g == G::Quadrilateral ? 2 : (g == G::Hexahedron ? 3 : 0))
          : 0;
        if (basis.tensorFactors > 0)
        {
          // Tensor conversion is a succession of small one-dimensional solves,
          // not an inverse of the full (k+1)^d by (k+1)^d lattice matrix.
          Math::Matrix<Real> vandermonde(k + 1, k + 1);
          for (size_t m = 0; m <= k; ++m)
            for (size_t c = 0; c <= k; ++c)
              vandermonde(m, c) = bernstein(k, c, static_cast<Real>(m) / k);
          const Eigen::FullPivLU<Math::Matrix<Real>> lu(vandermonde);
          if (lu.isInvertible())
          {
            basis.conversion = lu.inverse().transpose();
            basis.roundoffAmplification = ControlRoundoffFactor *
              static_cast<Real>((k + 1) * basis.tensorFactors) *
              std::numeric_limits<Real>::epsilon() /
              std::pow(lu.rcond(), static_cast<Real>(basis.tensorFactors));
            basis.valid = basis.conversion.allFinite() &&
              std::isfinite(basis.roundoffAmplification) &&
              basis.roundoffAmplification < Real(1);
          }
          return s_cache.emplace(key, std::move(basis)).first->second;
        }
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
            basis.roundoffAmplification = ControlRoundoffFactor * static_cast<Real>(n) *
              std::numeric_limits<Real>::epsilon() / lu.rcond();
            basis.valid = basis.conversion.allFinite() &&
              std::isfinite(basis.roundoffAmplification) &&
              basis.roundoffAmplification < Real(1);
          }
        }
        return s_cache.emplace(key, std::move(basis)).first->second;
      }

      void padBox(Bound& lo, Bound& hi, Real pad, size_t sdim) const
      {
        for (size_t i = 0; i < sdim; ++i)
        {
          lo[i] = std::nextafter(lo[i] - pad, -std::numeric_limits<Real>::infinity());
          hi[i] = std::nextafter(hi[i] + pad, std::numeric_limits<Real>::infinity());
        }
      }

      /**
       * @brief Adds conservative hull projections along mapped reference-face normals.
       *
       * Every image point is a convex combination of the control points, so its
       * projection cannot exceed their maximum. Axis-aligned directions are
       * already covered by the AABB. Normals select useful directions only;
       * they do not assert that curved faces themselves are planar.
       */
      template <class GetControlPoint>
      void makeProjections(const Geometry::Polytope& polytope, size_t count,
        GetControlPoint&& getControlPoint, Real conversionError,
        std::vector<ProjectionBound>& projections) const
      {
        const Geometry::Polytope::Traits traits(polytope.getGeometry());
        const size_t dimension = polytope.getDimension();
        const size_t sdim = m_mesh.get().getSpaceDimension();
        if (!m_projectionPruning || dimension != sdim || dimension < 2)
          return;
        Math::SpatialMatrix<Real> jac;
        polytope.getTransformation().jacobian(jac, traits.getCentroid());
        if (!isFinite(jac) || !std::isnormal(jac.determinant()))
          return;
        const auto& hs = traits.getHalfSpace();
        for (Eigen::Index face = 0; face < hs.vector.size(); ++face)
        {
          Math::SpatialPoint refNormal(dimension);
          for (size_t i = 0; i < dimension; ++i)
            refNormal[i] = hs.matrix(face, i);
          auto normal = jac.transpose().solve(refNormal);
          const Real norm = stableNorm(normal);
          if (!(norm > Real(0)) || !std::isfinite(norm))
            continue;
          normal /= norm;
          size_t nonzero = 0;
          for (size_t i = 0; i < sdim; ++i)
            nonzero += normal[i] != Real(0);
          if (nonzero <= 1)
            continue;
          Real upper = -std::numeric_limits<Real>::infinity();
          Real coordinateScale = 0;
          Math::SpatialPoint control;
          for (size_t j = 0; j < count; ++j)
          {
            getControlPoint(control, j);
            upper = std::max(upper, control.dot(normal));
            for (size_t i = 0; i < sdim; ++i)
              coordinateScale = std::max(coordinateScale, std::abs(control[i]));
          }
          const Real error = std::sqrt(static_cast<Real>(sdim)) * conversionError +
            ControlRoundoffFactor * std::numeric_limits<Real>::epsilon() *
              coordinateScale;
          upper = std::nextafter(
            upper + physicalTolerance() + error, std::numeric_limits<Real>::infinity());
          if (!std::isfinite(upper))
            continue;
          ProjectionBound bound{};
          for (size_t i = 0; i < sdim; ++i)
            bound.normal[i] = normal[i];
          bound.upper = upper;
          projections.push_back(bound);
        }
      }

      void makeBox(const Geometry::Polytope& polytope, Bound& lo, Bound& hi,
        std::vector<ProjectionBound>& projections) const
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
            lo[i] = hi[i] = x[static_cast<Eigen::Index>(i)];
          }
          padBox(lo, hi, physTol, sdim);
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
        const size_t degree = transformation.getFactorOrder();
        if (degree <= 1)
        {
          std::array<Math::SpatialPoint, RODIN_MAXIMUM_POLYTOPE_VERTICES> vertices;
          for (size_t i = 0; i < nv; ++i)
          {
            transformation.transform(vertices[i], traits.getVertex(i));
            add(vertices[i]);
          }
          makeProjections(
            polytope, nv, [&](Math::SpatialPoint& out, size_t j) { out = vertices[j]; },
            Real(0), projections);
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
          // Center before conversion: translating a cell far from the origin
          // must not amplify cancellation in recovered control-point offsets.
          const Math::Vector<Real> origin = sampled.col(0);
          sampled.colwise() -= origin;
          Math::Matrix<Real> control;
          if (basis.tensorFactors == 0)
            control = sampled * basis.conversion;
          else
          {
            control = sampled;
            Math::Matrix<Real> next(sdim, n);
            const size_t width = degree + 1;
            size_t stride = 1;
            for (size_t axis = 0; axis < basis.tensorFactors; ++axis)
            {
              for (size_t block = 0; block < n; block += width * stride)
                for (size_t offset = 0; offset < stride; ++offset)
                  for (size_t j = 0; j < width; ++j)
                  {
                    const size_t dst = block + offset + j * stride;
                    next.col(dst).setZero();
                    for (size_t c = 0; c < width; ++c)
                      next.col(dst) +=
                        control.col(block + offset + c * stride) * basis.conversion(c, j);
                  }
              control.swap(next);
              stride *= width;
            }
          }
          for (size_t i = 0; i < sdim; ++i)
          {
            const auto row = control.row(static_cast<Eigen::Index>(i));
            const Real roundoff = basis.roundoffAmplification *
              sampled.row(static_cast<Eigen::Index>(i)).cwiseAbs().maxCoeff();
            lo[i] = origin[i] + row.minCoeff() - roundoff;
            hi[i] = origin[i] + row.maxCoeff() + roundoff;
          }
          if (control.allFinite())
          {
            const Real error =
              basis.roundoffAmplification * sampled.cwiseAbs().maxCoeff();
            makeProjections(
              polytope, n,
              [&](Math::SpatialPoint& out, size_t j) {
                out.resize(sdim);
                for (size_t i = 0; i < sdim; ++i)
                  out[i] = origin[i] + control(i, j);
              },
              error, projections);
            padBox(lo, hi, physTol, sdim);
            return;
          }
        }

        // An unreliable conversion must not silently turn a sampled box into
        // a membership filter. Leave this entry unpruned so every query tests it.
        lo.fill(-std::numeric_limits<Real>::infinity());
        hi.fill(std::numeric_limits<Real>::infinity());
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

      bool containsReference(const Geometry::Polytope::Traits& traits,
        const Math::SpatialPoint& rc, bool& needsClip) const
      {
        needsClip = false;
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
          needsClip |= margin < Real(0);
        }
        return true;
      }

      void clipReference(Geometry::Polytope::Type geometry, Math::SpatialPoint& rc) const
      {
        using G = Geometry::Polytope::Type;
        auto clipSimplex = [&](size_t dimension) {
          Real sum = 0;
          for (size_t i = 0; i < dimension; ++i)
          {
            rc[i] = std::max(Real(0), rc[i]);
            sum += rc[i];
          }
          if (sum > Real(1))
          {
            for (size_t i = 0; i < dimension; ++i)
              rc[i] /= sum;
          }
        };

        switch (geometry)
        {
          case G::Point:
            break;
          case G::Segment:
            rc[0] = std::clamp(rc[0], Real(0), Real(1));
            break;
          case G::Triangle:
            clipSimplex(2);
            break;
          case G::Tetrahedron:
            clipSimplex(3);
            break;
          case G::Quadrilateral:
            rc[0] = std::clamp(rc[0], Real(0), Real(1));
            rc[1] = std::clamp(rc[1], Real(0), Real(1));
            break;
          case G::Hexahedron:
            for (size_t i = 0; i < 3; ++i)
              rc[i] = std::clamp(rc[i], Real(0), Real(1));
            break;
          case G::Wedge:
            clipSimplex(2);
            rc[2] = std::clamp(rc[2], Real(0), Real(1));
            break;
          case G::Pyramid:
            rc[2] = std::clamp(rc[2], Real(0), Real(1));
            rc[0] = std::clamp(rc[0], Real(0), Real(1) - rc[2]);
            rc[1] = std::clamp(rc[1], Real(0), Real(1) - rc[2]);
            break;
        }
      }

      /// @brief Inverts a valid, injective map from several seeds with controlled steps.
      bool invert(const Geometry::PolytopeTransformation& transformation,
        Geometry::Polytope::Type geometry, const Geometry::Polytope::Traits& traits,
        const Math::SpatialPoint& x, Math::SpatialPoint& rc) const
      {
        const Real physTol = physicalTolerance();

        const Real referenceAccuracy = std::max(
          m_tolerance, ReferenceRoundoffFactor * std::numeric_limits<Real>::epsilon());

        const size_t rdim = traits.getDimension();
        const size_t pdim = static_cast<size_t>(x.size());
        using G = Geometry::Polytope::Type;
        const bool affineSimplex = transformation.getOrder() <= 1 &&
          (geometry == G::Segment || geometry == G::Triangle ||
            geometry == G::Tetrahedron);
        Math::SpatialPoint mapped;
        Math::SpatialMatrix<Real> jac;
        Math::SpatialPoint residual;
        Math::SpatialPoint step;
        Math::SpatialPoint candidate;
        Math::SpatialPoint candidateMapped;
        auto accept = [&]() {
          bool needsClip;
          if (!containsReference(traits, rc, needsClip))
            return false;
          if (!needsClip)
            return true;
          clipReference(geometry, rc);
          Math::SpatialPoint clippedMapped;
          transformation.transform(clippedMapped, rc);
          return isFinite(clippedMapped) && stableNorm(x - clippedMapped) <= physTol;
        };

        for (size_t seed = 0;
             seed == 0 || seed < 1 + SeedsPerVertex * traits.getVertexCount(); ++seed)
        {
          if (seed == 0)
            rc = traits.getCentroid();
          else if (seed % SeedsPerVertex == 1)
            rc = traits.getVertex((seed - 1) / SeedsPerVertex);
          else
            rc = SeedCentroidWeight *
              (traits.getCentroid() +
                traits.getVertex((seed - SeedsPerVertex) / SeedsPerVertex));
          transformation.transform(mapped, rc);
          if (!isFinite(mapped))
            continue;
          residual = x - mapped;
          Real residualNorm = stableNorm(residual);
          for (size_t iteration = 0; iteration < m_maxNewtonIterations; ++iteration)
          {
            if (residualNorm == Real(0))
            {
              if (accept())
                return true;
              // Injectivity inside the reference cell says nothing about the
              // polynomial extension outside it. Try the next seed.
              break;
            }
            transformation.jacobian(jac, rc);
            if (!isFinite(jac))
              break;

            if (rdim == pdim)
            {
              Real determinant = jac.determinant();
              if (!std::isnormal(determinant))
              {
                // Rescale the system before the determinant test: small valid
                // cells can have determinants that underflow, and large ones
                // can overflow. The solution is unchanged by common scaling.
                Real scale = 0;
                for (size_t i = 0; i < pdim; ++i)
                  for (size_t j = 0; j < rdim; ++j)
                    scale = std::max(scale, std::abs(jac(i, j)));
                if (!(scale > Real(0)))
                  break;
                for (size_t i = 0; i < pdim; ++i)
                {
                  residual[i] /= scale;
                  for (size_t j = 0; j < rdim; ++j)
                    jac(i, j) /= scale;
                }
                determinant = jac.determinant();
                if (!std::isfinite(determinant) || determinant == Real(0))
                  break;
              }
              if (rdim == 1)
              {
                step.resize(1);
                step[0] = residual[0] / determinant;
              }
              else
              {
                step = jac.solve(residual);
              }
            }
            else
            {
              // Bounded storage keeps these <=3-dimensional QR solves off the heap.
              using Dense = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::ColMajor, MaxSpaceDimension, MaxSpaceDimension>;
              using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1, Eigen::ColMajor,
                MaxSpaceDimension, 1>;
              Dense dense(pdim, rdim);
              Vector rhs(pdim);
              for (size_t i = 0; i < pdim; ++i)
              {
                rhs[i] = residual[i];
                for (size_t j = 0; j < rdim; ++j)
                  dense(i, j) = jac(i, j);
              }
              const Real scale = dense.cwiseAbs().maxCoeff();
              if (scale > Real(0) && !std::isnormal(scale * scale))
              {
                dense /= scale;
                rhs /= scale;
              }
              Eigen::ColPivHouseholderQR<Dense> qr(dense);
              if (qr.rank() < static_cast<Eigen::Index>(rdim))
                break;
              step = qr.solve(rhs);
            }
            if (!isFinite(step))
              break;
            if (affineSimplex)
            {
              rc += step;
              if (rc.squaredNorm() > MaxReferenceNormSquared)
                return false;
              bool needsClip;
              if (!containsReference(traits, rc, needsClip))
                return false;
              if (needsClip)
                clipReference(geometry, rc);
              transformation.transform(mapped, rc);
              return isFinite(mapped) && stableNorm(x - mapped) <= physTol;
            }
            const Real stepNorm = stableNorm(step);
            // A small physical residual alone can hide a large coordinate
            // error on a thin element. Also require a small reference update.
            if (residualNorm <= physTol && stepNorm <= referenceAccuracy)
            {
              if (accept())
                return true;
              break;
            }
            if (stepNorm == Real(0))
              break;

            Real alpha = 1;
            bool advanced = false;
            bool retrySeed = false;
            candidate = rc + step;
            for (size_t trial = 0; trial < MaxBacktrackingTrials; ++trial)
            {
              if (candidate.squaredNorm() <= MaxReferenceNormSquared)
              {
                transformation.transform(candidateMapped, candidate);
                residual = x - candidateMapped;
                const Real candidateResidualNorm = stableNorm(residual);
                // Strict decrease also rejects NaN and infinite trial residuals.
                if (candidateResidualNorm < residualNorm)
                {
                  rc = candidate;
                  if (candidateResidualNorm == Real(0))
                  {
                    if (accept())
                      return true;
                    retrySeed = true;
                    break;
                  }
                  // Estimate the remaining reference correction using
                  // this step's residual ratio. Skip the next Jacobian
                  // only when the estimate is well below tolerance.
                  if (candidateResidualNorm <= physTol && alpha == Real(1) &&
                    candidateResidualNorm / residualNorm * stepNorm <=
                      CorrectionEstimateMargin * referenceAccuracy)
                  {
                    if (accept())
                      return true;
                    retrySeed = true;
                    break;
                  }
                  residualNorm = candidateResidualNorm;
                  advanced = true;
                  break;
                }
              }
              alpha *= BacktrackingContraction;
              candidate = rc + alpha * step;
            }
            if (retrySeed)
              break;
            if (!advanced)
            {
              // A floating-point stationary point can have a small residual
              // even when another Newton correction cannot reduce it.
              if (residualNorm <= physTol && accept())
                return true;
              if (residualNorm <= physTol)
                break;
              break;
            }
          }
          if (affineSimplex)
            return false;
        }
        return false;
      }

      Optional<Geometry::Point> narrowPhase(
        size_t dimension, Index polytopeIndex, const Math::SpatialPoint& x) const
      {
        const auto& mesh = m_mesh.get();
        const Geometry::Polytope polytope = *mesh.getPolytope(dimension, polytopeIndex);

        if (dimension == 0)
        {
          if (stableNorm(mesh.getVertexCoordinates(polytopeIndex) - x) <=
            physicalTolerance())
            return Geometry::Point(polytope, Math::SpatialPoint(0), x);
          return {};
        }

        const auto geometry = polytope.getGeometry();
        const Geometry::Polytope::Traits traits(geometry);
        Math::SpatialPoint rc;
        if (!invert(polytope.getTransformation(), geometry, traits, x, rc))
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
              bool outsideHull = false;
              const auto range = index.projectionRanges.empty()
                ? ProjectionRange{}
                : index.projectionRanges[k];
              for (size_t j = range.begin; j < range.end; ++j)
              {
                const auto& projection = index.projections[j];
                Real value = 0;
                for (size_t i = 0; i < sdim; ++i)
                  value += projection.normal[i] * x[i];
                if (value > projection.upper)
                {
                  outsideHull = true;
                  break;
                }
              }
              if (outsideHull)
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
      bool m_projectionPruning;
      Real m_scale;
      mutable std::vector<DimensionIndex> m_index;
  };
}

#endif
