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
#include <functional>

#include "Rodin/Types.h"
#include "Rodin/Geometry.h"

#include "ForwardDecls.h"

namespace Rodin::Location
{
  /**
   * @brief Bounding-volume-hierarchy point locator over axis-aligned boxes.
   *
   * Broad phase: a packed median-split BVH over per-polytope AABBs, built
   * lazily per queried dimension. Narrow phase: Newton inversion of the
   * polytope transformation (exact in one step for affine maps, iterative
   * for bilinear/curved maps) followed by the reference half-space test and
   * a physical residual check, so off-manifold points are rejected and a
   * diverged inversion can never be accepted.
   *
   * Boxes of curved polytopes are inflated by the sampled chord deviation,
   * which conservatively bounds quadratic (P2) geometries; for higher-order
   * transformations the inflation is heuristic and an exhaustive narrow-phase
   * fallback can be enabled with setExhaustiveFallback().
   *
   * Tolerances are relative to the mesh bounding-box diagonal. Queries are
   * thread-safe: the lazy per-dimension build is guarded by a mutex and
   * transformations are immutable after construction.
   */
  template <class MeshType>
  class AABB
  {
    public:
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

      Optional<Geometry::Point> locate(
          size_t dimension,
          const Math::SpatialPoint& x) const
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

        // Reorder entries to leaf order so leaves are contiguous.
        std::vector<Index> reordered(count);
        for (size_t i = 0; i < count; ++i)
          reordered[i] = index.entries[order[i]];
        index.entries = std::move(reordered);
      }

      int32_t buildNode(
          DimensionIndex& index,
          std::vector<uint32_t>& order,
          uint32_t begin, uint32_t end,
          const std::vector<Bound>& lo,
          const std::vector<Bound>& hi,
          const std::vector<Bound>& mid,
          size_t sdim) const
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
          if (e > extent) { extent = e; axis = i; }
        }
        if (!(extent > Real(0)))
          return self;

        const uint32_t midIdx = (begin + end) / 2;
        std::nth_element(
            order.begin() + begin, order.begin() + midIdx, order.begin() + end,
            [&](uint32_t a, uint32_t b) { return mid[a][axis] < mid[b][axis]; });

        const int32_t left = buildNode(index, order, begin, midIdx, lo, hi, mid, sdim);
        const int32_t right = buildNode(index, order, midIdx, end, lo, hi, mid, sdim);
        index.nodes[self].left = left;
        index.nodes[self].right = right;
        return self;
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
        const Geometry::Polytope::Traits traits(polytope.getGeometry());
        const size_t nv = traits.getVertexCount();

        Math::SpatialPoint x;
        std::vector<Math::SpatialPoint> vertex(nv);
        auto add = [&](const Math::SpatialPoint& p)
        {
          for (size_t i = 0; i < sdim; ++i)
          {
            lo[i] = std::min(lo[i], p[static_cast<Eigen::Index>(i)]);
            hi[i] = std::max(hi[i], p[static_cast<Eigen::Index>(i)]);
          }
        };

        for (size_t i = 0; i < nv; ++i)
        {
          transformation.transform(vertex[i], traits.getVertex(i));
          add(vertex[i]);
        }
        transformation.transform(x, traits.getCentroid());
        add(x);

        // Chord-deviation inflation: for each vertex pair, the deviation of
        // the mapped reference midpoint from the straight chord midpoint. For
        // quadratic geometry the image of an edge is a parabola contained in
        // the convex hull of its Bezier control points, and inflating by the
        // sampled deviation bounds that hull. Affine maps contribute zero.
        Real deviation = 0;
        for (size_t i = 0; i < nv; ++i)
        {
          for (size_t j = i + 1; j < nv; ++j)
          {
            transformation.transform(
                x, Real(0.5) * (traits.getVertex(i) + traits.getVertex(j)));
            add(x);
            deviation = std::max(deviation,
                (x - Real(0.5) * (vertex[i] + vertex[j])).norm());
          }
        }

        const Real pad = physTol + deviation;
        for (size_t i = 0; i < sdim; ++i)
        {
          lo[i] -= pad;
          hi[i] += pad;
        }
      }

      bool boxContains(const Node& node, const Math::SpatialPoint& x, size_t sdim) const
      {
        for (size_t i = 0; i < sdim; ++i)
        {
          const Real xi = x[static_cast<Eigen::Index>(i)];
          if (xi < node.lo[i] || xi > node.hi[i])
            return false;
        }
        return true;
      }

      bool containsReference(
          const Geometry::Polytope::Traits& traits,
          const Math::SpatialPoint& rc) const
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
      bool invert(
          const Geometry::PolytopeTransformation& transformation,
          const Geometry::Polytope::Traits& traits,
          const Math::SpatialPoint& x,
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
        const Geometry::Polytope polytope =
          *mesh.getPolytope(dimension, polytopeIndex);

        if (dimension == 0)
        {
          if ((mesh.getVertexCoordinates(polytopeIndex) - x).norm()
              <= physicalTolerance())
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
          const DimensionIndex& index,
          size_t dimension,
          const Math::SpatialPoint& x) const
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
          const DimensionIndex& index,
          size_t dimension,
          const Math::SpatialPoint& x) const
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
