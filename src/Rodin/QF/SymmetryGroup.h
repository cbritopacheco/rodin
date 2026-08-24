/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Defines SymmetryGroup, the symmetries of a reference element and the
 * orbit types they induce.
 */
#ifndef RODIN_QF_SYMMETRYGROUP_H
#define RODIN_QF_SYMMETRYGROUP_H

#include <algorithm>
#include <cmath>
#include <vector>
#include <functional>
#include <random>
#include <map>

#include <Eigen/Dense>

#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/Vector.h"
#include "Rodin/Types.h"

namespace Rodin::QF
{
  /**
   * @brief The affine symmetries of a reference element, and the kinds of
   * orbit they generate.
   *
   * A fully symmetric quadrature rule is a union of orbits of the element's
   * symmetry group, so a generator needs to know that group and the kinds of
   * orbit it has. Both are *derived* here rather than tabulated.
   *
   * The group is found by trying affine maps that send vertices to vertices
   * and keeping those that map the element onto itself. The orbit types follow
   * from it: a point lying on the mirror planes of some subgroup has a smaller
   * orbit than a generic one, and the possible orbit types are exactly the
   * subspaces such points can occupy --- the fixed subspaces of the subgroups.
   * Those are obtained by intersecting the fixed subspaces of individual
   * symmetries until the collection closes.
   *
   * Deriving this matters because tabulating it does not generalise. The
   * barycentric multiset patterns that describe simplex orbits are wrong for a
   * pyramid, whose group permutes the four base vertices and fixes the apex,
   * and a generator written around them needs a separate hand-written path per
   * element. Here every element is handled by the same code, and the pyramid is
   * not a special case.
   */
  class SymmetryGroup
  {
    public:
      /// @brief An affine map @f$ x \mapsto Ax + t @f$ of the element onto
      /// itself.
      struct Map
      {
          Math::Matrix<Real> linear;
          Math::Vector<Real> translation;

          Math::Vector<Real> operator()(const Math::Vector<Real>& x) const
          {
            return linear * x + translation;
          }
      };

      /**
       * @brief A kind of orbit: an affine subspace its seed may occupy,
       * together with the size of the resulting orbit.
       *
       * A seed is @f$ p + B\theta @f$, so @p dimension free parameters
       * describe it.
       */
      struct Stratum
      {
          Math::Vector<Real> point;   ///< A point of the subspace.
          Math::Matrix<Real> basis;   ///< Columns spanning its directions.
          size_t orbitSize = 1;       ///< Distinct images under the group.

          size_t getDimension() const
          {
            return static_cast<size_t>(basis.cols());
          }
      };

      /// @brief The symmetry group of @p g.
      static const std::vector<Map>& maps(Geometry::Polytope::Type g)
      {
        static std::map<Geometry::Polytope::Type, std::vector<Map>> cache;
        const auto found = cache.find(g);
        if (found != cache.end())
          return found->second;
        return cache.emplace(g, build(g)).first->second;
      }

      /// @brief The orbit types of @p g, coarsest orbit first.
      static const std::vector<Stratum>& strata(Geometry::Polytope::Type g)
      {
        static std::map<Geometry::Polytope::Type, std::vector<Stratum>> cache;
        const auto found = cache.find(g);
        if (found != cache.end())
          return found->second;
        return cache.emplace(g, buildStrata(g)).first->second;
      }

      /// @brief The distinct images of @p x under the group.
      static std::vector<Math::Vector<Real>> orbit(
        Geometry::Polytope::Type g, const Math::Vector<Real>& x, Real tolerance = 1e-10)
      {
        std::vector<Math::Vector<Real>> out;
        for (const auto& map : maps(g))
        {
          const Math::Vector<Real> y = map(x);
          bool seen = false;
          for (const auto& z : out)
            seen = seen || ((y - z).cwiseAbs().maxCoeff() < tolerance);
          if (!seen)
            out.push_back(y);
        }
        return out;
      }

    private:
      /// @brief Vertices of the reference element.
      static std::vector<Math::Vector<Real>> vertices(Geometry::Polytope::Type g)
      {
        const Geometry::Polytope::Traits traits(g);
        const size_t d = traits.getDimension();
        std::vector<Math::Vector<Real>> out;
        for (size_t i = 0; i < traits.getVertexCount(); ++i)
        {
          const auto& v = traits.getVertex(i);
          Math::Vector<Real> x(static_cast<Eigen::Index>(d));
          for (size_t k = 0; k < d; ++k)
            x(static_cast<Eigen::Index>(k)) = v[static_cast<Eigen::Index>(k)];
          out.push_back(std::move(x));
        }
        return out;
      }

      /**
       * @brief Finds every affine map of the element onto itself.
       *
       * An affine map is fixed by where it sends @f$ d + 1 @f$ affinely
       * independent vertices, so the search runs over the ordered tuples of
       * that many vertices and keeps the maps that carry the whole vertex set
       * onto itself. That is the symmetry group of the element, since an
       * affine bijection preserving the vertices preserves their convex hull.
       */
      static std::vector<Map> build(Geometry::Polytope::Type g)
      {
        constexpr Real tolerance = 1e-10;
        const auto vertex = vertices(g);
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        const size_t n = vertex.size();
        std::vector<Map> out;
        if (d == 0 || n == 0)
          return out;

        // An affinely independent frame to send elsewhere.
        std::vector<size_t> frame;
        Math::Matrix<Real> edges(
          static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d));
        for (size_t i = 0; i < n && frame.size() < d + 1; ++i)
        {
          frame.push_back(i);
          if (frame.size() < 2)
            continue;
          Math::Matrix<Real> trial(
            static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(frame.size() - 1));
          for (size_t k = 1; k < frame.size(); ++k)
            trial.col(static_cast<Eigen::Index>(k - 1)) =
              vertex[frame[k]] - vertex[frame[0]];
          Eigen::FullPivLU<Math::Matrix<Real>> lu(trial);
          if (static_cast<size_t>(lu.rank()) < frame.size() - 1)
            frame.pop_back();
        }
        if (frame.size() != d + 1)
          return out;
        for (size_t k = 1; k <= d; ++k)
          edges.col(static_cast<Eigen::Index>(k - 1)) =
            vertex[frame[k]] - vertex[frame[0]];
        const Math::Matrix<Real> inverse = edges.inverse();

        // Every ordered choice of images for the frame.
        std::vector<size_t> image(d + 1);
        std::vector<bool> used(n, false);
        const std::function<void(size_t)> rec = [&](size_t at) {
          if (at == d + 1)
          {
            Math::Matrix<Real> target(
              static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d));
            for (size_t k = 1; k <= d; ++k)
              target.col(static_cast<Eigen::Index>(k - 1)) =
                vertex[image[k]] - vertex[image[0]];
            Map candidate;
            candidate.linear = target * inverse;
            candidate.translation =
              vertex[image[0]] - candidate.linear * vertex[frame[0]];

            // Volume preserving, and carrying the vertex set onto itself.
            if (std::abs(std::abs(candidate.linear.determinant()) - 1) > tolerance)
              return;
            for (const auto& v : vertex)
            {
              const Math::Vector<Real> w = candidate(v);
              bool hit = false;
              for (const auto& u : vertex)
                hit = hit || ((w - u).cwiseAbs().maxCoeff() < tolerance);
              if (!hit)
                return;
            }
            out.push_back(std::move(candidate));
            return;
          }
          for (size_t i = 0; i < n; ++i)
          {
            if (used[i])
              continue;
            used[i] = true;
            image[at] = i;
            rec(at + 1);
            used[i] = false;
          }
        };
        rec(0);
        return out;
      }

      /// @brief The affine subspace fixed by @p map, if any.
      static bool fixedSpace(const Map& map, Stratum& out)
      {
        constexpr Real tolerance = 1e-9;
        const Eigen::Index d = map.linear.rows();
        const Math::Matrix<Real> shifted =
          map.linear - Math::Matrix<Real>::Identity(d, d);
        // Points with (A - I) x = -t, and directions in the null space.
        Eigen::JacobiSVD<Math::Matrix<Real>> svd(
          shifted, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Math::Vector<Real> particular = svd.solve(-map.translation);
        if ((shifted * particular + map.translation).cwiseAbs().maxCoeff() > tolerance)
          return false;   // no fixed point at all

        const auto& singular = svd.singularValues();
        std::vector<Eigen::Index> null;
        for (Eigen::Index i = 0; i < svd.matrixV().cols(); ++i)
          if (i >= singular.size() || singular(i) < tolerance)
            null.push_back(i);
        out.point = particular;
        out.basis.resize(d, static_cast<Eigen::Index>(null.size()));
        for (size_t k = 0; k < null.size(); ++k)
          out.basis.col(static_cast<Eigen::Index>(k)) = svd.matrixV().col(null[k]);
        return true;
      }

      /// @brief Intersection of two affine subspaces, if it is non-empty.
      static bool intersect(const Stratum& a, const Stratum& b, Stratum& out)
      {
        constexpr Real tolerance = 1e-9;
        const Eigen::Index d = a.point.size();
        // Solve a.point + a.basis s = b.point + b.basis t.
        Math::Matrix<Real> system(d, a.basis.cols() + b.basis.cols());
        system << a.basis, -b.basis;
        const Math::Vector<Real> rhs = b.point - a.point;
        Eigen::JacobiSVD<Math::Matrix<Real>> svd(
          system, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Math::Vector<Real> particular = svd.solve(rhs);
        if ((system * particular - rhs).cwiseAbs().maxCoeff() > tolerance)
          return false;

        out.point = a.point + a.basis * particular.head(a.basis.cols());
        // Directions common to both subspaces.
        const auto& singular = svd.singularValues();
        std::vector<Math::Vector<Real>> shared;
        for (Eigen::Index i = 0; i < svd.matrixV().cols(); ++i)
        {
          if (i < singular.size() && singular(i) >= tolerance)
            continue;
          const Math::Vector<Real> direction =
            a.basis * svd.matrixV().col(i).head(a.basis.cols());
          if (direction.norm() > tolerance)
            shared.push_back(direction.normalized());
        }
        // Orthonormalise what is left.
        Math::Matrix<Real> raw(d, static_cast<Eigen::Index>(shared.size()));
        for (size_t k = 0; k < shared.size(); ++k)
          raw.col(static_cast<Eigen::Index>(k)) = shared[k];
        if (raw.cols() == 0)
        {
          out.basis.resize(d, 0);
          return true;
        }
        Eigen::ColPivHouseholderQR<Math::Matrix<Real>> qr(raw);
        const Eigen::Index rank = qr.rank();
        const Math::Matrix<Real> q = qr.householderQ();
        out.basis = q.leftCols(rank);
        return true;
      }

      /// @brief An orthonormal basis for the column span of @p m.
      static Math::Matrix<Real> orthonormal(const Math::Matrix<Real>& m)
      {
        if (m.cols() == 0)
          return m;
        Eigen::ColPivHouseholderQR<Math::Matrix<Real>> qr(m);
        const Eigen::Index rank = qr.rank();
        const Math::Matrix<Real> q = qr.householderQ();
        return q.leftCols(rank);
      }

      /// @brief Whether two strata describe the same affine subspace.
      static bool same(const Stratum& a, const Stratum& b)
      {
        constexpr Real tolerance = 1e-8;
        if (a.basis.cols() != b.basis.cols())
          return false;
        // Same directions, and each point lying in the other subspace.
        const auto contains = [&](const Stratum& s, const Math::Vector<Real>& x) {
          const Math::Vector<Real> delta = x - s.point;
          if (s.basis.cols() == 0)
            return delta.cwiseAbs().maxCoeff() < tolerance;
          const Math::Vector<Real> residual =
            delta - s.basis * (s.basis.transpose() * delta);
          return residual.cwiseAbs().maxCoeff() < tolerance;
        };
        if (!contains(a, b.point) || !contains(b, a.point))
          return false;
        for (Eigen::Index k = 0; k < b.basis.cols(); ++k)
        {
          const Math::Vector<Real> direction = b.basis.col(k);
          const Math::Vector<Real> residual =
            direction - a.basis * (a.basis.transpose() * direction);
          if (residual.cwiseAbs().maxCoeff() > tolerance)
            return false;
        }
        return true;
      }

      /// @brief The orbit types, from the fixed subspaces of the subgroups.
      static std::vector<Stratum> buildStrata(Geometry::Polytope::Type g)
      {
        const auto& group = maps(g);
        const size_t d = Geometry::Polytope::Traits(g).getDimension();
        std::vector<Stratum> found;
        if (group.empty() || d == 0)
          return found;

        // The whole space, then what each symmetry fixes.
        Stratum whole;
        whole.point = Math::Vector<Real>::Zero(static_cast<Eigen::Index>(d));
        whole.basis = Math::Matrix<Real>::Identity(
          static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d));
        found.push_back(whole);
        for (const auto& map : group)
        {
          Stratum fixed;
          if (!fixedSpace(map, fixed))
            continue;
          bool seen = false;
          for (const auto& s : found)
            seen = seen || same(s, fixed);
          if (!seen)
            found.push_back(std::move(fixed));
        }

        // Close under intersection: a point may lie on several mirrors at once.
        for (size_t a = 0; a < found.size(); ++a)
        {
          for (size_t b = a + 1; b < found.size(); ++b)
          {
            Stratum meet;
            if (!intersect(found[a], found[b], meet))
              continue;
            bool seen = false;
            for (const auto& s : found)
              seen = seen || same(s, meet);
            if (!seen)
              found.push_back(std::move(meet));
          }
        }

        // The orbit size of each type, measured at a generic point of it.
        std::mt19937 rng(13579u);
        std::uniform_real_distribution<Real> uniform(0.11, 0.37);
        for (auto& stratum : found)
        {
          Math::Vector<Real> probe = stratum.point;
          for (Eigen::Index k = 0; k < stratum.basis.cols(); ++k)
            probe += uniform(rng) * stratum.basis.col(k);
          stratum.orbitSize = orbit(g, probe).size();
        }

        // Strata related by a symmetry describe the same kind of orbit --- the
        // three mirror lines of a triangle are one orbit type, not three --- so
        // only one representative of each class is kept. Enumerating the rest
        // would multiply the decompositions the search has to try without
        // offering it anything new.
        std::vector<Stratum> distinct;
        for (const auto& candidate : found)
        {
          bool equivalent = false;
          for (const auto& kept : distinct)
          {
            if (kept.orbitSize != candidate.orbitSize)
              continue;
            for (const auto& map : group)
            {
              Stratum moved;
              moved.point = map(candidate.point);
              // Re-orthonormalised: the symmetries of an element need not be
              // orthogonal maps --- the reference pyramid has its apex over a
              // corner, so its symmetries are shears --- and the subspace
              // comparison projects with this basis, which is only valid for
              // an orthonormal one.
              moved.basis = orthonormal(map.linear * candidate.basis);
              if (same(kept, moved))
              {
                equivalent = true;
                break;
              }
            }
            if (equivalent)
              break;
          }
          if (!equivalent)
            distinct.push_back(candidate);
        }
        found = std::move(distinct);

        std::stable_sort(found.begin(), found.end(),
          [](const Stratum& a, const Stratum& b) { return a.orbitSize < b.orbitSize; });
        return found;
      }
  };
}

#endif
