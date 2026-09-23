/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_THICKNESS_H
#define KELVIN_BALL_THICKNESS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Location.h>

#include "Common.h"
#include "SewedOutput.h"

namespace KelvinBall
{
  /**
   * @brief Minimum-thickness penalty of Allaire, Jouve and Michailidis.
   *
   * With \f$ n \f$ the normal from the body into the fluid,
   * @f[
   *   P(\Omega_s) = \int_\Gamma \int_0^{d_{\min}}
   *     \bigl[(d(s - \xi n(s)))_+\bigr]^2\,d\xi\,ds ,
   * @f]
   * where \f$ d \f$ is the signed distance to the complete interface, positive
   * in the fluid: an inward ray of length \f$ d_{\min} \f$ that leaves the body
   * reports how far it has left it. Its derivative in the normal velocity
   * \f$ w \f$ is
   * @f[
   *   dP(w) = \int_\Gamma \int_0^{d_{\min}}
   *     \bigl[2 d_+(x_m)\bigl(\nabla d(x_m) \cdot n(s)\,w(s) - w(y_m)\bigr)
   *       + H(s)d_+(x_m)^2w(s)\bigr]\,d\xi\,ds ,
   * @f]
   * with \f$ x_m = s - \xi n(s) \f$ and \f$ y_m \f$ its nearest point on the
   * interface. Here \f$ H=\operatorname{div}_{\Gamma}n \f$. The discrete normal used by
   * the rays is an H1 projection of the oriented interface normal, normalised
   * pointwise; its curvature is obtained by differentiating that normalisation.
   *
   * The chamber carries only one copy of the interface, and a thin part may lie
   * across a cut, so the distance is measured to the 24 rotated copies of the
   * chamber interface, found on a uniform grid of cell size \f$ d_{\min} \f$:
   * the nearest point of a ray point lies within the ray length, since the ray
   * starts on the interface. Whether a ray point lies in the fluid is read from
   * the label of the chamber cell containing its rotated image.
   */
  class ThicknessPenalty
  {
    public:
      struct Result
      {
          Real penalty = 0;
          Real deepest = 0;
          size_t rays = 0;
          size_t violating = 0;
          Real minimumCurvature = std::numeric_limits<Real>::infinity();
          Real maximumCurvature = -std::numeric_limits<Real>::infinity();
      };

      /**
       * @brief Extends and smooths the body normal in a vector H1 metric.
       *
       * The surface mass term anchors the extension to the piecewise-flat
       * interface normal. Rotational Nitsche matching is applied before the
       * solution is used on the chamber cuts. The returned field is unnormalised;
       * its point values are normalised after interpolation.
       */
      template <class Space, class Coupling>
      auto projectNormal(const Space& space, const Coupling& coupling,
        Real length, Real nitschePenalty) const
      {
        TrialFunction projected(space);
        TestFunction test(space);
        auto fluidNormal = FaceNormal(space.getMesh());
        fluidNormal.traceOf(Fluid);
        Problem projection(projected, test);
        projection = Integral(length * length * Jacobian(projected), Jacobian(test)) +
          Integral(projected, test) +
          length * FaceIntegral(projected, test).over(Gamma) +
          length * FaceIntegral(fluidNormal, test).over(Gamma);
        projection.assemble();
        coupling.assembleVector(space, projection.getLinearSystem(), length * length,
          nitschePenalty, FlatSet<Attribute>{});
        solveDirect(projection);
        const auto& system = projection.getLinearSystem();
        const Real residual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(residual) || residual > LinearResidualTolerance)
          throw std::runtime_error("The thickness normal projection did not converge.");
        return projected.getSolution();
      }

      template <class ChamberMesh>
      ThicknessPenalty(const ChamberMesh& mesh, Real minimum)
        : m_minimum(minimum),
          m_cell(minimum)
      {
        const auto& rotations = SewedOutput::getCubeRotations();
        const size_t D = mesh.getDimension();
        for (auto face = mesh.getPolytope(D - 1); face; ++face)
        {
          if (face->getAttribute() != Gamma)
            continue;
          const auto& vertices = face->getVertices();
          for (size_t r = 0; r < rotations.size(); ++r)
          {
            Triangle triangle;
            for (size_t k = 0; k < 3; ++k)
              triangle.x[k] = rotations[r] * mesh.getVertexCoordinates(vertices[k]);
            triangle.face = face->getIndex();
            triangle.rotation = r;
            const Index id = static_cast<Index>(m_triangles.size());
            m_triangles.push_back(triangle);
            const auto [lo, hi] = bounds(triangle);
            for (long i = lo[0]; i <= hi[0]; ++i)
              for (long j = lo[1]; j <= hi[1]; ++j)
                for (long k = lo[2]; k <= hi[2]; ++k)
                  m_grid[key(i, j, k)].push_back(id);
          }
        }
      }

      /**
       * @brief Evaluates the penalty over the chamber interface and adds
       * @f$ -\beta\,dP(\phi_i) @f$ to @p load for every basis field
       * @f$ \phi_i @f$ of @p space, whose degrees of freedom are nodal.
       */
      template <class ChamberMesh, class Locator, class Space, class Normal, class Load>
      Result evaluate(const ChamberMesh& mesh, const Locator& locator, const Space& space,
        const Normal& projectedNormal, Real weight, Load& load) const
      {
        static constexpr std::array<Real, 8> abscissae{0.0198550717512319,
          0.1016667612931866, 0.2372337950418355, 0.4082826787521751, 0.5917173212478249,
          0.7627662049581645, 0.8983332387068134, 0.9801449282487681};
        static constexpr std::array<Real, 8> weights{0.0506142681451881,
          0.1111905172266872, 0.1568533229389436, 0.1813418916891810, 0.1813418916891810,
          0.1568533229389436, 0.1111905172266872, 0.0506142681451881};
        static constexpr std::array<std::array<Real, 3>, 3> quadrature{
          {{Real(2) / 3, Real(1) / 6, Real(1) / 6},
            {Real(1) / 6, Real(2) / 3, Real(1) / 6},
            {Real(1) / 6, Real(1) / 6, Real(2) / 3}}};
        const auto& rotations = SewedOutput::getCubeRotations();
        const Real multiplicity = static_cast<Real>(rotations.size());
        const size_t D = mesh.getDimension();
        auto projectedJacobian = Jacobian(projectedNormal);
        projectedJacobian.traceOf(Fluid);
        const auto deposit = [&](Index face, const std::array<Real, 3>& barycentric,
                               const Math::SpatialVector<Real>& normal, Real value) {
          const auto& vertices = mesh.getPolytope(D - 1, face)->getVertices();
          for (size_t k = 0; k < 3; ++k)
          {
            const auto dofs = space.getDOFs(0, vertices[k]);
            for (Eigen::Index component = 0; component < 3; ++component)
              load(dofs(component)) += value * barycentric[k] * normal(component);
          }
        };

        Result result;
        for (auto face = mesh.getPolytope(D - 1); face; ++face)
        {
          if (face->getAttribute() != Gamma)
            continue;
          const auto& vertices = face->getVertices();
          const Math::SpatialVector<Real> a = mesh.getVertexCoordinates(vertices[0]);
          const Math::SpatialVector<Real> b = mesh.getVertexCoordinates(vertices[1]);
          const Math::SpatialVector<Real> c = mesh.getVertexCoordinates(vertices[2]);
          const Real area = cross(b - a, c - a).norm() / 2;
          for (const auto& point : quadrature)
          {
            const Math::SpatialVector<Real> s =
              point[0] * a + point[1] * b + point[2] * c;
            const Geometry::Point surfacePoint(*face, s);
            const Math::SpatialVector<Real> rawNormal =
              projectedNormal.getValue(surfacePoint);
            const Real normalMagnitude = rawNormal.norm();
            if (!(std::isfinite(normalMagnitude) && normalMagnitude > Real(1e-12)))
              throw std::runtime_error("The projected thickness normal vanishes.");
            const Math::SpatialVector<Real> normal = rawNormal / normalMagnitude;
            const auto gradientNormal = projectedJacobian.getValue(surfacePoint);
            const Real curvature =
              (gradientNormal.trace() - normal.dot(gradientNormal * normal)) /
              normalMagnitude;
            result.minimumCurvature = std::min(result.minimumCurvature, curvature);
            result.maximumCurvature = std::max(result.maximumCurvature, curvature);
            for (size_t j = 0; j < abscissae.size(); ++j)
            {
              const Real xi = abscissae[j] * m_minimum;
              const Real measure = multiplicity * (area / 3) * weights[j] * m_minimum;
              const Math::SpatialVector<Real> ray = s - xi * normal;
              ++result.rays;
              if (!inFluid(mesh, locator, ray))
                continue;
              const auto [distance, nearest, triangle, barycentric] = closest(ray, xi);
              if (!(distance > 0))
                continue;
              ++result.violating;
              result.deepest = std::max(result.deepest, distance);
              result.penalty += measure * distance * distance;
              const Math::SpatialVector<Real> gradient = (ray - nearest) / distance;
              // Local part at s: -beta * [2 d (grad d . n) + H d^2] w(s).
              deposit(face->getIndex(), point, normal,
                -weight * measure *
                  (2 * distance * gradient.dot(normal) + curvature * distance * distance));
              // Non-local part at the nearest point, returned to the chamber:
              // w(y) = theta(y_c) . n_c, so the load acts on the chamber face.
              const Triangle& hit = m_triangles[triangle];
              const auto& hitVertices = mesh.getPolytope(D - 1, hit.face)->getVertices();
              const Math::SpatialVector<Real> hitPoint =
                barycentric[0] * mesh.getVertexCoordinates(hitVertices[0]) +
                barycentric[1] * mesh.getVertexCoordinates(hitVertices[1]) +
                barycentric[2] * mesh.getVertexCoordinates(hitVertices[2]);
              Math::SpatialVector<Real> hitNormal =
                projectedNormal.getValue(
                  Geometry::Point(*mesh.getPolytope(D - 1, hit.face), hitPoint));
              if (!(hitNormal.norm() > Real(1e-12)))
                throw std::runtime_error("The projected nearest-point normal vanishes.");
              hitNormal.normalize();
              deposit(hit.face, barycentric, hitNormal, weight * measure * 2 * distance);
            }
          }
        }
        return result;
      }

    private:
      struct Triangle
      {
          std::array<Math::SpatialVector<Real>, 3> x;
          Index face;
          size_t rotation;
      };

      static Math::SpatialVector<Real> cross(
        const Math::SpatialVector<Real>& u, const Math::SpatialVector<Real>& v)
      {
        Math::SpatialVector<Real> w(3);
        w(0) = u(1) * v(2) - u(2) * v(1);
        w(1) = u(2) * v(0) - u(0) * v(2);
        w(2) = u(0) * v(1) - u(1) * v(0);
        return w;
      }

      static long long key(long i, long j, long k)
      {
        return (static_cast<long long>(i) * 73856093) ^
          (static_cast<long long>(j) * 19349663) ^ (static_cast<long long>(k) * 83492791);
      }

      long index(Real x) const
      {
        return static_cast<long>(std::floor(x / m_cell));
      }

      std::pair<std::array<long, 3>, std::array<long, 3>> bounds(const Triangle& t) const
      {
        std::array<long, 3> lo, hi;
        for (size_t d = 0; d < 3; ++d)
        {
          const Real min = std::min({t.x[0](d), t.x[1](d), t.x[2](d)});
          const Real max = std::max({t.x[0](d), t.x[1](d), t.x[2](d)});
          lo[d] = index(min);
          hi[d] = index(max);
        }
        return {lo, hi};
      }

      /// Whether a point lies in the fluid, read from the chamber cell that
      /// contains its rotated image; a point outside the container is fluid.
      template <class ChamberMesh, class Locator>
      static bool inFluid(const ChamberMesh& mesh, const Locator& locator,
        const Math::SpatialVector<Real>& x)
      {
        constexpr Real tolerance = 1e-12;
        for (const auto& rotation : SewedOutput::getCubeRotations())
        {
          const Math::SpatialVector<Real> z = rotation.transpose() * x;
          if (!(z(0) + tolerance >= z(1) && z(1) + tolerance >= std::abs(z(2))))
            continue;
          const auto located = locator.locate(z);
          if (!located)
            return true;
          return located->getPolytope().getAttribute() == Fluid;
        }
        return true;
      }

      /// Closest point of a triangle (Ericson, Real-Time Collision Detection).
      static std::pair<Math::SpatialVector<Real>, std::array<Real, 3>> closestOnTriangle(
        const Math::SpatialVector<Real>& p, const Triangle& t)
      {
        const auto& a = t.x[0];
        const auto& b = t.x[1];
        const auto& c = t.x[2];
        const Math::SpatialVector<Real> ab = b - a, ac = c - a, ap = p - a;
        const Real d1 = ab.dot(ap), d2 = ac.dot(ap);
        if (d1 <= 0 && d2 <= 0)
          return {a, {1, 0, 0}};
        const Math::SpatialVector<Real> bp = p - b;
        const Real d3 = ab.dot(bp), d4 = ac.dot(bp);
        if (d3 >= 0 && d4 <= d3)
          return {b, {0, 1, 0}};
        const Real vc = d1 * d4 - d3 * d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0)
        {
          const Real v = d1 / (d1 - d3);
          return {a + v * ab, {1 - v, v, 0}};
        }
        const Math::SpatialVector<Real> cp = p - c;
        const Real d5 = ab.dot(cp), d6 = ac.dot(cp);
        if (d6 >= 0 && d5 <= d6)
          return {c, {0, 0, 1}};
        const Real vb = d5 * d2 - d1 * d6;
        if (vb <= 0 && d2 >= 0 && d6 <= 0)
        {
          const Real w = d2 / (d2 - d6);
          return {a + w * ac, {1 - w, 0, w}};
        }
        const Real va = d3 * d6 - d5 * d4;
        if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
        {
          const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
          return {b + w * (c - b), {0, 1 - w, w}};
        }
        const Real denominator = 1 / (va + vb + vc);
        const Real v = vb * denominator, w = vc * denominator;
        return {a + v * ab + w * ac, {1 - v - w, v, w}};
      }

      /// Nearest point of the complete interface within @p radius of @p x.
      std::tuple<Real, Math::SpatialVector<Real>, Index, std::array<Real, 3>> closest(
        const Math::SpatialVector<Real>& x, Real radius) const
      {
        Real best = std::numeric_limits<Real>::infinity();
        Math::SpatialVector<Real> nearest = x;
        Index triangle = 0;
        std::array<Real, 3> barycentric{0, 0, 0};
        const Real r = radius * (1 + 1e-9) + 1e-12;
        for (long i = index(x(0) - r); i <= index(x(0) + r); ++i)
          for (long j = index(x(1) - r); j <= index(x(1) + r); ++j)
            for (long k = index(x(2) - r); k <= index(x(2) + r); ++k)
            {
              const auto it = m_grid.find(key(i, j, k));
              if (it == m_grid.end())
                continue;
              for (const Index id : it->second)
              {
                const auto [point, weights] = closestOnTriangle(x, m_triangles[id]);
                const Real distance = (point - x).norm();
                if (distance < best)
                {
                  best = distance;
                  nearest = point;
                  triangle = id;
                  barycentric = weights;
                }
              }
            }
        return {std::isfinite(best) ? best : Real(0), nearest, triangle, barycentric};
      }

      Real m_minimum;
      Real m_cell;
      std::vector<Triangle> m_triangles;
      std::unordered_map<long long, std::vector<Index>> m_grid;
  };
}

#endif
