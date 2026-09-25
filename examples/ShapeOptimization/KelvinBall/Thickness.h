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
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <Rodin/Geometry.h>

#include "Common.h"
#include "SewedOutput.h"

namespace KelvinBall
{
  /**
   * @brief A normal-ray minimum-thickness penalty for the sewn body.
   *
   * A surface H1 projection smooths the body-to-fluid normal. One inward
   * ray from each interface quadrature point finds its first outward crossing
   * of the complete surface, at distance @f$ t @f$. The penalty integrates
   * @f$ (d_{min}-t)_+^2 @f$. A single AABB tree of the chamber interface
   * is queried with all 24 inverse-rotated ray segments.
   *
   * The discrete variation differentiates source-face area and the
   * ray-triangle intersection, holding the regularized direction fixed.
   * The surface divergence of the normalized normal is evaluated as the
   * corresponding continuum curvature diagnostic.
   */
  class ThicknessPenalty
  {
    public:
      struct Result
      {
          Real penalty = 0;
          Real minimumExit = std::numeric_limits<Real>::infinity();
          Real maximumDeficit = 0;
          Real minimumTransversality = Real(1);
          size_t samples = 0;
          size_t violating = 0;
          size_t correctedNormals = 0;
          Real minimumCurvature = std::numeric_limits<Real>::infinity();
          Real maximumCurvature = -std::numeric_limits<Real>::infinity();
          Real minimumNormalMagnitude = std::numeric_limits<Real>::infinity();
          Real maximumNormalMagnitude = 0;
          Real minimumNormalAlignment = Real(1);
          Real meanNormalAlignment = 0;
          std::vector<Real> curvature;
      };

      /**
       * @brief Smooths the oriented normal along the interface alone.
       *
       * Opposing sides of a thin solid are not coupled through a bulk
       * extension. The resulting values are copied to the parent vertices;
       * only the interface trace is used by the thickness functional.
       */
      template <class Space>
      auto projectNormal(const Space& space, Real length) const
      {
        const auto& mesh = space.getMesh();
        SubMesh<Context::Local>::Builder builder;
        builder.initialize(mesh);
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          if (face->getAttribute() == Gamma)
            builder.include(mesh.getDimension() - 1, face->getIndex());
        SubMesh<Context::Local> surface = builder.finalize();
        P1 surfaceSpace(surface, 3);
        TrialFunction projected(surfaceSpace);
        TestFunction test(surfaceSpace);
        auto fluidNormal = FaceNormal(mesh);
        fluidNormal.traceOf(Fluid);
        const auto& parentFaces = surface.getPolytopeMap(surface.getDimension()).left;
        const VectorFunction faceNormal(static_cast<size_t>(3),
          [&](const Geometry::Point& point) {
            const Index parent = parentFaces[point.getPolytope().getIndex()];
            const Geometry::Point parentPoint(
              *mesh.getPolytope(mesh.getDimension() - 1, parent),
              point.getPhysicalCoordinates());
            const Math::SpatialVector<Real> normal = -fluidNormal.getValue(parentPoint);
            return normal;
          });
        Problem projection(projected, test);
        projection = Integral(length * length * Jacobian(projected), Jacobian(test)) +
          Integral(projected, test) - Integral(faceNormal, test);
        projection.assemble();
        solveDirect(projection);
        const auto& system = projection.getLinearSystem();
        const Real residual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(residual) || residual > LinearResidualTolerance)
          throw std::runtime_error("The thickness normal projection did not converge.");
        GridFunction normal(space);
        normal.getData().setZero();
        const auto& parentVertices = surface.getPolytopeMap(0).left;
        for (Index vertex = 0; vertex < surface.getVertexCount(); ++vertex)
        {
          const auto source = surfaceSpace.getDOFs(0, vertex);
          const auto target = space.getDOFs(0, parentVertices[vertex]);
          for (size_t component = 0; component < 3; ++component)
            normal.getData()(target(component)) =
              projected.getSolution().getData()(source(component));
        }
        return normal;
      }

      template <class ChamberMesh>
      ThicknessPenalty(const ChamberMesh& mesh, Real minimum)
        : m_minimum(minimum)
      {
        const size_t D = mesh.getDimension();
        auto fluidNormal = FaceNormal(mesh);
        fluidNormal.traceOf(Fluid);
        for (auto face = mesh.getPolytope(D - 1); face; ++face)
        {
          if (face->getAttribute() != Gamma)
            continue;
          const auto& vertices = face->getVertices();
          Triangle triangle;
          for (size_t k = 0; k < 3; ++k)
            triangle.x[k] = mesh.getVertexCoordinates(vertices[k]);
          const Math::SpatialVector<Real> midpoint =
            (triangle.x[0] + triangle.x[1] + triangle.x[2]) / 3;
          triangle.normal = -fluidNormal.getValue(Geometry::Point(*face, midpoint));
          triangle.normal.normalize();
          triangle.face = face->getIndex();
          for (size_t d = 0; d < 3; ++d)
          {
            triangle.lo[d] = std::min({triangle.x[0](d), triangle.x[1](d), triangle.x[2](d)});
            triangle.hi[d] = std::max({triangle.x[0](d), triangle.x[1](d), triangle.x[2](d)});
          }
          m_triangles.push_back(triangle);
        }
        m_order.resize(m_triangles.size());
        std::iota(m_order.begin(), m_order.end(), Index(0));
        if (!m_order.empty())
          build(0, m_order.size());
      }

      /**
       * @brief Adds the negative derivative of the first-exit penalty.
       *
       * The regularized directions are held fixed in this discrete
       * derivative; the source-face area and the ray-triangle crossing
       * are differentiated on the current triangulation.
       */
      template <class ChamberMesh, class Space, class Normal, class Load>
      Result evaluate(const ChamberMesh& mesh, const Space& space,
        const Normal& projectedNormal, Real weight, Load& load) const
      {
        static constexpr std::array<std::array<Real, 3>, 3> quadrature{
          {{Real(2) / 3, Real(1) / 6, Real(1) / 6},
            {Real(1) / 6, Real(2) / 3, Real(1) / 6},
            {Real(1) / 6, Real(1) / 6, Real(2) / 3}}};
        const auto& rotations = SewedOutput::getCubeRotations();
        const Real multiplicity = static_cast<Real>(rotations.size());
        const size_t D = mesh.getDimension();
        auto fluidNormal = FaceNormal(mesh);
        fluidNormal.traceOf(Fluid);
        const auto deposit = [&](Index face, const std::array<Real, 3>& barycentric,
                               const Math::SpatialVector<Real>& direction, Real value) {
          const auto& vertices = mesh.getPolytope(D - 1, face)->getVertices();
          for (size_t k = 0; k < 3; ++k)
          {
            const auto dofs = space.getDOFs(0, vertices[k]);
            for (Eigen::Index component = 0; component < 3; ++component)
              load(dofs(component)) += value * barycentric[k] * direction(component);
          }
        };

        Result result;
        result.curvature.assign(mesh.getVertexCount(), Real(0));
        std::vector<Real> curvatureMass(mesh.getVertexCount(), Real(0));
        for (auto face = mesh.getPolytope(D - 1); face; ++face)
        {
          if (face->getAttribute() != Gamma)
            continue;
          const auto& vertices = face->getVertices();
          const Math::SpatialVector<Real> a = mesh.getVertexCoordinates(vertices[0]);
          const Math::SpatialVector<Real> b = mesh.getVertexCoordinates(vertices[1]);
          const Math::SpatialVector<Real> c = mesh.getVertexCoordinates(vertices[2]);
          const Math::SpatialVector<Real> orientedArea = cross(b - a, c - a);
          const Real doubledArea = orientedArea.norm();
          if (!(doubledArea > 0))
            throw std::runtime_error("The thickness interface contains a degenerate face.");
          const Real area = doubledArea / 2;
          const Math::SpatialVector<Real> orientation = orientedArea / doubledArea;
          const std::array<Math::SpatialVector<Real>, 3> edgeOpposite{b - c, c - a, a - b};
          Math::SpatialMatrix<Real> gradientNormal(3, 3);
          gradientNormal.setZero();
          Math::SpatialMatrix<Real> gradientPosition(3, 3);
          gradientPosition.setZero();
          for (size_t k = 0; k < 3; ++k)
          {
            const auto dofs = space.getDOFs(0, vertices[k]);
            Math::SpatialVector<Real> value(3);
            for (size_t component = 0; component < 3; ++component)
              value(component) = projectedNormal.getData()(dofs(component));
            const Math::SpatialVector<Real> gradientBasis =
              cross(edgeOpposite[k], orientation) / doubledArea;
            for (size_t i = 0; i < 3; ++i)
              for (size_t j = 0; j < 3; ++j)
              {
                gradientNormal(i, j) += value(i) * gradientBasis(j);
                gradientPosition(i, j) +=
                  mesh.getVertexCoordinates(vertices[k])(i) * gradientBasis(j);
              }
          }
          if (std::abs(gradientPosition.trace() - Real(2)) > Real(1e-8))
            throw std::runtime_error("The surface barycentric gradient is inconsistent: " +
              std::to_string(gradientPosition.trace()));
          for (const auto& barycentric : quadrature)
          {
            const Math::SpatialVector<Real> s =
              barycentric[0] * a + barycentric[1] * b + barycentric[2] * c;
            const Geometry::Point surfacePoint(*face, s);
            // Use the P1 face trace of the fixed nodal field.
            Math::SpatialVector<Real> rawNormal =
              Math::SpatialVector<Real>::Zero(3);
            for (size_t k = 0; k < 3; ++k)
            {
              const auto dofs = space.getDOFs(0, vertices[k]);
              for (size_t component = 0; component < 3; ++component)
                rawNormal(component) +=
                  barycentric[k] * projectedNormal.getData()(dofs(component));
            }
            const Real magnitude = rawNormal.norm();
            if (!(std::isfinite(magnitude) && magnitude > Real(1e-12)))
              throw std::runtime_error("The smoothed thickness normal vanishes.");
            result.minimumNormalMagnitude = std::min(result.minimumNormalMagnitude, magnitude);
            result.maximumNormalMagnitude = std::max(result.maximumNormalMagnitude, magnitude);
            const Math::SpatialVector<Real> smoothNormal = rawNormal / magnitude;
            const Math::SpatialVector<Real> geometricNormal =
              -fluidNormal.getValue(surfacePoint);
            const Real alignment = smoothNormal.dot(geometricNormal);
            result.minimumNormalAlignment = std::min(result.minimumNormalAlignment, alignment);
            result.meanNormalAlignment += alignment;
            ++result.samples;
            const Real curvature =
              (gradientNormal.trace() -
                smoothNormal.dot(gradientNormal * smoothNormal)) / magnitude;
            result.minimumCurvature = std::min(result.minimumCurvature, curvature);
            result.maximumCurvature = std::max(result.maximumCurvature, curvature);
            for (size_t k = 0; k < 3; ++k)
            {
              const Real mass = area * barycentric[k] / 3;
              result.curvature[vertices[k]] += mass * curvature;
              curvatureMass[vertices[k]] += mass;
            }
            Math::SpatialVector<Real> normal = smoothNormal;
            if (alignment < Real(0.5))
            {
              ++result.correctedNormals;
              const Real blend = (Real(0.5) - alignment) / (Real(1) - alignment);
              normal = (Real(1) - blend) * smoothNormal + blend * geometricNormal;
              normal.normalize();
            }
            const Math::SpatialVector<Real> direction = -normal;
            const Exit exit = firstExit(s, direction);
            if (!exit.found || exit.distance >= m_minimum)
              continue;
            if (exit.transversality < Real(1e-3))
              throw std::runtime_error(
                "The first thickness exit is nearly tangent to the ray.");
            ++result.violating;
            result.minimumExit = std::min(result.minimumExit, exit.distance);
            result.minimumTransversality =
              std::min(result.minimumTransversality, exit.transversality);
            const Real deficit = m_minimum - exit.distance;
            result.maximumDeficit = std::max(result.maximumDeficit, deficit);
            const Real measure = multiplicity * area / Real(3);
            result.penalty += measure * deficit * deficit;
            const Triangle& target = m_triangles[exit.triangle];
            const Math::SpatialVector<Real> gradient =
              rotations[exit.rotation] * target.normal / exit.transversality;
            for (size_t k = 0; k < 3; ++k)
            {
              const Math::SpatialVector<Real> sourceGradient =
                Real(2) * deficit * barycentric[k] * gradient +
                (deficit * deficit / doubledArea) *
                  cross(edgeOpposite[k], orientation);
              const auto dofs = space.getDOFs(0, vertices[k]);
              for (Eigen::Index component = 0; component < 3; ++component)
                load(dofs(component)) -= weight * measure * sourceGradient(component);
            }
            deposit(target.face, exit.barycentric,
              rotations[exit.rotation].transpose() * gradient,
              weight * measure * Real(2) * deficit);
          }
        }
        if (result.samples)
          result.meanNormalAlignment /= result.samples;
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
          if (curvatureMass[vertex] > 0)
            result.curvature[vertex] /= curvatureMass[vertex];
        if (!std::isfinite(result.minimumExit))
          result.minimumExit = m_minimum;
        return result;
      }
    private:
      struct Triangle
      {
          std::array<Math::SpatialVector<Real>, 3> x;
          Math::SpatialVector<Real> normal;
          std::array<Real, 3> lo, hi;
          Index face;
      };

      struct Node
      {
          std::array<Real, 3> lo, hi;
          size_t first, count;
          Index left = 0, right = 0;
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

      Index build(size_t first, size_t end)
      {
        Node node;
        node.first = first;
        node.count = end - first;
        for (size_t d = 0; d < 3; ++d)
        {
          node.lo[d] = std::numeric_limits<Real>::infinity();
          node.hi[d] = -std::numeric_limits<Real>::infinity();
          for (size_t i = first; i < end; ++i)
          {
            node.lo[d] = std::min(node.lo[d], m_triangles[m_order[i]].lo[d]);
            node.hi[d] = std::max(node.hi[d], m_triangles[m_order[i]].hi[d]);
          }
        }
        const Index id = static_cast<Index>(m_nodes.size());
        m_nodes.push_back(node);
        if (node.count > 8)
        {
          size_t axis = 0;
          for (size_t d = 1; d < 3; ++d)
            if (node.hi[d] - node.lo[d] > node.hi[axis] - node.lo[axis])
              axis = d;
          const size_t middle = (first + end) / 2;
          std::nth_element(m_order.begin() + first, m_order.begin() + middle,
            m_order.begin() + end, [&](Index a, Index b) {
              return m_triangles[a].lo[axis] + m_triangles[a].hi[axis] <
                m_triangles[b].lo[axis] + m_triangles[b].hi[axis];
            });
          m_nodes[id].left = build(first, middle);
          m_nodes[id].right = build(middle, end);
        }
        return id;
      }

      struct Exit
      {
          Real distance = 0;
          Index triangle = 0;
          size_t rotation = 0;
          std::array<Real, 3> barycentric{0, 0, 0};
          Real transversality = 1;
          bool found = false;
      };

      static bool intersectsBox(const Math::SpatialVector<Real>& origin,
        const Math::SpatialVector<Real>& direction, const Node& node, Real length)
      {
        Real first = 0, last = length;
        for (size_t d = 0; d < 3; ++d)
        {
          if (std::abs(direction(d)) < Real(1e-14))
          {
            if (origin(d) < node.lo[d] - Real(1e-12) ||
                origin(d) > node.hi[d] + Real(1e-12))
              return false;
            continue;
          }
          const Real a = (node.lo[d] - Real(1e-12) - origin(d)) / direction(d);
          const Real b = (node.hi[d] + Real(1e-12) - origin(d)) / direction(d);
          first = std::max(first, std::min(a, b));
          last = std::min(last, std::max(a, b));
          if (first > last)
            return false;
        }
        return true;
      }

      static bool intersectsTriangle(const Math::SpatialVector<Real>& origin,
        const Math::SpatialVector<Real>& direction, const Triangle& triangle,
        Real minimumDistance, Real maximumDistance, Real& distance,
        std::array<Real, 3>& barycentric)
      {
        const Math::SpatialVector<Real> edge1 = triangle.x[1] - triangle.x[0];
        const Math::SpatialVector<Real> edge2 = triangle.x[2] - triangle.x[0];
        const Math::SpatialVector<Real> transverse = cross(direction, edge2);
        const Real determinant = edge1.dot(transverse);
        if (std::abs(determinant) <=
            Real(1e-12) * cross(edge1, edge2).norm())
          return false;
        const Real inverse = Real(1) / determinant;
        const Math::SpatialVector<Real> offset = origin - triangle.x[0];
        const Real u = offset.dot(transverse) * inverse;
        const Math::SpatialVector<Real> parallel = cross(offset, edge1);
        const Real v = direction.dot(parallel) * inverse;
        constexpr Real tolerance = 1e-10;
        if (u < -tolerance || v < -tolerance || u + v > Real(1) + tolerance)
          return false;
        const Real t = edge2.dot(parallel) * inverse;
        if (!(t > minimumDistance && t <= maximumDistance))
          return false;
        distance = t;
        barycentric = {Real(1) - u - v, u, v};
        return true;
      }

      /// The first sewn exit is the first chamber exit of an inverse-rotated ray.
      Exit firstExit(const Math::SpatialVector<Real>& origin,
        const Math::SpatialVector<Real>& direction) const
      {
        Exit exit;
        exit.distance = m_minimum;
        if (m_nodes.empty())
          return exit;
        const Real originTolerance =
          std::max(Real(1e-12), Real(1e-10) * m_minimum);
        std::vector<Index> pending;
        const auto& rotations = SewedOutput::getCubeRotations();
        for (size_t r = 0; r < rotations.size(); ++r)
        {
          const Math::SpatialVector<Real> localOrigin =
            rotations[r].transpose() * origin;
          const Math::SpatialVector<Real> localDirection =
            rotations[r].transpose() * direction;
          pending.push_back(0);
          while (!pending.empty())
          {
            const Node& node = m_nodes[pending.back()];
            pending.pop_back();
            if (!intersectsBox(localOrigin, localDirection, node, exit.distance))
              continue;
            if (node.count > 8)
            {
              pending.push_back(node.left);
              pending.push_back(node.right);
              continue;
            }
            for (size_t i = node.first; i < node.first + node.count; ++i)
            {
              const Index id = m_order[i];
              const Real transversality =
                m_triangles[id].normal.dot(localDirection);
              if (transversality <= Real(1e-12))
                continue;
              Real distance;
              std::array<Real, 3> barycentric;
              if (intersectsTriangle(localOrigin, localDirection, m_triangles[id],
                    originTolerance, exit.distance, distance, barycentric) &&
                  (!exit.found || distance < exit.distance - Real(1e-12) ||
                   (std::abs(distance - exit.distance) <= Real(1e-12) &&
                    transversality > exit.transversality)))
              {
                exit.found = true;
                exit.distance = distance;
                exit.triangle = id;
                exit.rotation = r;
                exit.barycentric = barycentric;
                exit.transversality = transversality;
              }
            }
          }
        }
        return exit;
      }
      Real m_minimum;
      std::vector<Triangle> m_triangles;
      std::vector<Index> m_order;
      std::vector<Node> m_nodes;
  };
}

#endif
