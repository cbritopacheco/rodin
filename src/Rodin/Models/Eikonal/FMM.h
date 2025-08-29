/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_EIKONAL_FMM_H
#define RODIN_MODELS_EIKONAL_FMM_H

/**
 * @file FMM.h
 * @brief Fast marching method on simplicial meshes
 */

#include <queue>
#include <limits>
#include <algorithm>
#include <functional>
#include <vector>
#include <cmath>
#include <cassert>

#include "Rodin/Context/Local.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/PolytopeIterator.h"
#include "Rodin/Variational/P1/P1.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Models::Eikonal
{
  /**
   * @brief Fast marching method for solving the eikonal equation
   */
  template <class Solution, class SpeedFunction>
  class FMM;

  template <class Data, class SpeedFunction>
  class FMM<Variational::GridFunction<Variational::P1<Real, Geometry::Mesh<Context::Local>>, Data>, SpeedFunction>
  {
    public:
      using ScalarType = Real;
      using Context = Context::Local;
      using Mesh = Geometry::Mesh<Context>;
      using FES = Variational::P1<ScalarType, Mesh>;
      using SolutionType = Variational::GridFunction<FES, Math::Vector<ScalarType>>;
      using SpeedFunctionType = SpeedFunction;

      enum class Label : uint8_t
      {
        Far, Considered, Accepted
      };

    private:
      struct PQItem
      {
        Index nodeIndex;
        Real value;
        bool operator>(const PQItem& other) const { return value > other.value; }
      };

      using PriorityQueue =
        std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>>;

    public:
      template <class Callable>
      FMM(SolutionType& u, Callable&& speed)
        : m_u(u), m_speed(std::forward<Callable>(speed))
      {}

      template <class Seed>
      FMM& seed(Seed&& seed)
      {
        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& conn = mesh.getConnectivity();
        const Index nV = fes.getSize();
        m_labels.assign(nV, Label::Far);
        u = std::numeric_limits<Real>::infinity();
        for (Index s : std::forward<Seed>(seed))
        {
          if (s >= nV)
            continue;
          m_labels[s] = Label::Accepted;
          u[s] = 0.0;
          assert(!std::isnan(u[s]));

          const auto& N = conn.getIncidence(0, 0).at(s);
          for (Index nb : N)
          {
            if (m_labels[nb] != Label::Far)
              continue;
            const Real arr = local(nb, u, mesh);
            if (std::isnan(arr))
              continue;
            assert(!std::isnan(arr));
            if (arr < u[nb])
            {
              u[nb] = arr;
              assert(!std::isnan(u[nb]));
              m_labels[nb] = Label::Considered;
              m_pq.push({ nb, arr });
            }
          }
        }
        return *this;
      }

      void solve()
      {
        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const int D = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();
        const Index nV = fes.getSize();

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 0, 0);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 0, D);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, 0);

        // March
        while (!m_pq.empty())
        {
          const auto cur = m_pq.top(); m_pq.pop();
          const Index i = cur.nodeIndex;

          if (i >= nV)
            continue;
          if (std::isnan(cur.value))
            continue;
          if (cur.value != u[i])
            continue; // stale-key skip
          if (m_labels[i] == Label::Accepted)
            continue;

          m_labels[i] = Label::Accepted;

          const auto& N = conn.getIncidence(0, 0).at(i);
          for (Index j : N)
          {
            if (m_labels[j] == Label::Accepted)
              continue;
            const Real arr = local(j, u, mesh);
            if (std::isnan(arr))
              continue;
            assert(!std::isnan(arr));
            if (arr < u[j]) {
              u[j] = arr;
              assert(!std::isnan(u[j]));
              if (m_labels[j] == Label::Far) m_labels[j] = Label::Considered;
              m_pq.push({ j, arr });
            }
          }
        }
      }

    private:
      Real local(Index p, const SolutionType& u, const Mesh& mesh) const
      {
        static thread_local Math::SpatialPoint s_dummy(mesh.getSpaceDimension());

        const int D = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();
        const auto& v2c  = conn.getIncidence(0, D);
        const auto& c2v  = conn.getIncidence(D, 0);

        Real best = std::numeric_limits<Real>::infinity();

        for (Index cid : v2c.at(p))
        {
          const auto& verts = c2v.at(cid);

          Index A[3]; int k = 0;
          for (Index q : verts)
          {
            if (q == p)
              continue;

            if (m_labels[q] == Label::Accepted)
              A[k++] = q;

            if (k == D)
              break;
          }

          if (k == 0)
            continue;

          const Geometry::VertexIterator vertex = mesh.getVertex(p);
          const Real s =
            m_speed(Geometry::Point(*vertex, s_dummy, vertex->getCoordinates()));
          const Real F =
            std::isnan(s) ? std::numeric_limits<Real>::infinity() : std::max(std::numeric_limits<Real>::min(), s);

          if (std::isnan(F) || F <= 0)
            continue;

          assert(!std::isnan(F));

          if (k == 1)
          {
            // Use improved distance computation
            const Real d = computeEdgeDistance(p, A[0], mesh);
            if (std::isnan(d))
              continue;
            assert(!std::isnan(d));
            best = std::min(best, u[A[0]] + d / F);
          }
          else if (D == 1) // 1D curve mesh (possibly embedded in higher dimensions)
          {
            // For 1D meshes, only use single neighbor updates
            const Real d = computeEdgeDistance(p, A[0], mesh);
            if (!std::isnan(d))
              best = std::min(best, u[A[0]] + d / F);
          }
          else if (D == 2) // 2D surface mesh (possibly embedded in 3D)
          {
            const Real t = surfaceTriangleUpdate(p, A[0], A[1], u, mesh, F);
            if (!std::isnan(t))
              best = std::min(best, t);
          }
          else if (D == 3) // 3D volume mesh
          {
            const Real t2 = volumeTetrahedronUpdate(p, A[0], A[1], u, mesh, F);
            if (!std::isnan(t2))
              best = std::min(best, t2);
            if (k == 3)
            {
              const Real t3 = volumeTetrahedronUpdate(p, A[0], A[1], A[2], u, mesh, F);
              if (!std::isnan(t3))
                best = std::min(best, t3);
            }
          }
        }
        assert(!std::isnan(best));
        return best;
      }

      static Real computeGeometricDistance(Index a, Index b, const Mesh& mesh)
      {
        // For basic geometric distance, use the existing vertex coordinates approach
        // but let the mesh handle any embedding properly
        const auto xa = mesh.getVertexCoordinates(a);
        const auto xb = mesh.getVertexCoordinates(b);

        // This automatically handles surface meshes and embedded coordinates
        const Real distance = (xa - xb).norm();
        assert(!std::isnan(distance));
        return distance;
      }

      static Real computeEdgeDistance(Index a, Index b, const Mesh& mesh)
      {
        // Check if there's a direct edge between vertices a and b
        // If so, compute the geodesic distance along the edge
        const auto& conn = mesh.getConnectivity();
        const auto& v2e = conn.getIncidence(0, 1); // vertex to edge

        // Find common edge between vertices a and b
        if (v2e.size() > a && v2e.size() > b)
        {
          const auto& edges_a = v2e.at(a);
          const auto& edges_b = v2e.at(b);

          for (Index edge_a : edges_a)
          {
            for (Index edge_b : edges_b)
            {
              if (edge_a == edge_b)
              {
                // Found common edge, compute geodesic distance along it
                return computeEdgeGeodesicDistance(a, b, edge_a, mesh);
              }
            }
          }
        }
        // No direct edge, fallback to geometric distance
        return computeGeometricDistance(a, b, mesh);
      }

      static Real computeEdgeGeodesicDistance(
          Index a, Index b, Index edge_idx, const Mesh& mesh)
      {
        // For now, use geometric distance as approximation
        // In the future, this could integrate along curved edges using PolytopeTransformation
        return computeGeometricDistance(a, b, mesh);
      }

      // Surface triangle update - works for 2D triangles embedded in any dimension
      Real surfaceTriangleUpdate(Index p, Index i, Index j,
                        const SolutionType& u, const Mesh& mesh, Real F) const
      {
        // Get vertex coordinates using the mesh's coordinate access
        const auto xp = mesh.getVertexCoordinates(p);
        const auto xi = mesh.getVertexCoordinates(i);
        const auto xj = mesh.getVertexCoordinates(j);

        // Compute edge vectors and lengths
        const auto vi = xi - xp;  // Vector from p to i
        const auto vj = xj - xp;  // Vector from p to j

        const Real a = vi.norm();  // Distance from p to i
        const Real b = vj.norm();  // Distance from p to j
        if (!(a > 0 && b > 0)) return std::numeric_limits<Real>::infinity();
        assert(!std::isnan(a) && !std::isnan(b));

        // Compute cosine of angle at p using dot product
        const Real cos_angle = vi.dot(vj) / (a * b);
        assert(!std::isnan(cos_angle));

        const Real ui = u[i], uj = u[j];
        assert(!std::isnan(ui) && !std::isnan(uj));

        // Solve quadratic equation for travel time
        // This formulation works for surfaces embedded in any dimension
        const Real A = 1 / (a * a) + 1 / (b * b) - 2 * cos_angle / (a * b);
        const Real B = -2.0 * (ui / (a * a) + uj / (b * b) - (cos_angle / (a * b)) * (ui + uj));
        const Real C = (ui * ui) / (a * a) + (uj * uj) / (b * b) - 2 * cos_angle * ui * uj / (a * b) - 1.0 / (F * F);
        assert(!std::isnan(A) && !std::isnan(B) && !std::isnan(C));

        const Real disc = B * B - 4 * A * C;
        if (disc < 0 || std::isnan(disc))
          return std::numeric_limits<Real>::infinity();

        Real t = (-B + std::sqrt(disc)) / (2 * A);
        if (std::isnan(t))
          return std::numeric_limits<Real>::infinity();

        const Real umax = std::max(ui, uj);
        if (t < umax)
          t = std::min(ui + a / F, uj + b / F);
        assert(!std::isnan(t));
        return t;
      }

      // 3D volume element 2-neighbor update using geometric infrastructure
      Real volumeTetrahedronUpdate(
          Index p, Index i, Index j,
          const SolutionType& u, const Mesh& mesh, Real F) const
      {
        // Get vertex coordinates using the mesh's coordinate access
        const auto xp = mesh.getVertexCoordinates(p);
        const auto xi = mesh.getVertexCoordinates(i);
        const auto xj = mesh.getVertexCoordinates(j);

        const auto d1 = xi - xp;
        const auto d2 = xj - xp;

        const Real g11 = d1.dot(d1);
        const Real g22 = d2.dot(d2);
        const Real g12 = d1.dot(d2);
        if (!(g11 > 0 && g22 > 0))
          return std::numeric_limits<Real>::infinity();
        assert(!std::isnan(g11) && !std::isnan(g22) && !std::isnan(g12));

        const Real det = g11 * g22 - g12 * g12;
        if (det <= 0 || std::isnan(det))
          return std::numeric_limits<Real>::infinity();

        const Real inv11 =  g22 / det;
        const Real inv22 =  g11 / det;
        const Real inv12 = -g12 / det;
        assert(!std::isnan(inv11) && !std::isnan(inv22) && !std::isnan(inv12));

        const Real ui = u[i], uj = u[j];
        assert(!std::isnan(ui) && !std::isnan(uj));

        const Real alpha = (inv11 + 2 * inv12 + inv22);
        const Real beta  = (inv11 * ui + inv12 * uj + inv12 * ui + inv22 * uj);
        const Real gamma = (inv11 * ui * ui + 2 * inv12 * ui * uj + inv22 * uj * uj);
        assert(!std::isnan(alpha) && !std::isnan(beta) && !std::isnan(gamma));

        const Real A = alpha;
        const Real B = -2 * beta;
        const Real C = gamma - 1.0 / (F * F);
        assert(!std::isnan(A) && !std::isnan(B) && !std::isnan(C));

        const Real disc = B * B - 4 * A * C;
        if (disc < 0 || std::isnan(disc))
          return std::numeric_limits<Real>::infinity();

        Real t = (-B + std::sqrt(disc)) / (2 * A);
        if (std::isnan(t))
          return std::numeric_limits<Real>::infinity();

        const Real umax = std::max(ui, uj);
        if (t < umax)
        {
          const Real a = std::sqrt(g11);
          const Real b = std::sqrt(g22);
          assert(!std::isnan(a) && !std::isnan(b));
          t = std::min(ui + a / F, uj + b / F);
        }
        assert(!std::isnan(t));
        return t;
      }

      // 3D volume element 3-neighbor update using geometric infrastructure  
      Real volumeTetrahedronUpdate(
          Index p, Index i, Index j, Index k,
          const SolutionType& u, const Mesh& mesh, Real F) const
      {
        // Get vertex coordinates using the mesh's coordinate access
        const auto xp = mesh.getVertexCoordinates(p);
        const auto xi = mesh.getVertexCoordinates(i);
        const auto xj = mesh.getVertexCoordinates(j);
        const auto xk = mesh.getVertexCoordinates(k);

        const auto d1 = xi - xp;
        const auto d2 = xj - xp;
        const auto d3 = xk - xp;

        const Real g11 = d1.dot(d1);
        const Real g22 = d2.dot(d2);
        const Real g33 = d3.dot(d3);
        const Real g12 = d1.dot(d2);
        const Real g13 = d1.dot(d3);
        const Real g23 = d2.dot(d3);
        assert(!std::isnan(g11) && !std::isnan(g22) && !std::isnan(g33));
        assert(!std::isnan(g12) && !std::isnan(g13) && !std::isnan(g23));

        const Real det =
          g11 * (g22 * g33 - g23 * g23)
          - g12 * (g12 * g33 - g13 * g23)
          + g13 * (g12 * g23 - g13 * g22);

        if (det <= 0 || std::isnan(det))
          return std::numeric_limits<Real>::infinity();

        const Real c11 =  (g22 * g33 - g23 * g23);
        const Real c22 =  (g11 * g33 - g13 * g13);
        const Real c33 =  (g11 * g22 - g12 * g12);
        const Real c12 = -(g12 * g33 - g13 * g23);
        const Real c13 =  (g12 * g23 - g13 * g22);
        const Real c23 = -(g11 * g23 - g12 * g13);

        const Real inv11 = c11 / det, inv22 = c22 / det, inv33 = c33 / det;
        const Real inv12 = c12 / det, inv13 = c13 / det, inv23 = c23 / det;
        assert(!std::isnan(inv11) && !std::isnan(inv22) && !std::isnan(inv33));
        assert(!std::isnan(inv12) && !std::isnan(inv13) && !std::isnan(inv23));

        const Real ui = u[i], uj = u[j], uk = u[k];
        assert(!std::isnan(ui) && !std::isnan(uj) && !std::isnan(uk));

        const Real alpha = inv11 + 2 * inv12 + 2 * inv13 + inv22 + 2 * inv23 + inv33;
        const Real beta  = (inv11 * ui + inv12 * uj + inv13 * uk)
                         + (inv12 * ui + inv22 * uj + inv23 * uk)
                         + (inv13 * ui + inv23 * uj + inv33 * uk);
        const Real gamma = inv11 * ui * ui + 2 * inv12 * ui * uj + 2 * inv13 * ui * uk
                         + inv22 * uj * uj + 2 * inv23 * uj * uk + inv33 * uk * uk;
        assert(!std::isnan(alpha) && !std::isnan(beta) && !std::isnan(gamma));

        const Real A = alpha;
        const Real B = -2 * beta;
        const Real C = gamma - 1.0 / (F * F);
        assert(!std::isnan(A) && !std::isnan(B) && !std::isnan(C));

        const Real disc = B * B - 4 * A * C;
        if (disc < 0 || std::isnan(disc)) return std::numeric_limits<Real>::infinity();

        Real t = (-B + std::sqrt(disc)) / (2 * A);
        if (std::isnan(t)) return std::numeric_limits<Real>::infinity();

        const Real umax = std::max(ui, std::max(uj, uk));
        if (t < umax)
        {
          Real best = std::numeric_limits<Real>::infinity();
          const Real ni = d1.norm();
          const Real nj = d2.norm();
          const Real nk = d3.norm();
          assert(!std::isnan(ni) && !std::isnan(nj) && !std::isnan(nk));
          best = std::min(best, ui + ni / F);
          best = std::min(best, uj + nj / F);
          best = std::min(best, uk + nk / F);
          t = best;
        }
        assert(!std::isnan(t));
        return t;
      }

    private:
      std::reference_wrapper<SolutionType> m_u;
      SpeedFunctionType m_speed;

      PriorityQueue m_pq;
      std::vector<Label> m_labels;
  };

  template <class Solution, class SpeedFunction>
  FMM(Solution& u, SpeedFunction&& speed) -> FMM<Solution, SpeedFunction>;
}

#endif
