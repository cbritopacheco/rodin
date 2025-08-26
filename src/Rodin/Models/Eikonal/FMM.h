/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_EIKONAL_FMM_H
#define RODIN_MODELS_EIKONAL_FMM_H

#include <queue>
#include <limits>
#include <algorithm>
#include <functional>
#include <vector>
#include <cmath>
#include <cassert>

#include "Rodin/Context/Local.h"
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Variational/P1/P1.h"
#include "Rodin/Variational/GridFunction.h"

namespace Rodin::Models::Eikonal
{
  template <class Solution, class SpeedFunction>
  class FMM;

  template <class Data, class SpeedFunction>
  class FMM<Variational::GridFunction<Variational::P1<Real, Geometry::Mesh<Context::Local>>, Data>, SpeedFunction>
  {
    public:
      using ScalarType        = Real;
      using Context           = Context::Local;
      using Mesh              = Geometry::Mesh<Context>;
      using FES               = Variational::P1<ScalarType, Mesh>;
      using SolutionType      = Variational::GridFunction<FES, Math::Vector<ScalarType>>;
      using SpeedFunctionType = SpeedFunction;

      enum class Label { Far, Considered, Accepted };

    private:
      struct PQItem
      {
        Index nodeIndex;
        Real  value;
        bool operator>(const PQItem& other) const { return value > other.value; }
      };

      using PriorityQueue =
        std::priority_queue<PQItem, std::vector<PQItem>, std::greater<PQItem>>;

    public:
      template <class Callable>
      FMM(SolutionType& u, Callable&& speed)
        : m_u(u), m_speed(std::forward<Callable>(speed))
      {}

      template <class Container>
      FMM& setInterface(Container&& seeds)
      {
        m_interface.assign(std::begin(seeds), std::end(seeds));
        return *this;
      }

      void solve()
      {
        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const int D = mesh.getDimension();
        const auto& conn = mesh.getConnectivity();

        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 0, 0);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, 0, D);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, 0);

        const Index nV = fes.getSize();
        m_labels.assign(nV, Label::Far);
        u = std::numeric_limits<Real>::infinity();

        PriorityQueue pq;

        // Seeds
        for (Index s : m_interface)
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
            if (std::isnan(arr)) continue;
            assert(!std::isnan(arr));
            if (arr < u[nb])
            {
              u[nb] = arr;
              assert(!std::isnan(u[nb]));
              m_labels[nb] = Label::Considered;
              pq.push({ nb, arr });
            }
          }
        }

        // March
        while (!pq.empty())
        {
          const auto cur = pq.top(); pq.pop();
          const Index i = cur.nodeIndex;
          if (i >= nV) continue;

          if (std::isnan(cur.value)) continue; 
          if (cur.value != u[i]) continue; // stale-key skip
          if (m_labels[i] == Label::Accepted) continue;
          m_labels[i] = Label::Accepted;

          const auto& N = conn.getIncidence(0,0).at(i);
          for (Index j : N)
          {
            if (m_labels[j] == Label::Accepted) continue;
            const Real arr = local(j, u, mesh);
            if (std::isnan(arr)) continue;
            assert(!std::isnan(arr));
            if (arr < u[j]) {
              u[j] = arr;
              assert(!std::isnan(u[j]));
              if (m_labels[j] == Label::Far) m_labels[j] = Label::Considered;
              pq.push({ j, arr });
            }
          }
        }
      }

    private:
      Real local(Index p, const SolutionType& u, const Mesh& mesh) const
      {
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
            if (q == p) continue;
            if (m_labels[q] == Label::Accepted) A[k++] = q;
            if (k == D) break;
          }
          if (k == 0) continue;

          const Real F = speedAtVertex(p, mesh);
          if (std::isnan(F) || F <= 0) continue;
          assert(!std::isnan(F));

          if (k == 1)
          {
            const Real d = dist(mesh, p, A[0]);
            if (std::isnan(d)) continue;
            assert(!std::isnan(d));
            best = std::min(best, u[A[0]] + d / F);
          }
          else if (D == 2)
          {
            const Real t = tri2D_update(p, A[0], A[1], u, mesh, F);
            if (!std::isnan(t)) best = std::min(best, t);
          }
          else // D == 3
          {
            const Real t2 = gram_update_3D(p, A[0], A[1], u, mesh, F);
            if (!std::isnan(t2)) best = std::min(best, t2);
            if (k == 3)
            {
              const Real t3 = gram_update_3D(p, A[0], A[1], A[2], u, mesh, F);
              if (!std::isnan(t3)) best = std::min(best, t3);
            }
          }
        }
        assert(!std::isnan(best));
        return best;
      }

      static Real dist(const Mesh& mesh, Index a, Index b)
      {
        const auto& xa = mesh.getVertexCoordinates(a);
        const auto& xb = mesh.getVertexCoordinates(b);
        const Real n = (xa - xb).norm();
        assert(!std::isnan(n));
        return n;
      }

      Real speedAtVertex(Index p, const Mesh& mesh) const
      {
        const auto& x = mesh.getVertexCoordinates(p);
        const Real s = m_speed(x);
        if (std::isnan(s)) return std::numeric_limits<Real>::infinity();
        assert(!std::isnan(s));
        return std::max(std::numeric_limits<Real>::min(), s);
      }

      // 2D triangle update
      Real tri2D_update(Index p, Index i, Index j,
                        const SolutionType& u, const Mesh& mesh, Real F) const
      {
        const auto xp = mesh.getVertexCoordinates(p);
        const auto xi = mesh.getVertexCoordinates(i);
        const auto xj = mesh.getVertexCoordinates(j);

        const Real a = (xp - xi).norm();
        const Real b = (xp - xj).norm();
        if (!(a > 0 && b > 0)) return std::numeric_limits<Real>::infinity();
        assert(!std::isnan(a) && !std::isnan(b));

        const Real c = (xp - xi).dot(xp - xj) / (a * b);
        assert(!std::isnan(c));

        const Real ui = u[i], uj = u[j];
        assert(!std::isnan(ui) && !std::isnan(uj));

        const Real A = 1/(a*a) + 1/(b*b) + 2*c/(a*b);
        const Real B = -2.0 * ( ui/(a*a) + uj/(b*b) + (c/(a*b))*(ui + uj) );
        const Real C = (ui*ui)/(a*a) + (uj*uj)/(b*b) + 2*c*ui*uj/(a*b) - 1.0/(F*F);
        assert(!std::isnan(A) && !std::isnan(B) && !std::isnan(C));

        const Real disc = B*B - 4*A*C;
        if (disc < 0 || std::isnan(disc)) return std::numeric_limits<Real>::infinity();

        Real t = (-B + std::sqrt(disc)) / (2*A);
        if (std::isnan(t)) return std::numeric_limits<Real>::infinity();

        const Real umax = std::max(ui, uj);
        if (t < umax) t = std::min(ui + a/F, uj + b/F);
        assert(!std::isnan(t));
        return t;
      }

      // 3D 2-neighbor update
      Real gram_update_3D(Index p, Index i, Index j,
                          const SolutionType& u, const Mesh& mesh, Real F) const
      {
        const auto xp = mesh.getVertexCoordinates(p);
        const auto xi = mesh.getVertexCoordinates(i);
        const auto xj = mesh.getVertexCoordinates(j);

        const auto d1 = xp - xi;
        const auto d2 = xp - xj;

        const Real g11 = d1.dot(d1);
        const Real g22 = d2.dot(d2);
        const Real g12 = d1.dot(d2);
        if (!(g11 > 0 && g22 > 0)) return std::numeric_limits<Real>::infinity();
        assert(!std::isnan(g11) && !std::isnan(g22) && !std::isnan(g12));

        const Real det = g11*g22 - g12*g12;
        if (det <= 0 || std::isnan(det)) return std::numeric_limits<Real>::infinity();

        const Real inv11 =  g22 / det;
        const Real inv22 =  g11 / det;
        const Real inv12 = -g12 / det;
        assert(!std::isnan(inv11) && !std::isnan(inv22) && !std::isnan(inv12));

        const Real ui = u[i], uj = u[j];
        assert(!std::isnan(ui) && !std::isnan(uj));

        const Real alpha = (inv11 + 2*inv12 + inv22);
        const Real beta  = (inv11*ui + inv12*uj + inv12*ui + inv22*uj);
        const Real gamma = (inv11*ui*ui + 2*inv12*ui*uj + inv22*uj*uj);
        assert(!std::isnan(alpha) && !std::isnan(beta) && !std::isnan(gamma));

        const Real A = alpha;
        const Real B = -2*beta;
        const Real C = gamma - 1.0/(F*F);
        assert(!std::isnan(A) && !std::isnan(B) && !std::isnan(C));

        const Real disc = B*B - 4*A*C;
        if (disc < 0 || std::isnan(disc)) return std::numeric_limits<Real>::infinity();

        Real t = (-B + std::sqrt(disc)) / (2*A);
        if (std::isnan(t)) return std::numeric_limits<Real>::infinity();

        const Real umax = std::max(ui, uj);
        if (t < umax) {
          const Real a = std::sqrt(g11);
          const Real b = std::sqrt(g22);
          assert(!std::isnan(a) && !std::isnan(b));
          t = std::min(ui + a/F, uj + b/F);
        }
        assert(!std::isnan(t));
        return t;
      }

      // 3D 3-neighbor update
      Real gram_update_3D(Index p, Index i, Index j, Index k,
                          const SolutionType& u, const Mesh& mesh, Real F) const
      {
        const auto xp = mesh.getVertexCoordinates(p);
        const auto xi = mesh.getVertexCoordinates(i);
        const auto xj = mesh.getVertexCoordinates(j);
        const auto xk = mesh.getVertexCoordinates(k);

        const auto d1 = xp - xi;
        const auto d2 = xp - xj;
        const auto d3 = xp - xk;

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

        if (det <= 0 || std::isnan(det)) return std::numeric_limits<Real>::infinity();

        const Real c11 =  (g22 * g33 - g23*g23);
        const Real c22 =  (g11 * g33 - g13*g13);
        const Real c33 =  (g11 * g22 - g12*g12);
        const Real c12 = -(g12 * g33 - g13*g23);
        const Real c13 =  (g12 * g23 - g13*g22);
        const Real c23 = -(g11 * g23 - g12*g13);

        const Real inv11 = c11 / det, inv22 = c22 / det, inv33 = c33 / det;
        const Real inv12 = c12 / det, inv13 = c13 / det, inv23 = c23 / det;
        assert(!std::isnan(inv11) && !std::isnan(inv22) && !std::isnan(inv33));
        assert(!std::isnan(inv12) && !std::isnan(inv13) && !std::isnan(inv23));

        const Real ui = u[i], uj = u[j], uk = u[k];
        assert(!std::isnan(ui) && !std::isnan(uj) && !std::isnan(uk));

        const Real alpha = inv11 + 2*inv12 + 2*inv13 + inv22 + 2*inv23 + inv33;
        const Real beta  = (inv11*ui + inv12*uj + inv13*uk)
                         + (inv12*ui + inv22*uj + inv23*uk)
                         + (inv13*ui + inv23*uj + inv33*uk);
        const Real gamma = inv11*ui*ui + 2*inv12*ui*uj + 2*inv13*ui*uk
                         + inv22*uj*uj + 2*inv23*uj*uk + inv33*uk*uk;
        assert(!std::isnan(alpha) && !std::isnan(beta) && !std::isnan(gamma));

        const Real A = alpha;
        const Real B = -2*beta;
        const Real C = gamma - 1.0/(F*F);
        assert(!std::isnan(A) && !std::isnan(B) && !std::isnan(C));

        const Real disc = B*B - 4*A*C;
        if (disc < 0 || std::isnan(disc)) return std::numeric_limits<Real>::infinity();

        Real t = (-B + std::sqrt(disc)) / (2*A);
        if (std::isnan(t)) return std::numeric_limits<Real>::infinity();

        const Real umax = std::max(ui, std::max(uj, uk));
        if (t < umax)
        {
          Real best = std::numeric_limits<Real>::infinity();
          const Real ni = (xp - xi).norm();
          const Real nj = (xp - xj).norm();
          const Real nk = (xp - xk).norm();
          assert(!std::isnan(ni) && !std::isnan(nj) && !std::isnan(nk));
          best = std::min(best, ui + ni/F);
          best = std::min(best, uj + nj/F);
          best = std::min(best, uk + nk/F);
          t = best;
        }
        assert(!std::isnan(t));
        return t;
      }

    private:
      std::reference_wrapper<SolutionType> m_u;
      SpeedFunctionType m_speed;

      std::vector<Label> m_labels;
      std::vector<Index> m_interface;
  };

  template <class Solution, class SpeedFunction>
  FMM(Solution& u, SpeedFunction&& speed) -> FMM<Solution, SpeedFunction>;
}

#endif
