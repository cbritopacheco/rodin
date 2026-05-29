/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MinSTCut.h"

#include <algorithm>
#include <deque>
#include <stdexcept>

namespace Rodin::Geometry
{
  namespace
  {
    class Dinic
    {
      public:
        struct Arc
        {
          Index to;
          Index rev;
          Real cap;
        };

        explicit Dinic(Index n)
          : m_graph(n), m_level(n), m_next(n)
        {}

        void addDirected(Index from, Index to, Real capacity)
        {
          if (capacity < 0)
            throw std::invalid_argument("MinSTCut received a negative graph capacity.");
          Arc fwd{ to, m_graph[to].size(), capacity };
          Arc rev{ from, m_graph[from].size(), 0 };
          m_graph[from].push_back(fwd);
          m_graph[to].push_back(rev);
        }

        void addUndirected(Index a, Index b, Real capacity)
        {
          addDirected(a, b, capacity);
          addDirected(b, a, capacity);
        }

        Real maxFlow(Index source, Index sink)
        {
          Real flow = 0;
          while (buildLevelGraph(source, sink))
          {
            std::fill(m_next.begin(), m_next.end(), 0);
            while (true)
            {
              const Real pushed =
                augment(source, sink, std::numeric_limits<Real>::infinity());
              if (pushed <= 0)
                break;
              flow += pushed;
            }
          }
          return flow;
        }

        std::vector<Boolean> reachableFrom(Index source) const
        {
          std::vector<Boolean> seen(m_graph.size(), false);
          std::deque<Index> queue;
          seen[source] = true;
          queue.push_back(source);
          while (!queue.empty())
          {
            const Index current = queue.front();
            queue.pop_front();
            for (const auto& arc : m_graph[current])
            {
              if (arc.cap > 0 && !seen[arc.to])
              {
                seen[arc.to] = true;
                queue.push_back(arc.to);
              }
            }
          }
          return seen;
        }

      private:
        Boolean buildLevelGraph(Index source, Index sink)
        {
          std::fill(m_level.begin(), m_level.end(), -1);
          std::deque<Index> queue;
          m_level[source] = 0;
          queue.push_back(source);
          while (!queue.empty())
          {
            const Index current = queue.front();
            queue.pop_front();
            for (const auto& arc : m_graph[current])
            {
              if (arc.cap > 0 && m_level[arc.to] < 0)
              {
                m_level[arc.to] = m_level[current] + 1;
                queue.push_back(arc.to);
              }
            }
          }
          return m_level[sink] >= 0;
        }

        Real augment(Index current, Index sink, Real amount)
        {
          if (current == sink)
            return amount;
          for (Index& i = m_next[current]; i < m_graph[current].size(); ++i)
          {
            Arc& arc = m_graph[current][i];
            if (arc.cap <= 0 || m_level[arc.to] != m_level[current] + 1)
              continue;
            const Real pushed = augment(arc.to, sink, std::min(amount, arc.cap));
            if (pushed > 0)
            {
              arc.cap -= pushed;
              m_graph[arc.to][arc.rev].cap += pushed;
              return pushed;
            }
          }
          return 0;
        }

        std::vector<std::vector<Arc>> m_graph;
        std::vector<Integer> m_level;
        std::vector<Index> m_next;
    };
  }

  Real MinSTCut::getInsideCost(Real volume, Real moment) noexcept
  {
    return volume * std::max<Real>(0, moment);
  }

  Real MinSTCut::getOutsideCost(Real volume, Real moment) noexcept
  {
    return volume * std::max<Real>(0, -moment);
  }

  MinSTCut::Result MinSTCut::classify(
      const std::vector<Real>& volumes,
      const std::vector<Real>& moments,
      const std::vector<Edge>& edges) const
  {
    if (volumes.size() != moments.size())
      throw std::invalid_argument("MinSTCut volumes and moments have different sizes.");

    std::vector<Real> insideCosts(volumes.size());
    std::vector<Real> outsideCosts(volumes.size());
    for (Index i = 0; i < volumes.size(); ++i)
    {
      if (volumes[i] < 0)
        throw std::invalid_argument("MinSTCut received a negative cell volume.");
      insideCosts[i] = getInsideCost(volumes[i], moments[i]);
      outsideCosts[i] = getOutsideCost(volumes[i], moments[i]);
    }
    return solve(insideCosts, outsideCosts, edges);
  }

  MinSTCut::Result MinSTCut::solve(
      const std::vector<Real>& insideCosts,
      const std::vector<Real>& outsideCosts,
      const std::vector<Edge>& edges) const
  {
    if (insideCosts.size() != outsideCosts.size())
      throw std::invalid_argument("MinSTCut unary cost arrays have different sizes.");

    const Index cellCount = insideCosts.size();
    const Index source = cellCount;
    const Index sink = cellCount + 1;

    Dinic graph(cellCount + 2);
    for (Index i = 0; i < cellCount; ++i)
    {
      if (insideCosts[i] < 0 || outsideCosts[i] < 0)
        throw std::invalid_argument("MinSTCut received a negative unary cost.");

      graph.addDirected(source, i, outsideCosts[i]);
      graph.addDirected(i, sink, insideCosts[i]);
    }

    for (const Edge& edge : edges)
    {
      if (edge.first >= cellCount || edge.second >= cellCount)
        throw std::out_of_range("MinSTCut edge references a cell outside the graph.");
      graph.addUndirected(edge.first, edge.second, edge.capacity);
    }

    graph.maxFlow(source, sink);
    const auto reachable = graph.reachableFrom(source);

    Result result;
    result.labels.assign(cellCount, Outside);
    for (Index i = 0; i < cellCount; ++i)
    {
      if (reachable[i])
      {
        result.labels[i] = Inside;
        result.insideCells.push_back(i);
        result.energy += insideCosts[i];
      }
      else
      {
        result.labels[i] = Outside;
        result.outsideCells.push_back(i);
        result.energy += outsideCosts[i];
      }
    }

    for (const Edge& edge : edges)
    {
      if (result.labels[edge.first] != result.labels[edge.second])
      {
        result.cutEdges.push_back(edge);
        result.energy += edge.capacity;
      }
    }
    return result;
  }
}
