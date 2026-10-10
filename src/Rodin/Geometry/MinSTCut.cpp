/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "MinSTCut.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>

#include "Rodin/Alert/MemberFunctionException.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
  MinSTCut<LocalMesh>::Dinic::Dinic(Index count)
    : m_graph(count),
      m_level(count, -1),
      m_next(count)
  {}

  void MinSTCut<LocalMesh>::Dinic::add(Index first, Index second, Real capacity, Type type)
  {
    assert(first < m_graph.size() && second < m_graph.size());
    assert(std::isfinite(capacity) && capacity >= 0);
    Arc forward{second, m_graph[second].size(), capacity};
    Arc reverse{first, m_graph[first].size(), 0};
    m_graph[first].push_back(forward);
    m_graph[second].push_back(reverse);
    if (type == Type::Undirected)
      add(second, first, capacity, Type::Directed);
  }

  Real MinSTCut<LocalMesh>::Dinic::getMaximumFlow(Index source, Index sink)
  {
    assert(source < m_graph.size() && sink < m_graph.size() && source != sink);
    Real flow = 0;
    while (build(source, sink))
    {
      std::fill(m_next.begin(), m_next.end(), 0);
      while (true)
      {
        const Real pushed = augment(source, sink, std::numeric_limits<Real>::infinity());
        if (pushed <= 0)
          break;
        flow += pushed;
      }
    }
    return flow;
  }

  Boolean MinSTCut<LocalMesh>::Dinic::isReachable(Index vertex) const
  {
    assert(vertex < m_level.size());
    return m_level[vertex] >= 0;
  }

  Boolean MinSTCut<LocalMesh>::Dinic::build(Index source, Index sink)
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
        if (arc.capacity > 0 && m_level[arc.to] < 0)
        {
          m_level[arc.to] = m_level[current] + 1;
          queue.push_back(arc.to);
        }
      }
    }
    return m_level[sink] >= 0;
  }

  Real MinSTCut<LocalMesh>::Dinic::augment(Index current, Index sink, Real amount)
  {
    if (current == sink)
      return amount;
    for (Index& i = m_next[current]; i < m_graph[current].size(); ++i)
    {
      Arc& arc = m_graph[current][i];
      if (arc.capacity <= 0 || m_level[arc.to] != m_level[current] + 1)
        continue;
      const Real pushed = augment(arc.to, sink, std::min(amount, arc.capacity));
      if (pushed > 0)
      {
        arc.capacity -= pushed;
        m_graph[arc.to][arc.reverse].capacity += pushed;
        return pushed;
      }
    }
    return 0;
  }

  MinSTCut<LocalMesh>::MinSTCut(const MeshType& mesh)
    : m_mesh(mesh)
  {}

  MinSTCut<LocalMesh>::Result MinSTCut<LocalMesh>::classify(
    const std::function<Real(const Polytope&)>& average) const
  {
    const auto& mesh = getMesh();
    const auto& parameters = getParameters();
    const size_t d = mesh.getDimension();
    const size_t cellCount = mesh.getCellCount();
    if (d == 0 || cellCount == 0)
      return {};
    if (!average)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Cell average must be callable." << Alert::Raise;
    }
    if (!std::isfinite(parameters.fidelity) || parameters.fidelity < 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Fidelity must be finite and nonnegative." << Alert::Raise;
    }
    if (!parameters.smoothing)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "Smoothing must be callable." << Alert::Raise;
    }
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d - 1, d);
    const auto& incidence = mesh.getConnectivity().getIncidence(d - 1, d);
    std::vector<Real> insideCosts(cellCount);
    std::vector<Real> outsideCosts(cellCount);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Index i = cell->getIndex();
      const Real volume = cell->getMeasure();
      const Real value = average(*cell);
      if (!std::isfinite(volume) || volume <= 0 || !std::isfinite(value))
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Positive finite cell measures and finite averages are required."
          << Alert::Raise;
      }
      insideCosts[i] = parameters.fidelity * volume * std::max(Real(0), value);
      outsideCosts[i] = parameters.fidelity * volume * std::max(Real(0), -value);
      if (!std::isfinite(insideCosts[i]) || !std::isfinite(outsideCosts[i]))
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Nonfinite unary cost." << Alert::Raise;
      }
    }

    std::vector<Edge> edges;
    edges.reserve(mesh.getFaceCount());
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto& cells = incidence.at(face->getIndex());
      if (cells.size() != 1 && cells.size() != 2)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "One or two cells must be incident to each facet."
          << Alert::Raise;
      }
      if (cells.size() == 1)
        continue;
      const Real smoothing = parameters.smoothing(*face);
      if (!std::isfinite(smoothing) || smoothing < 0)
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Smoothing must be finite and nonnegative." << Alert::Raise;
      }
      // Codimension-one measure in 1D is counting measure, not vertex measure.
      const Real measure = d == 1 ? Real(1) : face->getMeasure();
      const Real capacity = measure * smoothing;
      if (!std::isfinite(measure) || measure <= 0 || !std::isfinite(capacity))
      {
        Alert::MemberFunctionException(*this, __func__)
          << "Invalid facet measure or capacity."
          << Alert::Raise;
      }
      edges.push_back({cells[0], cells[1], capacity, face->getIndex()});
    }

    return solve(insideCosts, outsideCosts, edges);
  }

  MinSTCut<LocalMesh>::Result MinSTCut<LocalMesh>::solve(const std::vector<Real>& insideCosts,
    const std::vector<Real>& outsideCosts, const std::vector<Edge>& edges) const
  {
    assert(insideCosts.size() == outsideCosts.size());

    const Index cellCount = insideCosts.size();
    const Index source = cellCount;
    const Index sink = cellCount + 1;

    Dinic graph(cellCount + 2);
    for (Index i = 0; i < cellCount; ++i)
    {
      assert(insideCosts[i] >= 0 && outsideCosts[i] >= 0);

      graph.add(source, i, outsideCosts[i], Dinic::Type::Directed);
      graph.add(i, sink, insideCosts[i], Dinic::Type::Directed);
    }

    for (const Edge& edge : edges)
    {
      assert(edge.first < cellCount && edge.second < cellCount);
      graph.add(edge.first, edge.second, edge.capacity, Dinic::Type::Undirected);
    }

    graph.getMaximumFlow(source, sink);

    Result result;
    for (Index i = 0; i < cellCount; ++i)
    {
      if (graph.isReachable(i))
      {
        result.inside.push_back(i);
        result.energy += insideCosts[i];
      }
      else
      {
        result.outside.push_back(i);
        result.energy += outsideCosts[i];
      }
    }

    for (const Edge& edge : edges)
    {
      if (graph.isReachable(edge.first) != graph.isReachable(edge.second))
      {
        result.cut.push_back(edge.face);
        result.energy += edge.capacity;
      }
    }
    return result;
  }
}
