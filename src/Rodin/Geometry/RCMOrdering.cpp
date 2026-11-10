/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <numeric>

#include "RCMOrdering.h"

namespace Rodin::Geometry
{
  IndexVector RCMOrdering::operator()(const Mesh<Context::Local>& mesh) const
  {
    const size_t count = mesh.getVertexCount();
    std::vector<IndexVector> graph(count);
    const auto& connectivity = mesh.getConnectivity();
    const size_t dim = mesh.getDimension();
    for (Index cell = 0; cell < mesh.getCellCount(); ++cell)
    {
      const auto& vertices = connectivity.getPolytope(dim, cell);
      for (size_t i = 0; i < vertices.size(); ++i)
      {
        for (size_t j = i + 1; j < vertices.size(); ++j)
        {
          if (vertices[i] != vertices[j])
          {
            graph[vertices[i]].push_back(vertices[j]);
            graph[vertices[j]].push_back(vertices[i]);
          }
        }
      }
    }
    for (auto& neighbors : graph)
    {
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
    const auto less = [&](Index a, Index b) {
      return graph[a].size() == graph[b].size() ? a < b
                                                : graph[a].size() < graph[b].size();
    };
    for (auto& neighbors : graph)
      std::sort(neighbors.begin(), neighbors.end(), less);
    IndexVector seeds(count);
    std::iota(seeds.begin(), seeds.end(), Index(0));
    std::sort(seeds.begin(), seeds.end(), less);
    IndexVector order;
    order.reserve(count);
    IndexVector queue;
    queue.reserve(count);
    std::vector<bool> visited(count, false), explored(count, false);
    for (const Index seed : seeds)
    {
      if (visited[seed])
        continue;
      Index root = seed;
      size_t previousDepth = 0;
      // Repeated level structures move the root toward the component periphery.
      while (true)
      {
        queue.clear();
        queue.push_back(root);
        explored[root] = true;
        size_t begin = 0, depth = 0;
        Index candidate = root;
        while (begin < queue.size())
        {
          const size_t end = queue.size();
          candidate = queue[begin];
          for (size_t i = begin; i < end; ++i)
          {
            const Index vertex = queue[i];
            if (less(vertex, candidate))
              candidate = vertex;
            for (const Index neighbor : graph[vertex])
            {
              if (!explored[neighbor])
              {
                explored[neighbor] = true;
                queue.push_back(neighbor);
              }
            }
          }
          begin = end;
          ++depth;
        }
        for (const Index vertex : queue)
          explored[vertex] = false;
        if (depth <= previousDepth || candidate == root)
          break;
        previousDepth = depth;
        root = candidate;
      }
      queue.clear();
      queue.push_back(root);
      visited[root] = true;
      for (size_t i = 0; i < queue.size(); ++i)
      {
        for (const Index neighbor : graph[queue[i]])
        {
          if (!visited[neighbor])
          {
            visited[neighbor] = true;
            queue.push_back(neighbor);
          }
        }
      }
      order.insert(order.end(), queue.rbegin(), queue.rend());
    }
    return order;
  }
}
