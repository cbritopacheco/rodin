/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "Rodin/Alert/MemberFunctionException.h"

#include "BalancedCompactPartitioner.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
  /// @brief Local mesh type partitioned by BalancedCompactPartitioner.
  using MeshType = Mesh<Context::Local>;

  BalancedCompactPartitioner::BalancedCompactPartitioner(const MeshType& mesh)
    : m_mesh(mesh)
  {
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, mesh.getDimension(), mesh.getDimension());
  }

  const Mesh<Context::Local>& BalancedCompactPartitioner::getMesh() const
  {
    return m_mesh.get();
  }

  void BalancedCompactPartitioner::partition(size_t count, size_t d)
  {
    if (count == 0)
    {
      Alert::MemberFunctionException(*this, __func__)
        << "The number of partitions must be positive." << Alert::Raise;
    }
    m_count = count;
    const MeshBase& mesh = getMesh();
    const auto& conn = mesh.getConnectivity();
    const size_t n = conn.getCount(d);
    m_partition.assign(n, 0);
    if (n == 0)
      return;
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, d, d);
    const auto& adjacency = conn.getIncidence(d, d);
    const size_t sdim = mesh.getSpaceDimension();

    // Centroids of the polytopes, stored once.
    std::vector<Real> centroids(n * sdim, 0);
    for (Index i = 0; i < n; i++)
    {
      const auto& vertices = conn.getPolytope(d, i);
      Real* c = centroids.data() + i * sdim;
      for (std::uint8_t v = 0; v < vertices.size(); v++)
      {
        const auto x = mesh.getVertexCoordinates(vertices[v]);
        for (size_t k = 0; k < sdim; k++)
          c[k] += x(k);
      }
      for (size_t k = 0; k < sdim; k++)
        c[k] /= vertices.size();
    }
    const auto distance = [&](Index i, const Real* c)
    {
      const Real* x = centroids.data() + i * sdim;
      Real r = 0;
      for (size_t k = 0; k < sdim; k++)
        r += (x[k] - c[k]) * (x[k] - c[k]);
      return r;
    };

    // Recursive bisection. A set of polytopes for the parts [first, first +
    // count) splits into one region for the first count / 2 parts and the rest
    // for the others. The region starts from a polytope at the periphery of
    // the set and grows through the adjacency, lowest first along the axis to
    // the polytope furthest from it, until it has its exact share. Every
    // polytope of a set carries the first part of the set, so that the growth
    // stays inside it.
    struct Task
    {
      std::vector<Index> cells;
      size_t first;
      size_t count;
    };
    std::vector<Task> tasks;
    {
      std::vector<Index> cells(n);
      for (Index i = 0; i < n; i++)
        cells[i] = i;
      tasks.push_back({ std::move(cells), 0, count });
    }
    // The height of a polytope does not change, so a polytope enters the front
    // once.
    std::vector<bool> queued(n, false);
    using Candidate = std::pair<Real, Index>;
    std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> front;
    while (!tasks.empty())
    {
      Task task = std::move(tasks.back());
      tasks.pop_back();
      if (task.count == 1 || task.cells.empty())
        continue;
      const size_t half = task.count / 2;
      const size_t second = task.first + half;
      const size_t target = task.cells.size() * half / task.count;
      for (const Index i : task.cells)
        m_partition[i] = second;

      // The polytope furthest from the first one of the set.
      const auto furthest = [&](Index from)
      {
        const Real* c = centroids.data() + from * sdim;
        Index best = from;
        Real bestDistance = -1;
        for (const Index i : task.cells)
        {
          if (m_partition[i] != second)
            continue;
          const Real r = distance(i, c);
          if (r > bestDistance)
          {
            bestDistance = r;
            best = i;
          }
        }
        return best;
      };
      const Index seed = furthest(task.cells.front());
      const Index end = furthest(seed);
      const Real* origin = centroids.data() + seed * sdim;
      std::vector<Real> axis(sdim);
      for (size_t k = 0; k < sdim; k++)
        axis[k] = centroids[end * sdim + k] - origin[k];
      const auto height = [&](Index i)
      {
        const Real* x = centroids.data() + i * sdim;
        Real r = 0;
        for (size_t k = 0; k < sdim; k++)
          r += (x[k] - origin[k]) * axis[k];
        return r;
      };

      size_t taken = 0;
      while (taken < target)
      {
        if (front.empty())
        {
          // A new component of the set, from its lowest polytope.
          Index next = n;
          Real lowest = std::numeric_limits<Real>::max();
          for (const Index i : task.cells)
          {
            if (m_partition[i] != second)
              continue;
            const Real h = height(i);
            if (h < lowest)
            {
              lowest = h;
              next = i;
            }
          }
          front.push({ lowest, next });
          queued[next] = true;
        }
        const Index i = front.top().second;
        front.pop();
        m_partition[i] = task.first;
        taken++;
        for (const Index j : adjacency[i])
        {
          if (m_partition[j] == second && !queued[j])
          {
            front.push({ height(j), j });
            queued[j] = true;
          }
        }
      }
      front = {};
      for (const Index i : task.cells)
        queued[i] = false;

      std::vector<Index> region, rest;
      region.reserve(target);
      rest.reserve(task.cells.size() - target);
      for (const Index i : task.cells)
        (m_partition[i] == task.first ? region : rest).push_back(i);
      task.cells = {};
      tasks.push_back({ std::move(rest), second, task.count - half });
      tasks.push_back({ std::move(region), task.first, half });
    }
  }

  size_t BalancedCompactPartitioner::getPartition(Index index) const
  {
    return m_partition[index];
  }

  size_t BalancedCompactPartitioner::getCount() const
  {
    return m_count;
  }
}
