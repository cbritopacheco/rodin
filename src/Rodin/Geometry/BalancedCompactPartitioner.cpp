#include "BalancedCompactPartitioner.h"
#include "Mesh.h"
#include "Polytope.h"
#include "PolytopeIterator.h"
#include <queue>
#include <vector>
#include <cmath>
#include <limits>

namespace Rodin::Geometry
{
  using MeshType = Mesh<Context::Local>;

  // Helper: Compute the centroid of a polytope (cell) of dimension d.
  // The centroid is computed as the average of the vertex coordinates.
  static Math::SpatialVector<Real> computePolytopeCentroid(const MeshBase& mesh, size_t d, Index polyIndex)
  {
    auto poly = mesh.getPolytope(d, polyIndex);
    const auto& vertices = poly->getVertices();
    size_t spaceDim = mesh.getSpaceDimension();
    Math::SpatialVector<Real> centroid(spaceDim);
    centroid.setZero();
    for (size_t i = 0; i < static_cast<size_t>(vertices.size()); ++i)
    {
      centroid += mesh.getVertexCoordinates(vertices[i]);
    }
    centroid /= vertices.size();
    return centroid;
  }

  // Structure to hold cluster (seed) information.
  struct ClusterInfo
  {
    size_t size;                        // Number of cells already assigned.
    Math::SpatialVector<Real> centroid; // Running centroid of the cluster.
    Real radius;                        // Maximum distance of any assigned cell from centroid.
  };

  BalancedCompactPartitioner::BalancedCompactPartitioner(const MeshType& mesh)
    : m_mesh(mesh)
  {}

  const Mesh<Context::Local>& BalancedCompactPartitioner::getMesh() const
  {
    return m_mesh.get();
  }

  // Partition the mesh into exactly numClusters clusters along dimension d.
  // The algorithm selects numClusters seeds (via furthest–point sampling),
  // computes the ideal cluster size, and then grows regions via a best–first search.
  // To enforce compactness, a candidate cell is rejected if its distance from
  // the cluster's current centroid is more than compactFactor times the cluster's radius.
  void BalancedCompactPartitioner::partition(size_t numClusters, size_t d)
  {
    m_count = numClusters;
    const MeshBase& mesh = getMesh();
    size_t n = mesh.getPolytopeCount(d);
    if (n == 0)
      return;

    // --- Seed Selection: Furthest–Point Sampling ---
    struct SeedItem {
      Index polyIndex;
      Math::SpatialVector<Real> centroid;
    };
    std::vector<SeedItem> seeds;
    std::vector<bool> seedChosen(n, false);

    // Choose the first seed arbitrarily (cell 0).
    seeds.push_back({0, computePolytopeCentroid(mesh, d, 0)});
    seedChosen[0] = true;

    for (size_t s = 1; s < numClusters; s++)
    {
      Index bestCandidate = 0;
      Real bestDist = -1;
      for (Index i = 0; i < n; i++)
      {
        if (seedChosen[i])
          continue;
        Math::SpatialVector<Real> cent = computePolytopeCentroid(mesh, d, i);
        Real minDist = std::numeric_limits<Real>::max();
        for (const auto& seed : seeds)
        {
          Real dVal = (cent - seed.centroid).norm();
          if (dVal < minDist)
            minDist = dVal;
        }
        if (minDist > bestDist)
        {
          bestDist = minDist;
          bestCandidate = i;
        }
      }
      seeds.push_back({bestCandidate, computePolytopeCentroid(mesh, d, bestCandidate)});
      seedChosen[bestCandidate] = true;
    }

    // --- Region Growing with Compactness Constraint ---
    // Initialize the partition mapping.
    m_partition.resize(n, n); // Default invalid id = n.
    // Cluster info: for each seed, maintain its size, current centroid, and current radius.
    std::vector<ClusterInfo> clusters(numClusters);
    for (size_t s = 0; s < numClusters; s++)
    {
      clusters[s].size = 0;
      clusters[s].centroid = seeds[s].centroid;
      clusters[s].radius = 0; // Initially, no cell has been added beyond the seed.
    }
    std::vector<bool> assigned(n, false);

    // Define the compactness factor. Only cells within compactFactor * (current cluster radius)
    // from the cluster centroid are accepted. (For a seed with zero radius, we always accept.)
    const Real compactFactor = 2.0;

    // PQ item holds a candidate cell assignment.
    struct PQItem {
      Index polyIndex;
      size_t seedId;
      Real distance; // Distance from the cluster's current centroid.
    };
    struct ComparePQItem {
      bool operator()(const PQItem& a, const PQItem& b)
      {
        return a.distance > b.distance;
      }
    };
    std::priority_queue<PQItem, std::vector<PQItem>, ComparePQItem> pq;

    // Initialize PQ with the seed cells.
    for (size_t s = 0; s < seeds.size(); s++)
    {
      pq.push({seeds[s].polyIndex, s, 0.0});
    }

    // Grow regions until all cells in n are assigned.
    while (!pq.empty())
    {
      PQItem item = pq.top();
      pq.pop();
      Index cell = item.polyIndex;
      size_t seedId = item.seedId;
      if (assigned[cell])
        continue;

      // Check compactness: if the cluster is non-empty, reject cell if too far.
      Math::SpatialVector<Real> cellCentroid = computePolytopeCentroid(mesh, d, cell);
      Real dist = (cellCentroid - clusters[seedId].centroid).norm();
      if (clusters[seedId].size > 0 && clusters[seedId].radius > 0)
      {
        if (dist > compactFactor * clusters[seedId].radius)
          continue;
      }

      // Accept the cell into the cluster.
      m_partition[cell] = seedId;
      assigned[cell] = true;

      // Update cluster info.
      size_t oldSize = clusters[seedId].size;
      Math::SpatialVector<Real> oldCentroid = clusters[seedId].centroid;
      clusters[seedId].size++;
      // New centroid is the weighted average.
      clusters[seedId].centroid = (oldCentroid * oldSize + cellCentroid) / clusters[seedId].size;
      // Update radius: use the maximum distance from the new centroid.
      Real newDist = (cellCentroid - clusters[seedId].centroid).norm();
      clusters[seedId].radius = std::max(clusters[seedId].radius, newDist);

      // Expand: push all unassigned adjacent cells.
      auto poly = mesh.getPolytope(d, cell);
      for (auto adj = poly->getAdjacent(); !adj.end(); ++adj)
      {
        Index adjIndex = (*adj).getIndex();
        if (assigned[adjIndex])
          continue;
        Math::SpatialVector<Real> adjCentroid = computePolytopeCentroid(mesh, d, adjIndex);
        Real newDistance = (adjCentroid - clusters[seedId].centroid).norm();
        pq.push({adjIndex, seedId, newDistance});
      }
    }

    // For any cell still unassigned, assign it to the nearest seed.
    for (Index i = 0; i < n; i++)
    {
      if (!assigned[i])
      {
        Math::SpatialVector<Real> cent = computePolytopeCentroid(mesh, d, i);
        size_t bestSeed = 0;
        Real bestD = std::numeric_limits<Real>::max();
        for (size_t s = 0; s < numClusters; s++)
        {
          Real dVal = (cent - clusters[s].centroid).norm();
          if (dVal < bestD)
          {
            bestD = dVal;
            bestSeed = s;
          }
        }
        m_partition[i] = bestSeed;
        clusters[bestSeed].size++;
        assigned[i] = true;
      }
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
