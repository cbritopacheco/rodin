/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "LaplacianSmoothing.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace KelvinBall
{
  size_t LaplacianSmoothing::execute(Mesh& mesh) const
  {
    // Conservative heuristic budget: two half-relaxed Jacobi passes, with at
    // most eight halvings to preserve validity and the existing quality floor.
    constexpr size_t passes = 2, backtracks = 8;
    constexpr Real relaxation = Real(0.5);
    constexpr Real orthogonalityTolerance = 64 * std::numeric_limits<Real>::epsilon();
    const size_t count = mesh.getVertexCount();
    std::vector<FlatSet<Index>> neighbours(count);
    std::vector<bool> fixed(count, false);
    std::vector<std::vector<Math::SpatialVector<Real>>> normals(count);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& vertices = cell->getVertices();
      if (vertices.size() != 4)
        throw std::runtime_error("KelvinBall Laplacian smoothing requires tetrahedra.");
      for (const Index i : vertices)
      {
        for (const Index j : vertices)
        {
          if (i != j)
            neighbours[i].insert(j);
        }
      }
    }
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (attribute == Gamma || attribute == Outer)
      {
        for (const Index vertex : face->getVertices())
          fixed[vertex] = true;
      }
      Math::SpatialVector<Real> normal(3);
      normal.setZero();
      if (attribute == SigmaPlus)
      {
        normal(1) = 1;
        normal(2) = -1;
      }
      else if (attribute == SigmaMinus)
      {
        normal(1) = 1;
        normal(2) = 1;
      }
      else if (attribute == SigmaXYPlus || attribute == SigmaXYMinus)
      {
        normal(0) = 1;
        normal(1) = -1;
      }
      else
        continue;
      for (const Index vertex : face->getVertices())
      {
        auto independent = normal;
        for (const auto& basis : normals[vertex])
          independent -= independent.dot(basis) * basis;
        if (independent.norm() > orthogonalityTolerance)
          normals[vertex].push_back(independent.normalized());
      }
    }
    // Native point transformations supply signed Jacobians, including cache
    // invalidation after coordinate changes. The quality matches the CSV metric.
    const auto qualities = [](const Mesh& input) {
      std::vector<std::pair<Real, Real>> values;
      values.reserve(input.getCellCount());
      for (auto cell = input.getCell(); cell; ++cell)
      {
        Math::SpatialPoint reference(3);
        reference.setConstant(Real(0.25));
        Geometry::Point point(*cell, reference);
        const Real determinant = point.getJacobian().determinant();
        Real squaredEdges = 0;
        const auto& vertices = cell->getVertices();
        for (size_t i = 0; i < 4; ++i)
        {
          for (size_t j = i + 1; j < 4; ++j)
          {
            squaredEdges += (input.getVertexCoordinates(vertices[i]) -
              input.getVertexCoordinates(vertices[j])).squaredNorm();
          }
        }
        const Real quality = 12 * std::pow(std::abs(determinant) / 2,
          Real(2) / 3) / squaredEdges;
        values.emplace_back(determinant, quality);
      }
      return values;
    };
    size_t accepted = 0;
    for (size_t pass = 0; pass < passes; ++pass)
    {
      std::vector<Math::SpatialPoint> coordinates(count), displacement(count);
      bool moving = false;
      for (Index i = 0; i < count; ++i)
      {
        coordinates[i] = mesh.getVertexCoordinates(i);
        displacement[i] = Math::SpatialPoint::Zero(3);
        if (fixed[i] || neighbours[i].empty())
          continue;
        for (const Index j : neighbours[i])
          displacement[i] += mesh.getVertexCoordinates(j);
        displacement[i] /= static_cast<Real>(neighbours[i].size());
        displacement[i] -= coordinates[i];
        for (const auto& normal : normals[i])
          displacement[i] -= displacement[i].dot(normal) * normal;
        moving = moving || displacement[i].squaredNorm() > 0;
      }
      if (!moving)
        break;
      const auto before = qualities(mesh);
      Real minimum = 1, sum = 0;
      for (const auto& value : before)
      {
        minimum = std::min(minimum, value.second);
        sum += value.second;
      }
      Real alpha = relaxation;
      bool committed = false;
      for (size_t attempt = 0; attempt <= backtracks; ++attempt, alpha /= 2)
      {
        Mesh candidate(mesh);
        for (Index i = 0; i < count; ++i)
          candidate.setVertexCoordinates(i, coordinates[i] + alpha * displacement[i]);
        const auto after = qualities(candidate);
        bool valid = true;
        Real candidateMinimum = 1, candidateSum = 0;
        for (size_t i = 0; i < after.size(); ++i)
        {
          valid = valid && std::isfinite(after[i].first) &&
            std::isfinite(after[i].second) && after[i].first * before[i].first > 0;
          candidateMinimum = std::min(candidateMinimum, after[i].second);
          candidateSum += after[i].second;
        }
        if (valid && candidateMinimum >= minimum && candidateSum >= sum)
        {
          // Update coordinates only: preserve MMG feature tags on the caller.
          for (Index i = 0; i < count; ++i)
            mesh.setVertexCoordinates(i, candidate.getVertexCoordinates(i));
          ++accepted;
          committed = true;
          break;
        }
      }
      if (!committed)
        break;
    }
    return accepted;
  }
}
