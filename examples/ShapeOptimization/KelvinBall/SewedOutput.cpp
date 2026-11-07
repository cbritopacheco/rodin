/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "SewedOutput.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <tuple>

namespace KelvinBall
{
  void SewedOutput::addStates(IO::XDMF::Grid& grid,
    const std::array<Velocity*, 3>& translations,
    const std::array<Velocity*, 3>& rotations,
    const std::array<Pressure*, 3>& translationPressures,
    const std::array<Pressure*, 3>& rotationPressures,
    const Math::SpatialVector<Real>& translation,
    const Math::SpatialVector<Real>& angular)
  {
    grid.clear();
    grid.setMesh(m_mesh, IO::XDMF::MeshPolicy::Transient);
    for (size_t load = 0; load < 3; ++load)
    {
      grid.add("Translation_" + std::to_string(load), 3,
        [this, translations, load](Index vertex, std::span<Real> value) {
          const auto sample = getVectorLoad(translations, load, vertex);
          for (size_t component = 0; component < 3; ++component)
            value[component] = sample(component);
        });
    }
    for (size_t load = 0; load < 3; ++load)
    {
      grid.add("Rotation_" + std::to_string(load), 3,
        [this, rotations, load](Index vertex, std::span<Real> value) {
          const auto sample = getVectorLoad(rotations, load, vertex);
          for (size_t component = 0; component < 3; ++component)
            value[component] = sample(component);
        });
    }
    for (size_t load = 0; load < 3; ++load)
    {
      grid.add("Pressure_Translation_" + std::to_string(load), 1,
        [this, translationPressures, load](Index vertex, std::span<Real> value) {
          value[0] = getScalarLoad(translationPressures, load, vertex);
        });
    }
    for (size_t load = 0; load < 3; ++load)
    {
      grid.add("Pressure_Rotation_" + std::to_string(load), 1,
        [this, rotationPressures, load](Index vertex, std::span<Real> value) {
          value[0] = getScalarLoad(rotationPressures, load, vertex);
        });
    }
    grid.add("Motion", 3,
      [this, translations, rotations, translation, angular](
        Index vertex, std::span<Real> value) {
        const auto t0 = getVectorLoad(translations, 0, vertex);
        const auto t1 = getVectorLoad(translations, 1, vertex);
        const auto t2 = getVectorLoad(translations, 2, vertex);
        const auto r0 = getVectorLoad(rotations, 0, vertex);
        const auto r1 = getVectorLoad(rotations, 1, vertex);
        const auto r2 = getVectorLoad(rotations, 2, vertex);
        for (size_t component = 0; component < 3; ++component)
        {
          // Preserve the separate products of the original Eigen field expression.
          const Real t0Value = translation(0) * t0(component);
          const Real t1Value = translation(1) * t1(component);
          const Real t2Value = translation(2) * t2(component);
          const Real r0Value = angular(0) * r0(component);
          const Real r1Value = angular(1) * r1(component);
          const Real r2Value = angular(2) * r2(component);
          value[component] = t0Value + t1Value + t2Value + r0Value + r1Value + r2Value;
        }
      });
  }

  const std::vector<Math::SpatialMatrix<Real>>& SewedOutput::getCubeRotations()
  {
    static const std::vector<Math::SpatialMatrix<Real>> rotations = [] {
      std::vector<Math::SpatialMatrix<Real>> result;
      std::array<size_t, 3> permutation{0, 1, 2};
      do
      {
        for (const int sx : {-1, 1})
        {
          for (const int sy : {-1, 1})
          {
            for (const int sz : {-1, 1})
            {
              Math::SpatialMatrix<Real> rotation(3, 3);
              rotation.setZero();
              rotation(0, permutation[0]) = sx;
              rotation(1, permutation[1]) = sy;
              rotation(2, permutation[2]) = sz;
              if (std::abs(rotation.determinant() - 1) < 1e-12)
                result.push_back(rotation);
            }
          }
        }
      } while (std::next_permutation(permutation.begin(), permutation.end()));
      return result;
    }();
    return rotations;
  }

  SewedOutput::SewedOutput(
    const Mesh& chamber, const FlatSet<Attribute>& boundaryAttributes, Real tolerance)
    : m_rotations(getCubeRotations())
  {
    Mesh::Builder builder;
    {
      using Key = std::tuple<long long, long long, long long>;
      std::vector<Math::SpatialPoint> coordinates;
    std::vector<std::vector<Index>> copies(
      m_rotations.size(), std::vector<Index>(chamber.getVertexCount()));
      std::vector<Index> counts;
    {
        std::map<Key, Index> vertices;
        for (size_t r = 0; r < m_rotations.size(); ++r)
      {
          for (Index vertex = 0; vertex < chamber.getVertexCount(); ++vertex)
        {
            const Math::SpatialPoint x =
              m_rotations[r] * chamber.getVertexCoordinates(vertex);
            const Key key{std::llround(x(0) / tolerance), std::llround(x(1) / tolerance),
              std::llround(x(2) / tolerance)};
            const auto [position, inserted] =
              vertices.emplace(key, static_cast<Index>(coordinates.size()));
            if (inserted)
            {
              coordinates.push_back(x);
              counts.push_back(0);
            }
            ++counts[position->second];
            copies[r][vertex] = position->second;
          }
        }
      }
      // Flat representatives avoid an allocation and vector header per sewn node.
      m_sourceOffsets.resize(counts.size() + 1, 0);
      for (Index vertex = 0; vertex < counts.size(); ++vertex)
        m_sourceOffsets[vertex + 1] = m_sourceOffsets[vertex] + counts[vertex];
      m_sources.resize(m_sourceOffsets.back());
      std::fill(counts.begin(), counts.end(), 0);
      for (size_t r = 0; r < m_rotations.size(); ++r)
      {
        for (Index vertex = 0; vertex < chamber.getVertexCount(); ++vertex)
        {
          const Index target = copies[r][vertex];
          m_sources[m_sourceOffsets[target] + counts[target]++] = {r, vertex};
        }
      }
      builder.initialize(3).nodes(coordinates.size());
    for (const auto& x : coordinates)
      builder.vertex(x);

    for (size_t r = 0; r < m_rotations.size(); ++r)
    {
      for (auto cell = chamber.getCell(); cell; ++cell)
      {
        IndexArray vertices(cell->getVertices().size());
        for (size_t local = 0; local < vertices.size(); ++local)
          vertices(local) = copies[r][cell->getVertices()(local)];
        Index index;
        builder.polytope(cell->getGeometry(), std::move(vertices), index);
        builder.attribute({chamber.getDimension(), index}, cell->getAttribute());
      }
      for (auto face = chamber.getPolytope(chamber.getDimension() - 1); face; ++face)
      {
        const auto attribute = face->getAttribute();
        if (!attribute || !boundaryAttributes.contains(*attribute))
          continue;
        IndexArray vertices(face->getVertices().size());
        for (size_t local = 0; local < vertices.size(); ++local)
          vertices(local) = copies[r][face->getVertices()(local)];
        Index index;
        builder.polytope(face->getGeometry(), std::move(vertices), index);
        builder.attribute({chamber.getDimension() - 1, index}, attribute);
      }
    }
    }
    m_mesh = builder.finalize();
  }

  const SewedOutput::Mesh& SewedOutput::getMesh() const
  {
    return m_mesh;
  }

  SewedOutput::Mesh& SewedOutput::getMesh()
  {
    return m_mesh;
  }

}
