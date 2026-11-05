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
  struct SewedOutput::StateFields
  {
    explicit StateFields(const Mesh& mesh)
      : velocitySpace(mesh, 3),
        pressureSpace(mesh),
        motion(velocitySpace)
    {
      for (size_t load = 0; load < 3; ++load)
      {
        translations[load] = std::make_unique<Velocity>(velocitySpace);
        rotations[load] = std::make_unique<Velocity>(velocitySpace);
        translationPressures[load] = std::make_unique<Pressure>(pressureSpace);
        rotationPressures[load] = std::make_unique<Pressure>(pressureSpace);
      }
    }

    VelocitySpace velocitySpace;
    PressureSpace pressureSpace;
    std::array<std::unique_ptr<Velocity>, 3> translations, rotations;
    std::array<std::unique_ptr<Pressure>, 3> translationPressures, rotationPressures;
    Velocity motion;
  };

  SewedOutput::~SewedOutput() = default;

  void SewedOutput::addStates(IO::XDMF::Grid& grid,
    const std::array<Velocity*, 3>& translations,
    const std::array<Velocity*, 3>& rotations,
    const std::array<Pressure*, 3>& translationPressures,
    const std::array<Pressure*, 3>& rotationPressures,
    const Math::SpatialVector<Real>& translation,
    const Math::SpatialVector<Real>& angular)
  {
    grid.clear();
    m_states = std::make_unique<StateFields>(m_mesh);
    auto& fields = *m_states;
    grid.setMesh(m_mesh, IO::XDMF::MeshPolicy::Transient);
    for (size_t load = 0; load < 3; ++load)
    {
      setVectorLoad(*fields.translations[load], translations, load);
      setVectorLoad(*fields.rotations[load], rotations, load);
      setScalarLoad(*fields.translationPressures[load], translationPressures, load);
      setScalarLoad(*fields.rotationPressures[load], rotationPressures, load);
    }
    // Preserve the field ordering of the chamber output.
    for (size_t load = 0; load < 3; ++load)
      grid.add("Translation_" + std::to_string(load), *fields.translations[load]);
    for (size_t load = 0; load < 3; ++load)
      grid.add("Rotation_" + std::to_string(load), *fields.rotations[load]);
    for (size_t load = 0; load < 3; ++load)
    {
      grid.add("Pressure_Translation_" + std::to_string(load),
        *fields.translationPressures[load]);
    }
    for (size_t load = 0; load < 3; ++load)
    {
      grid.add(
        "Pressure_Rotation_" + std::to_string(load), *fields.rotationPressures[load]);
    }
    fields.motion.getData() = translation(0) * fields.translations[0]->getData() +
      translation(1) * fields.translations[1]->getData() +
      translation(2) * fields.translations[2]->getData() +
      angular(0) * fields.rotations[0]->getData() +
      angular(1) * fields.rotations[1]->getData() +
      angular(2) * fields.rotations[2]->getData();
    grid.add("Motion", fields.motion);
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
    using Key = std::tuple<long long, long long, long long>;
    std::map<Key, Index> vertices;
    std::vector<Math::SpatialPoint> coordinates;
    std::vector<std::vector<Index>> copies(
      m_rotations.size(), std::vector<Index>(chamber.getVertexCount()));

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
          m_sources.emplace_back();
        }
        m_sources[position->second].push_back({r, vertex});
        copies[r][vertex] = position->second;
      }
    }

    Mesh::Builder builder;
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
