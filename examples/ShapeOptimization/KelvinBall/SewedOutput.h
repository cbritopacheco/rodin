/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_SEWED_OUTPUT_H
#define KELVIN_BALL_SEWED_OUTPUT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

namespace KelvinBall
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  inline std::vector<Math::SpatialMatrix<Real>> properOctahedralRotations()
  {
    std::vector<Math::SpatialMatrix<Real>> rotations;
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
            if (std::abs(rotation.determinant() - 1.0) < 1e-12)
              rotations.push_back(rotation);
          }
        }
      }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return rotations;
  }

  struct SewedMesh
  {
      struct Source
      {
          size_t rotation;
          Index vertex;
      };

      Geometry::Mesh<Context::Local> mesh;
      std::vector<std::vector<Source>> sources;
      std::vector<Math::SpatialMatrix<Real>> rotations;
  };

  template <class Mesh>
  SewedMesh sew(const Mesh& chamber, const FlatSet<Attribute>& boundaryAttributes = {},
    Real tolerance = 1e-10)
  {
    using Key = std::tuple<long long, long long, long long>;
    SewedMesh result;
    result.rotations = properOctahedralRotations();
    std::map<Key, Index> vertices;
    std::vector<Math::SpatialPoint> coordinates;
    std::vector<std::vector<Index>> copies(
      result.rotations.size(), std::vector<Index>(chamber.getVertexCount()));

    for (size_t r = 0; r < result.rotations.size(); ++r)
    {
      for (Index vertex = 0; vertex < chamber.getVertexCount(); ++vertex)
      {
        const Math::SpatialPoint x =
          result.rotations[r] * chamber.getVertexCoordinates(vertex);
        const Key key{std::llround(x(0) / tolerance), std::llround(x(1) / tolerance),
          std::llround(x(2) / tolerance)};
        const auto [position, inserted] =
          vertices.emplace(key, static_cast<Index>(coordinates.size()));
        if (inserted)
        {
          coordinates.push_back(x);
          result.sources.emplace_back();
        }
        result.sources[position->second].push_back({r, vertex});
        copies[r][vertex] = position->second;
      }
    }

    Geometry::Mesh<Context::Local>::Builder builder;
    builder.initialize(3).nodes(coordinates.size());
    for (const auto& x : coordinates)
      builder.vertex(x);

    for (size_t r = 0; r < result.rotations.size(); ++r)
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
    result.mesh = builder.finalize();
    return result;
  }

  template <class Output, class Input>
  void sewScalar(Output& output, const Input& input, const SewedMesh& sewed)
  {
    for (Index vertex = 0; vertex < sewed.mesh.getVertexCount(); ++vertex)
    {
      const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
      Real value = 0;
      for (const auto& source : sewed.sources[vertex])
      {
        const auto sourceDOFs = input.getFiniteElementSpace().getDOFs(0, source.vertex);
        value += input.getData()(sourceDOFs(0));
      }
      output.getData()(targetDOFs(0)) = value / sewed.sources[vertex].size();
    }
  }

  template <class Output, class Input>
  void sewVector(Output& output, const Input& input, const SewedMesh& sewed)
  {
    for (Index vertex = 0; vertex < sewed.mesh.getVertexCount(); ++vertex)
    {
      const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
      Math::SpatialVector<Real> value(3);
      value.setZero();
      for (const auto& source : sewed.sources[vertex])
      {
        const auto sourceDOFs = input.getFiniteElementSpace().getDOFs(0, source.vertex);
        Math::SpatialVector<Real> chamberValue(3);
        for (size_t component = 0; component < 3; ++component)
          chamberValue(component) = input.getData()(sourceDOFs(component));
        value += sewed.rotations[source.rotation] * chamberValue;
      }
      value /= sewed.sources[vertex].size();
      for (size_t component = 0; component < 3; ++component)
        output.getData()(targetDOFs(component)) = value(component);
    }
  }

  template <class Output, class Family>
  void sewVectorLoad(
    Output& output, const Family& family, size_t load, const SewedMesh& sewed)
  {
    for (Index vertex = 0; vertex < sewed.mesh.getVertexCount(); ++vertex)
    {
      const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
      Math::SpatialVector<Real> value(3);
      value.setZero();
      for (const auto& source : sewed.sources[vertex])
      {
        const auto& rotation = sewed.rotations[source.rotation];
        for (size_t sourceLoad = 0; sourceLoad < 3; ++sourceLoad)
        {
          const auto sourceDOFs =
            family[sourceLoad]->getFiniteElementSpace().getDOFs(0, source.vertex);
          Math::SpatialVector<Real> chamberValue(3);
          for (size_t component = 0; component < 3; ++component)
          {
            chamberValue(component) =
              family[sourceLoad]->getData()(sourceDOFs(component));
          }
          value += rotation(load, sourceLoad) * rotation * chamberValue;
        }
      }
      value /= sewed.sources[vertex].size();
      for (size_t component = 0; component < 3; ++component)
        output.getData()(targetDOFs(component)) = value(component);
    }
  }

  template <class Output, class Family>
  void sewScalarLoad(
    Output& output, const Family& family, size_t load, const SewedMesh& sewed)
  {
    for (Index vertex = 0; vertex < sewed.mesh.getVertexCount(); ++vertex)
    {
      const auto targetDOFs = output.getFiniteElementSpace().getDOFs(0, vertex);
      Real value = 0;
      for (const auto& source : sewed.sources[vertex])
      {
        const auto& rotation = sewed.rotations[source.rotation];
        for (size_t sourceLoad = 0; sourceLoad < 3; ++sourceLoad)
        {
          const auto sourceDOFs =
            family[sourceLoad]->getFiniteElementSpace().getDOFs(0, source.vertex);
          value +=
            rotation(load, sourceLoad) * family[sourceLoad]->getData()(sourceDOFs(0));
        }
      }
      output.getData()(targetDOFs(0)) = value / sewed.sources[vertex].size();
    }
  }
}

#endif
