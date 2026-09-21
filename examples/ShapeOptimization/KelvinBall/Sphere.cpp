/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Sphere.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <Rodin/Variational.h>

namespace KelvinBall
{
  Sphere::Sphere(const Configuration& configuration)
    : m_configuration(configuration)
  {}

  Mesh Sphere::makeUniformChamber() const
  {
    Mesh cube = Mesh::UniformGrid(Polytope::Type::Tetrahedron,
      {m_configuration.points, m_configuration.points, m_configuration.points});
    cube.scale(m_configuration.getH());

    using Key = std::tuple<long long, long long, long long>;
    constexpr Real tolerance = 1e-12;
    std::map<Key, Index> indices;
    std::vector<Math::SpatialPoint> coordinates;
    std::vector<std::pair<IndexArray, Attribute>> cells;
    const auto insertVertex = [&](const Math::SpatialPoint& x) {
      const Key key{std::llround(x(0) / tolerance), std::llround(x(1) / tolerance),
        std::llround(x(2) / tolerance)};
      const auto [position, inserted] =
        indices.emplace(key, static_cast<Index>(coordinates.size()));
      if (inserted)
        coordinates.push_back(x);
      return position->second;
    };

    for (auto cell = cube.getCell(); cell; ++cell)
    {
      bool inside = true;
      for (const Index vertex : cell->getVertices())
      {
        const auto& x = cube.getVertexCoordinates(vertex);
        inside = inside && x(0) + tolerance >= x(1) && x(1) + tolerance >= x(2);
      }
      if (!inside)
        continue;
      for (const Real sign : {Real(1), Real(-1)})
      {
        IndexArray vertices(cell->getVertices().size());
        for (size_t local = 0; local < vertices.size(); ++local)
        {
          Math::SpatialPoint x = cube.getVertexCoordinates(cell->getVertices()(local));
          x(2) *= sign;
          vertices(local) = insertVertex(x);
        }
        if (sign < 0)
          std::swap(vertices(0), vertices(1));
        cells.emplace_back(std::move(vertices), Fluid);
      }
    }

    Mesh::Builder builder;
    builder.initialize(3).nodes(coordinates.size());
    for (const auto& x : coordinates)
      builder.vertex(x);
    for (auto& [vertices, attribute] : cells)
    {
      Index index;
      builder.polytope(Polytope::Type::Tetrahedron, std::move(vertices), index);
      builder.attribute({3, index}, attribute);
    }
    Mesh chamber = builder.finalize();
    chamber.getConnectivity().compute(2, 3);
    for (auto face = chamber.getBoundary(); face; ++face)
    {
      const auto c = centroid(chamber, *face);
      Attribute attribute = 0;
      if (std::abs(c(0) - m_configuration.outerRadius) < tolerance)
        attribute = Outer;
      else if (std::abs(c(0) - c(1)) < tolerance)
        attribute = c(2) < 0 ? SigmaXYMinus : SigmaXYPlus;
      else if (std::abs(c(1) - c(2)) < tolerance)
        attribute = SigmaPlus;
      else if (std::abs(c(1) + c(2)) < tolerance)
        attribute = SigmaMinus;
      if (attribute == 0)
        throw std::runtime_error("The uniform chamber has an unclassified boundary.");
      chamber.setAttribute({2, face->getIndex()}, attribute);
    }
    return chamber;
  }

  size_t Sphere::protectFixedGeometry(MMG::Mesh& mesh) const
  {
    mesh.getRequiredTriangles().clear();
    size_t count = 0;
    const FlatSet<Attribute> fixed{
      Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
    for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (attribute && fixed.contains(*attribute))
      {
        mesh.setRequiredTriangle(face->getIndex());
        ++count;
      }
    }
    return count;
  }

  SphereDiscretization Sphere::discretize() const
  {
    const Real h = m_configuration.getH();
    MMG::Mesh mesh(makeUniformChamber());
    const size_t cellsBefore = mesh.getCellCount();
    protectFixedGeometry(mesh);
    const Real hmin = 0.1 * h;
    const Real hmax = 10 * h;
    const Real hausdorff = 0.1 * h * h;
    P1 levelSetSpace(mesh);
    GridFunction sphere(levelSetSpace);
    sphere = RealFunction([](const Geometry::Point& point) {
      return point.getPhysicalCoordinates().norm() - Real(1);
    });
    MMG::LevelSetDiscretizer discretizer;
    discretizer.split(Fluid, {Obstacle, Fluid})
      .setHMin(hmin)
      .setHMax(hmax)
      .setHausdorff(hausdorff)
      .setGradation(2)
      .setBaseReferences(FlatSet<Attribute>{Fluid})
      .setBoundaryReference(Gamma)
      .setAngleDetection(false);
    mesh = discretizer.discretize(sphere);
    splitSelfPairedCut(mesh);
    protectFixedGeometry(mesh);
    MMG::Optimizer()
      .setHMin(hmin)
      .setHMax(hmax)
      .setHausdorff(hausdorff)
      .setGradation(2)
      .setAngleDetection(false)
      .optimize(mesh);
    splitSelfPairedCut(mesh);
    const size_t requiredTriangles = protectFixedGeometry(mesh);
    const size_t cellsAfter = mesh.getCellCount();
    return {std::move(mesh),
      {hmin, hmax, hausdorff, requiredTriangles, cellsBefore, cellsAfter}};
  }

  SphereDiscretization Sphere::prepareWNGIRBackground() const
  {
    const Real h = m_configuration.getH();
    const Real hmin = m_configuration.backgroundHMin * h;
    const Real hmax = m_configuration.backgroundHMax * h;
    const Real hausdorff = m_configuration.backgroundHausdorff * h;
    MMG::Mesh mesh(makeUniformChamber());
    const size_t cellsBefore = mesh.getCellCount();
    protectFixedGeometry(mesh);
    MMG::Optimizer()
      .setHMin(hmin)
      .setHMax(hmax)
      .setHausdorff(hausdorff)
      .setGradation(m_configuration.backgroundGradation)
      .setAngleDetection(false)
      .optimize(mesh);
    splitSelfPairedCut(mesh);
    const size_t requiredTriangles = protectFixedGeometry(mesh);

    const size_t cellsAfter = mesh.getCellCount();
    return {std::move(mesh),
      {hmin, hmax, hausdorff, requiredTriangles, cellsBefore, cellsAfter}};
  }
}
