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

#include <Rodin/Distance/Eikonal.h>
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

  size_t Sphere::protectFixedGeometry(MMG::Mesh& mesh, bool cuts) const
  {
    mesh.getRequiredTriangles().clear();
    mesh.getRidges().clear();
    mesh.getCorners().clear();
    const FlatSet<Attribute> fixed{
      Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
    const size_t faceDimension = mesh.getDimension() - 1;
    size_t count = 0;
    for (auto face = mesh.getPolytope(faceDimension); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (attribute && (*attribute == Outer || (cuts && fixed.contains(*attribute))))
      {
        mesh.setRequiredTriangle(face->getIndex());
        ++count;
      }
    }
    if (cuts)
      return count;

    // Edge labels returned by a previous MMG call describe the previous
    // interface; only the features derived below are handed on.
    mesh.getConnectivity().compute(faceDimension, 1);
    for (auto edge = mesh.getPolytope(1); edge; ++edge)
    {
      if (edge->getAttribute())
        mesh.setAttribute({1, edge->getIndex()}, {});
    }
    const auto& faceEdges = mesh.getConnectivity().getIncidence(faceDimension, 1);
    std::map<Index, FlatSet<Attribute>> edgeLabels;
    std::map<Index, FlatSet<Attribute>> vertexLabels;
    for (auto face = mesh.getPolytope(faceDimension); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (!attribute)
        continue;
      for (const Index edge : faceEdges.at(face->getIndex()))
        edgeLabels[edge].insert(*attribute);
      for (const Index vertex : face->getVertices())
        vertexLabels[vertex].insert(*attribute);
    }
    const auto onFixedFace = [&](const FlatSet<Attribute>& labels) {
      return std::any_of(labels.begin(), labels.end(),
        [&](Attribute label) { return fixed.contains(label); });
    };
    // The two halves of the x = y cut lie in one plane. The line between them
    // is a reference edge, which MMG keeps as a line of edges; a ridge there
    // would join two identical normals and leave its tangent undefined.
    const auto plane = [](Attribute label) {
      return label == SigmaXYMinus ? SigmaXYPlus : label;
    };
    for (const auto& [edge, labels] : edgeLabels)
    {
      if (labels.size() < 2 || !onFixedFace(labels))
        continue;
      FlatSet<Attribute> planes;
      for (const Attribute label : labels)
        planes.insert(plane(label));
      mesh.setAttribute({1, edge}, Ridge);
      if (planes.size() >= 2)
        mesh.setRidge(edge);
    }
    for (const auto& [vertex, labels] : vertexLabels)
    {
      if (labels.size() >= 3 && onFixedFace(labels))
        mesh.setCorner(vertex);
    }
    return count;
  }

  SphereDiscretization Sphere::discretize() const
  {
    const Real h = m_configuration.getH();
    MMG::Mesh mesh(makeUniformChamber());
    const size_t cellsBefore = mesh.getCellCount();
    protectFixedGeometry(mesh, false);
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
    if (m_configuration.adapt)
    {
      adapt(mesh, h);
    }
    else
    {
      protectFixedGeometry(mesh, false);
      MMG::Optimizer()
        .setHMin(hmin)
        .setHMax(hmax)
        .setHausdorff(hausdorff)
        .setGradation(2)
        .setAngleDetection(false)
        .optimize(mesh);
      splitSelfPairedCut(mesh);
    }
    const size_t requiredTriangles = protectFixedGeometry(mesh, false);
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
    protectFixedGeometry(mesh, true);
    MMG::Optimizer()
      .setHMin(hmin)
      .setHMax(hmax)
      .setHausdorff(hausdorff)
      .setGradation(m_configuration.backgroundGradation)
      .setAngleDetection(false)
      .optimize(mesh);
    splitSelfPairedCut(mesh);
    const size_t requiredTriangles = protectFixedGeometry(mesh, true);

    const size_t cellsAfter = mesh.getCellCount();
    return {std::move(mesh),
      {hmin, hmax, hausdorff, requiredTriangles, cellsBefore, cellsAfter}};
  }

  void Sphere::adapt(MMG::Mesh& mesh, Real h) const
  {
    const Real interfaceSize = m_configuration.adaptInterfaceSize * h;
    const Real farSize = m_configuration.adaptFarSize * h;
    const Real width = m_configuration.adaptWidth * h;

    P1<Real, Mesh> sizeSpace(mesh);
    MMG::RealGridFunction size(sizeSpace);
    mesh.getConnectivity().compute(0, 0);
    mesh.getConnectivity().compute(0, mesh.getDimension());
    GridFunction distance(sizeSpace);
    Distance::Eikonal(distance).setInterior(Obstacle).setInterface(Gamma).solve();
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const Real ratio = std::min(Real(1), std::abs(distance[vertex]) / width);
      size[vertex] = interfaceSize + (farSize - interfaceSize) * ratio;
    }

    protectFixedGeometry(mesh, false);
    MMG::Adapt()
      .setHMin(0.1 * h)
      .setHMax(std::max(interfaceSize, farSize))
      .setHausdorff(0.1 * h * h)
      .setGradation(m_configuration.adaptGradation)
      .setAngleDetection(false)
      .adapt(mesh, size);
    splitSelfPairedCut(mesh);
  }
}
