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
#include <Rodin/Adaptation/WNGIR/Loss.h>
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
    cube.scale(m_configuration.getGridSpacing());

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

  void Sphere::projectFixedGeometry(MMG::Mesh& mesh) const
  {
    constexpr unsigned outer = 1, xy = 2, plus = 4, minus = 8;
    // Relative roundoff floor, not an additional mesh-quality budget.
    constexpr Real roundoff = 64 * std::numeric_limits<Real>::epsilon();
    std::vector<unsigned> constraints(mesh.getVertexCount(), 0);
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (!attribute)
        continue;
      const unsigned constraint = *attribute == Outer               ? outer
        : (*attribute == SigmaXYPlus || *attribute == SigmaXYMinus) ? xy
        : *attribute == SigmaPlus                                   ? plus
        : *attribute == SigmaMinus                                  ? minus
                                                                    : 0;
      for (const Index vertex : face->getVertices())
        constraints[vertex] |= constraint;
    }

    std::vector<Math::SpatialPoint> coordinates;
    coordinates.reserve(mesh.getVertexCount());
    Real maximumCorrection = 0;
    size_t correctedVertices = 0;
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const auto original = mesh.getVertexCoordinates(vertex);
      Math::SpatialPoint projected = original;
      std::array<Math::SpatialPoint, 3> basis;
      std::array<Real, 3> offsets{};
      size_t rank = 0;
      // Orthogonalize the affine equations together. Successive projections
      // onto the original planes would not preserve their intersections.
      for (const unsigned constraint : {outer, xy, plus, minus})
      {
        if (!(constraints[vertex] & constraint))
          continue;
        Math::SpatialPoint normal(3);
        normal(0) = (constraint == outer || constraint == xy) ? 1 : 0;
        normal(1) = constraint == outer ? 0 : constraint == xy ? -1 : 1;
        normal(2) = constraint == plus ? -1 : constraint == minus ? 1 : 0;
        Real offset = constraint == outer ? m_configuration.outerRadius : 0;
        for (size_t i = 0; i < rank; ++i)
        {
          const Real coefficient = normal.dot(basis[i]);
          normal -= coefficient * basis[i];
          offset -= coefficient * offsets[i];
        }
        const Real norm = normal.norm();
        if (norm <= roundoff)
        {
          if (std::abs(offset) > roundoff * m_configuration.outerRadius)
            throw std::runtime_error(
              "Incompatible fixed-boundary labels at a chamber vertex.");
          continue;
        }
        basis[rank] = normal / norm;
        offsets[rank] = offset / norm;
        ++rank;
      }
      for (size_t i = 0; i < rank; ++i)
        projected += (offsets[i] - basis[i].dot(original)) * basis[i];
      const Real correction = (projected - original).norm();
      maximumCorrection = std::max(maximumCorrection, correction);
      correctedVertices += correction > roundoff * m_configuration.outerRadius;
      coordinates.push_back(std::move(projected));
    }

    // Validate the proposed geometry before mutating the mesh, including
    // cells incident to vertices shared by cuts and the material interface.
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& vertices = cell->getVertices();
      Math::SpatialMatrix<Real> before(3, 3), after(3, 3);
      for (size_t i = 0; i < 3; ++i)
      {
        const auto oldEdge = mesh.getVertexCoordinates(vertices(i + 1)) -
          mesh.getVertexCoordinates(vertices(0));
        const auto newEdge = coordinates[vertices(i + 1)] - coordinates[vertices(0)];
        for (size_t j = 0; j < 3; ++j)
        {
          before(j, i) = oldEdge(j);
          after(j, i) = newEdge(j);
        }
      }
      const Real determinant = after.determinant();
      const Real scale = after.col(0).norm() * after.col(1).norm() * after.col(2).norm();
      if (!std::isfinite(determinant) || determinant * before.determinant() <= 0 ||
        std::abs(determinant) <= roundoff * scale)
        throw std::runtime_error(
          "Fixed-boundary projection would invert or degenerate a tetrahedron.");
    }
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      if (constraints[vertex])
        mesh.setVertexCoordinates(vertex, coordinates[vertex]);
    Alert::Text<Alert::YellowT> heading(Alert::Yellow, "Fixed-boundary projection");
    Alert::Info() << heading.setBold() << Alert::NewLine
                  << "Corrected vertices:                 "
                  << Alert::Notation::Number(correctedVertices) << Alert::NewLine
                  << "Maximum displacement:               "
                  << Alert::Notation::Number(maximumCorrection) << Alert::Raise;
  }

  SphereDiscretization Sphere::discretize(bool conformingCuts,
    Real requestedWelschScale) const
  {
    MMG::Mesh mesh(makeUniformChamber());
    const Real h = m_configuration.getGridSpacing();
    const size_t cellsBefore = mesh.getCellCount();
    protectFixedGeometry(mesh, conformingCuts);
    const Real hmin = m_configuration.hmin;
    const Real hmax = m_configuration.hmax;
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
      .setBoundaryReference(Gamma)
      .setAngleDetection(false);
    mesh = discretizer.discretize(sphere);
    projectFixedGeometry(mesh);
    splitSelfPairedCut(mesh);
    if (m_configuration.adapt && !conformingCuts)
    {
      adapt(mesh, requestedWelschScale);
    }
    else
    {
      protectFixedGeometry(mesh, conformingCuts);
      MMG::Optimizer()
        .setHMin(hmin)
        .setHMax(hmax)
        .setHausdorff(hausdorff)
        .setGradation(2)
        .setAngleDetection(false)
        .optimize(mesh);
      projectFixedGeometry(mesh);
      splitSelfPairedCut(mesh);
    }
    const size_t requiredTriangles = protectFixedGeometry(mesh, conformingCuts);
    const size_t cellsAfter = mesh.getCellCount();
    return {std::move(mesh),
      {hmin, hmax, hausdorff, requiredTriangles, cellsBefore, cellsAfter}};
  }

  SphereDiscretization Sphere::prepareWNGIRBackground(Real requestedWelschScale) const
  {
    MMG::Mesh mesh(makeUniformChamber());
    const Real h = m_configuration.getGridSpacing();
    const Real interfaceSize = m_configuration.hmin;
    const Real farSize = m_configuration.hmax;
    const Real hmin = m_configuration.hmin;
    const Real hmax = m_configuration.hmax;
    const Real hausdorff = m_configuration.backgroundHausdorff * h;
    const Real backgroundMeanElementSize = meanElementSize(mesh);
    const Real welschScale = requestedWelschScale > 0
      ? requestedWelschScale : Real(3) * h;
    const size_t cellsBefore = mesh.getCellCount();
    if (m_configuration.adapt)
    {
      P1<Real, Mesh> sizeSpace(mesh);
      MMG::RealGridFunction size(sizeSpace);
      const Adaptation::WNGIRLoss welsch(welschScale);
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        const Real distance =
          std::abs(mesh.getVertexCoordinates(vertex).norm() - Real(1));
        size[vertex] =
          farSize - (farSize - interfaceSize) * welsch.getWeight(distance);
      }
      protectFixedGeometry(mesh, false);
      MMG::Adapt()
        .setHMin(hmin)
        .setHMax(hmax)
        .setHausdorff(hausdorff)
        .setGradation(m_configuration.adaptGradation)
        .setAngleDetection(false)
        .adapt(mesh, size);
    }
    else
    {
      protectFixedGeometry(mesh, true);
      MMG::Optimizer()
        .setHMin(hmin)
        .setHMax(hmax)
        .setHausdorff(hausdorff)
        .setGradation(m_configuration.backgroundGradation)
        .setAngleDetection(false)
        .optimize(mesh);
    }
    projectFixedGeometry(mesh);
    splitSelfPairedCut(mesh);
    const size_t requiredTriangles = protectFixedGeometry(mesh, !m_configuration.adapt);

    const size_t cellsAfter = mesh.getCellCount();
    ReconstructionDiagnostics diagnostics{
      hmin, hmax, hausdorff, requiredTriangles, cellsBefore, cellsAfter};
    diagnostics.backgroundMeanElementSize = backgroundMeanElementSize;
    diagnostics.welschScale = welschScale;
    return {std::move(mesh), diagnostics};
  }

  void Sphere::adapt(MMG::Mesh& fitted, Real requestedWelschScale) const
  {
    // MMG and the subsequent exact-boundary projection operate on a candidate.
    // Any exception leaves the fitted geometry, topology and MMG tags intact.
    MMG::Mesh mesh = fitted;
    const Real h = m_configuration.getGridSpacing();
    const Real interfaceSize = m_configuration.hmin;
    const Real farSize = m_configuration.hmax;
    const Real welschScale = requestedWelschScale > 0
      ? requestedWelschScale : Real(3) * h;
    const Adaptation::WNGIRLoss welsch(welschScale);

    P1<Real, Mesh> sizeSpace(mesh);
    MMG::RealGridFunction size(sizeSpace);
    mesh.getConnectivity().compute(0, 0);
    mesh.getConnectivity().compute(0, mesh.getDimension());
    GridFunction distance(sizeSpace);
    Distance::Eikonal(distance).setInterior(Obstacle).setInterface(Gamma).solve();
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      size[vertex] = farSize - (farSize - interfaceSize)
        * welsch.getWeight(std::abs(distance[vertex]));
    }

    protectFixedGeometry(mesh, false);
    MMG::Adapt()
      .setHMin(interfaceSize)
      .setHMax(farSize)
      .setHausdorff(0.1 * h * h)
      .setGradation(m_configuration.adaptGradation)
      .setAngleDetection(false)
      .adapt(mesh, size);
    projectFixedGeometry(mesh);
    splitSelfPairedCut(mesh);
    fitted = std::move(mesh);
  }
}
