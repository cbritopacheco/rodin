/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "ProblemData.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <Eigen/Eigenvalues>
#include <Rodin/Alert.h>
#include <Rodin/Solver/SparseLU.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples::Heart::AtriumContraction
{
  Options readOptions(int argc, char** argv)
  {
    Options out;
    for (int i = 1; i < argc; ++i)
    {
      const std::string key(argv[i]);
      if (key == "--help")
      {
        std::cout << "AtriumContraction [--mesh path] [--output stem] [--steps 24]\n"
                     "  [--contraction .03] [--speed 200] [--width 6]\n"
                     "  [--young 50] [--poisson .3] [--active-stiffness 20]\n"
                     "  [--density 1060 (kg/m^3)]\n"
                     "  [--stop-after n] [--base-label 6] [--apex x y z]\n"
                     "  [--apex-radius 3] [--fibers-only]\n";
        std::exit(0);
      }
      if (key == "--fibers-only")
      {
        out.fibersOnly = true;
        continue;
      }
      if (key == "--apex")
      {
        if (i + 3 >= argc)
        {
          Alert::Exception() << "--apex requires three coordinates in mm" << Alert::Raise;
        }
        for (size_t j = 0; j < 3; ++j)
          out.appendageApex(j) = std::stod(argv[++i]);
        continue;
      }
      if (++i == argc)
      {
        Alert::Exception() << "Missing value for " + key << Alert::Raise;
      }
      const std::string value(argv[i]);
      if (key == "--mesh")
        out.mesh = value;
      else if (key == "--output")
        out.output = value;
      else if (key == "--steps")
        out.steps = std::stoul(value);
      else if (key == "--stop-after")
        out.stopAfter = std::stoul(value);
      else if (key == "--base-label")
        out.baseLabel = std::stoul(value);
      else if (key == "--apex-radius")
        out.apexRadius = std::stod(value);
      else if (key == "--contraction")
        out.contraction = std::stod(value);
      else if (key == "--speed")
        out.speed = std::stod(value);
      else if (key == "--width")
        out.width = std::stod(value);
      else if (key == "--young")
        out.young = std::stod(value);
      else if (key == "--poisson")
        out.poisson = std::stod(value);
      else if (key == "--density")
        out.density = std::stod(value);
      else if (key == "--active-stiffness")
        out.activeStiffness = std::stod(value);
      else
      {
        Alert::Exception() << "Unknown option " + key << Alert::Raise;
      }
    }
    if (out.steps == 0 || !std::isfinite(out.contraction) || out.contraction <= 0 ||
      out.contraction >= 0.5 || !std::isfinite(out.speed) || out.speed <= 0 ||
      !std::isfinite(out.width) || out.width <= 0 || !std::isfinite(out.young) ||
      out.young <= 0 || !std::isfinite(out.poisson) || out.poisson <= -1 ||
      out.poisson >= 0.5 || !std::isfinite(out.activeStiffness) ||
      out.activeStiffness <= 0 || !out.appendageApex.allFinite() ||
      !std::isfinite(out.apexRadius) || out.apexRadius <= 0 ||
      !std::isfinite(out.density) || out.density <= 0)
    {
      Alert::Exception() << "Invalid material, wave, or time discretization parameters"
                         << Alert::Raise;
    }
    return out;
  }

  Vec point3(const LocalMesh& mesh, Index vertex)
  {
    const auto x = mesh.getVertexCoordinates(vertex);
    return Vec(x(0), x(1), x(2));
  }

  // The wedge determinant is affine in triangle coordinates and quadratic
  // through the thickness. Its global minimum occurs at a triangle vertex,
  // at a thickness endpoint or at an interior minimum of that quadratic.
  Real minimumJacobian(const std::array<Vec, 6>& q)
  {
    Real minimum = std::numeric_limits<Real>::infinity();
    for (size_t vertex = 0; vertex < 3; ++vertex)
    {
      const Vec dt = q[vertex + 3] - q[vertex];
      const auto determinant = [&](Real t) {
        const Vec dr = (1 - t) * (q[1] - q[0]) + t * (q[4] - q[3]);
        const Vec ds = (1 - t) * (q[2] - q[0]) + t * (q[5] - q[3]);
        return dr.cross(ds).dot(dt);
      };
      const Real c = determinant(0), end = determinant(1), mid = determinant(0.5);
      const Real a = 2 * (end + c - 2 * mid), b = end - c - a;
      minimum = std::min({minimum, c, end});
      if (a > 0)
      {
        const Real t = -b / (2 * a);
        if (t > 0 && t < 1)
          minimum = std::min(minimum, c + t * (b + a * t));
      }
    }
    return minimum;
  }
  LocalMesh loadWall(const Options& cfg)
  {
    LocalMesh mesh;
    mesh.load(cfg.mesh, IO::FileFormat::MEDIT);
    if (mesh.getDimension() != 3 || mesh.getSpaceDimension() != 3)
    {
      Alert::Exception() << "Expected a 3D wall volume mesh" << Alert::Raise;
    }
    mesh.getConnectivity().compute(2, 3);
    return mesh;
  }

  ProblemData::ProblemData(LocalMesh& wall, const Options& options)
    : mesh(wall),
      cfg(options),
      scalarSpace(mesh),
      fiberSpace(mesh, 3),
      diagnosticSpace(mesh),
      fiberCoordinate(scalarSpace),
      waveCoordinate(scalarSpace),
      fiberField(fiberSpace),
      normalField(fiberSpace),
      longitudinalField(fiberSpace),
      fallbackField(diagnosticSpace),
      directions(mesh.getCellCount()),
      normals(mesh.getCellCount()),
      longitudinal(mesh.getCellCount()),
      fallback(mesh.getCellCount(), 0)
  {
    inspectGeometry();
    generateFibers();
    minimumProjection = waveCoordinate.getData().minCoeff();
    maximumProjection = waveCoordinate.getData().maxCoeff();
    constexpr Real PulsePadding = 3;
    endTime = (maximumProjection + PulsePadding * cfg.width) / cfg.speed;
  }

  void ProblemData::inspectGeometry()
  {
    mesh.getConnectivity().compute(2, 3);
    std::vector<bool> innerVertex(mesh.getVertexCount(), false);
    std::set<Index> baseVertices;
    size_t badCells = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (cell->getGeometry() != Polytope::Type::Wedge ||
        cell->getAttribute() != Attribute{4})
      {
        Alert::Exception() << "The mechanics input must contain only material-4 wedges"
                           << Alert::Raise;
      }
      std::array<Vec, 6> q;
      for (size_t j = 0; j < 6; ++j)
      {
        const auto x = mesh.getVertexCoordinates(cell->getVertices()(j));
        for (size_t k = 0; k < 3; ++k)
          q[j](k) = x(k);
      }
      const Real determinant = minimumJacobian(q);
      const Real area = (q[1] - q[0]).cross(q[2] - q[0]).norm();
      minimumReferenceJacobian = std::min(minimumReferenceJacobian, determinant / area);
      if (!(determinant > 0) || !std::isfinite(determinant))
        ++badCells;
      normals[cell->getIndex()] =
        ((q[3] - q[0]) + (q[4] - q[1]) + (q[5] - q[2])).normalized();
      for (size_t j = 0; j < 3; ++j)
        innerVertex[cell->getVertices()(j)] = true;
    }
    if (badCells != 0)
    {
      Alert::Exception() << std::to_string(badCells) +
          " folded/degenerate wedges; repair before mechanics"
                         << Alert::Raise;
    }
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      const auto tag = face->getAttribute();
      if (tag == cfg.baseLabel)
      {
        for (Index vertex : face->getVertices())
        {
          baseAllVertices.insert(vertex);
          if (innerVertex[vertex])
            baseVertices.insert(vertex);
        }
      }

      if (tag == Attribute{1} || tag == Attribute{2} || tag == Attribute{3} ||
        tag == Attribute{5} || tag == Attribute{6})
        ++rimFaces;
    }
    if (rimFaces == 0)
    {
      Alert::Exception() << "No tagged inlet/outlet rim support faces found"
                         << Alert::Raise;
    }
    if (baseVertices.size() < 3)
    {
      Alert::Exception() << "Base label " << cfg.baseLabel
                         << " needs at least three inner-rim vertices" << Alert::Raise;
    }
    Vec wallCenter = Vec::Zero();
    for (Index v : baseVertices)
      baseOrigin += point3(mesh, v);
    baseOrigin /= baseVertices.size();
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (Index v : baseVertices)
    {
      const Vec dx = point3(mesh, v) - baseOrigin;
      covariance += dx * dx.transpose();
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> plane(covariance);
    if (plane.info() != Eigen::Success || plane.eigenvalues()(1) <= 0)
    {
      Alert::Exception() << "The basal rim does not define a plane" << Alert::Raise;
    }
    waveDirection = plane.eigenvectors().col(0);
    for (auto v = mesh.getVertex(); v; ++v)
      wallCenter += point3(mesh, v->getIndex());
    wallCenter /= mesh.getVertexCount();
    if (waveDirection.dot(wallCenter - baseOrigin) < 0)
      waveDirection *= -1;
  }

  void ProblemData::generateFibers()
  {
    // One harmonic base-to-appendage coordinate supplies a simple anatomical
    // direction field. This omits the paper's regional bundle classifications.
    constexpr Attribute ApexPatchLabel = 8;
    std::set<std::array<Index, 3>> apexKeys;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      Vec center = Vec::Zero();
      for (size_t j = 0; j < 3; ++j)
        center += point3(mesh, cell->getVertices()(j));
      center /= 3;
      if ((center - cfg.appendageApex).norm() > cfg.apexRadius)
        continue;
      for (size_t layer = 0; layer < 2; ++layer)
      {
        std::array<Index, 3> key;
        for (size_t j = 0; j < 3; ++j)
          key[j] = cell->getVertices()(3 * layer + j);
        std::sort(key.begin(), key.end());
        apexKeys.insert(key);
      }
    }
    std::vector<std::pair<Index, Optional<Attribute>>> restoredLabels;
    std::set<Index> apexVertices;
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      if (face->getAttribute() == ApexPatchLabel)
      {
        Alert::Exception() << "Boundary label 8 is reserved for the temporary apex patch"
                           << Alert::Raise;
      }
      if (face->getGeometry() != Polytope::Type::Triangle)
        continue;
      std::array<Index, 3> key;
      for (size_t j = 0; j < 3; ++j)
        key[j] = face->getVertices()(j);
      std::sort(key.begin(), key.end());
      if (!apexKeys.contains(key))
        continue;
      restoredLabels.emplace_back(face->getIndex(), face->getAttribute());
      apexVertices.insert(key.begin(), key.end());
    }
    if (restoredLabels.empty())
    {
      Alert::Exception() << "No wall triangles lie in the appendage-apex patch; check "
                            "--apex and --apex-radius"
                         << Alert::Raise;
    }
    for (Index vertex : apexVertices)
    {
      if (baseAllVertices.contains(vertex))
      {
        Alert::Exception() << "The appendage patch overlaps the base rim" << Alert::Raise;
      }
    }
    for (const auto& item : restoredLabels)
      mesh.setAttribute({2, item.first}, ApexPatchLabel);
    TrialFunction coordinate(scalarSpace);
    TestFunction coordinateTest(scalarSpace);
    Problem laplace(coordinate, coordinateTest);
    laplace = Integral(Grad(coordinate), Grad(coordinateTest)) +
      DirichletBC(coordinate, RealFunction(0)).on(cfg.baseLabel) +
      DirichletBC(coordinate, RealFunction(1)).on(ApexPatchLabel);
    Solver::SparseLU coordinateSolver(laplace);
    coordinateSolver.solve();
    const auto& coordinateSystem = laplace.getLinearSystem();
    coordinateResidual =
      (coordinateSystem.getOperator() * coordinateSystem.getSolution() -
        coordinateSystem.getVector())
        .norm() /
      std::max<Real>(1, coordinateSystem.getVector().norm());
    if (!std::isfinite(coordinateResidual) || coordinateResidual > 1e-8)
    {
      Alert::Exception() << "The fiber-coordinate Laplace solve failed: "
                         << coordinateResidual << Alert::Raise;
    }
    for (const auto& item : restoredLabels)
      mesh.setAttribute({2, item.first}, item.second);
    fiberCoordinate = coordinate.getSolution();
    fiberCoordinate.setName("FiberCoordinate");
    apexPatchFaces = restoredLabels.size();
    const auto coordinateGradient = Grad(fiberCoordinate);
    constexpr Real DirectionTolerance = 1e-10;
    const Real gradientTolerance = DirectionTolerance /
      std::max(cfg.apexRadius, (cfg.appendageApex - baseOrigin).norm());
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const auto gradient = coordinateGradient.getValue(point);
      Vec direction(gradient(0), gradient(1), gradient(2));
      const Vec normal = normals[cell->getIndex()];
      if (direction.norm() > gradientTolerance)
        direction.normalize();
      else
        direction.setZero();
      direction -= direction.dot(normal) * normal;
      if (direction.norm() < DirectionTolerance)
      {
        fallback[cell->getIndex()] = 1;
        ++fallbackCells;
        const auto x = point.getPhysicalCoordinates();
        direction = cfg.appendageApex - Vec(x(0), x(1), x(2));
        if (direction.norm() > 0)
          direction.normalize();
        direction -= direction.dot(normal) * normal;
      }
      if (direction.norm() < DirectionTolerance)
      {
        fallback[cell->getIndex()] = 2;
        Eigen::Index axis;
        normal.cwiseAbs().minCoeff(&axis);
        direction = Vec::Unit(axis) - normal(axis) * normal;
      }
      longitudinal[cell->getIndex()] = direction.normalized();
      directions[cell->getIndex()] = direction.normalized().cross(normal).normalized();
    }
    const auto vectorField = [&](const std::vector<Vec>& values) {
      return VectorFunction(3, [&values](const Point& point) {
        Math::SpatialVector<Real> result(3);
        result = values.at(point.getPolytope().getIndex());
        return result;
      });
    };
    fiberField.setName("FiberDirection");
    fiberField = vectorField(directions);
    normalField.setName("WallNormal");
    normalField = vectorField(normals);
    longitudinalField.setName("LongitudinalDirection");
    longitudinalField = vectorField(longitudinal);
    fallbackField.setName("FiberFallback");
    fallbackField = [&](const Point& point) {
      return fallback.at(point.getPolytope().getIndex());
    };
    waveCoordinate.setName("WaveCoordinate");
    waveCoordinate = [&](const Point& point) { return waveCoordinateAt(point); };
  }

  Real ProblemData::waveCoordinateAt(const Point& point) const
  {
    const auto x = point.getPhysicalCoordinates();
    return waveDirection.dot(Vec(x(0), x(1), x(2)) - baseOrigin);
  }

  Real ProblemData::activationAt(const Point& point, Real time) const
  {
    const Real argument = (waveCoordinateAt(point) - cfg.speed * time) / cfg.width;
    return std::exp(-argument * argument);
  }

  Math::SpatialVector<Real> ProblemData::fiberAt(const Point& point) const
  {
    Math::SpatialVector<Real> fiber(3);
    fiber = directions.at(point.getPolytope().getIndex());
    return fiber;
  }
}
