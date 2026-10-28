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
#include <boost/mpi/collectives.hpp>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples::Heart::Atrium
{
  Options readOptions(int argc, char** argv, std::string& meshPath)
  {
    Options out;
    for (int i = 1; i < argc; ++i)
    {
      const std::string key(argv[i]);
      // PETSc consumes single-dash solver options; double-dash model options
      // remain owned by this parser.
      if (key.size() > 1 && key[0] == '-' && key[1] != '-')
      {
        while (i + 1 < argc && argv[i + 1][0] != '-')
          ++i;
        continue;
      }
      if (key == "--help")
      {
        std::cout << "AtriumContraction [--mesh path] [--output stem] [--steps 24]\n"
                     "  [--activation 20] [--relaxation 10] (1/s)\n"
                     "  [--matrix-scale 2.92] [--matrix-exponent 5.6]\n"
                     "  [--fiber-scale 11.84] [--fiber-exponent 17.95] [--volume-penalty 500]\n"
                     "  [--viscosity .2] [--active-damping 3] (kPa s)\n"
                     "  [--active-stiffness 20] [--cross-bridge-stiffness 100]\n"
                     "  [--contractility 20] [--destruction-rate .4]\n"
                     "  [--density 1060] (kg/m^3) [--speed 200] [--width 6]\n"
                     "  [--stop-after n] [--dt seconds] [--base-label 6] [--apex x y z]\n"
                     "  [--apex-radius 3] [--fibers-only] [--cache-constitutive]\n";
        std::exit(0);
      }
      if (key == "--cache-constitutive")
      {
        out.cacheConstitutive = true;
        continue;
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
        meshPath = value;
      else if (key == "--output")
        out.output = value;
      else if (key == "--steps")
        out.steps = std::stoul(value);
      else if (key == "--stop-after")
        out.stopAfter = std::stoul(value);
      else if (key == "--dt")
        out.timeStep = std::stod(value);
      else if (key == "--base-label")
        out.baseLabel = std::stoul(value);
      else if (key == "--apex-radius")
        out.apexRadius = std::stod(value);
      else if (key == "--activation")
        out.activation = std::stod(value);
      else if (key == "--relaxation")
        out.relaxation = std::stod(value);
      else if (key == "--speed")
        out.speed = std::stod(value);
      else if (key == "--width")
        out.width = std::stod(value);
      else if (key == "--matrix-scale")
        out.passive.matrixScale = std::stod(value);
      else if (key == "--matrix-exponent")
        out.passive.matrixExponent = std::stod(value);
      else if (key == "--fiber-scale")
        out.passive.fiberScale = std::stod(value);
      else if (key == "--fiber-exponent")
        out.passive.fiberExponent = std::stod(value);
      else if (key == "--volume-penalty")
        out.passive.volumePenalty = std::stod(value);
      else if (key == "--viscosity")
        out.viscosity = std::stod(value);
      else if (key == "--active-damping")
        out.activeDamping = std::stod(value);
      else if (key == "--cross-bridge-stiffness")
        out.crossBridgeStiffness = std::stod(value);
      else if (key == "--contractility")
        out.contractility = std::stod(value);
      else if (key == "--destruction-rate")
        out.destructionRate = std::stod(value);
      else if (key == "--density")
        out.density = std::stod(value);
      else if (key == "--active-stiffness")
        out.activeStiffness = std::stod(value);
      else
      {
        Alert::Exception() << "Unknown option " + key << Alert::Raise;
      }
    }
    return validateOptions(out);
  }

  Options validateOptions(const Options& cfg)
  {
    (void)Rodin::Solid::DemirayHolzapfelGasserOgden(cfg.passive);
    if (cfg.steps == 0 || !std::isfinite(cfg.timeStep) || cfg.timeStep < 0 ||
      !std::isfinite(cfg.activation) || cfg.activation < 0 ||
      !std::isfinite(cfg.relaxation) || cfg.relaxation < 0 ||
      !std::isfinite(cfg.viscosity) || cfg.viscosity < 0 ||
      !std::isfinite(cfg.activeDamping) || cfg.activeDamping < 0 ||
      !std::isfinite(cfg.crossBridgeStiffness) || cfg.crossBridgeStiffness <= 0 ||
      !std::isfinite(cfg.contractility) || cfg.contractility < 0 ||
      !std::isfinite(cfg.destructionRate) || cfg.destructionRate < 0 ||
      !std::isfinite(cfg.speed) || cfg.speed <= 0 || !std::isfinite(cfg.width) ||
      cfg.width <= 0 || !std::isfinite(cfg.activeStiffness) || cfg.activeStiffness <= 0 ||
      !cfg.appendageApex.allFinite() || !std::isfinite(cfg.apexRadius) ||
      cfg.apexRadius <= 0 || !std::isfinite(cfg.density) || cfg.density <= 0)
    {
      Alert::Exception() << "Invalid material, wave, or time discretization parameters"
                         << Alert::Raise;
    }
    return cfg;
  }

  Vec point3(const MeshType& mesh, Index vertex)
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

  ProblemData::ProblemData(MeshType& wall, const Options& options)
    : mesh(wall),
      cfg(options),
      scalarSpace(mesh),
      fiberSpace(mesh, 3),
      diagnosticSpace(mesh),
      fiberCoordinate(scalarSpace),
      transmuralCoordinate(scalarSpace),
      waveCoordinate(scalarSpace),
      fiberField(fiberSpace),
      normalField(fiberSpace),
      longitudinalField(fiberSpace),
      fallbackField(diagnosticSpace),
      directions(mesh.getShard().getPolytopeCount(3)),
      normals(mesh.getShard().getPolytopeCount(3)),
      longitudinal(mesh.getShard().getPolytopeCount(3)),
      fallback(mesh.getShard().getPolytopeCount(3), 0),
      m_coordinate(scalarSpace),
      m_coordinateTest(scalarSpace),
      m_laplace(m_coordinate, m_coordinateTest),
      m_coordinateSolver(m_laplace)
  {
    inspectGeometry();
    generateTransmuralCoordinate();
    generateFibers();
    Index selected;
    minimumProjection = waveCoordinate.min(selected);
    maximumProjection = waveCoordinate.max(selected);
    constexpr Real PulsePadding = 3;
    endTime = (maximumProjection + PulsePadding * cfg.width) / cfg.speed;
  }

  void ProblemData::inspectGeometry()
  {
    const auto& comm = mesh.getContext().getCommunicator();
    const auto& shard = mesh.getShard();
    mesh.getConnectivity().compute(2, 3);
    std::vector<bool> innerVertex(shard.getVertexCount(), false);
    std::set<Index> baseVertices;
    size_t badCells = 0, invalidCells = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (cell->getGeometry() != Polytope::Type::Wedge)
      {
        invalidCells += shard.isOwned(3, cell->getIndex());
        continue;
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
        badCells += shard.isOwned(3, cell->getIndex());
      for (size_t j = 0; j < 3; ++j)
        innerVertex[cell->getVertices()(j)] = true;
    }
    invalidCells = boost::mpi::all_reduce(comm, invalidCells, std::plus<size_t>());
    if (invalidCells)
      Alert::Exception() << "The solid submesh must contain only wedges" << Alert::Raise;
    badCells = boost::mpi::all_reduce(comm, badCells, std::plus<size_t>());
    minimumReferenceJacobian =
      boost::mpi::all_reduce(comm, minimumReferenceJacobian, boost::mpi::minimum<Real>());
    if (badCells != 0)
    {
      Alert::Exception() << std::to_string(badCells) +
          " folded/degenerate wedges; repair before mechanics"
                         << Alert::Raise;
    }
    for (auto face = mesh.getFace(); face; ++face)
    {
      if (!shard.isBoundary(face->getIndex()))
        continue;
      const auto tag = face->getAttribute();
      if (tag == cfg.baseLabel)
      {
        for (Index vertex : face->getVertices())
        {
          baseAllVertices.insert(vertex);
          if (innerVertex[vertex] && shard.isOwned(0, vertex))
            baseVertices.insert(vertex);
        }
      }

      if (tag == Attribute{1} || tag == Attribute{2} || tag == Attribute{3} ||
        tag == Attribute{5} || tag == Attribute{6})
      {
        if (shard.isOwned(2, face->getIndex()))
          ++rimFaces;
      }
    }
    rimFaces = boost::mpi::all_reduce(comm, rimFaces, std::plus<size_t>());
    if (rimFaces == 0)
    {
      Alert::Exception() << "No tagged inlet/outlet rim support faces found"
                         << Alert::Raise;
    }
    const size_t baseCount =
      boost::mpi::all_reduce(comm, baseVertices.size(), std::plus<size_t>());
    if (baseCount < 3)
    {
      Alert::Exception() << "Base label " << cfg.baseLabel
                         << " needs at least three inner-rim vertices" << Alert::Raise;
    }
    Vec wallCenter = Vec::Zero();
    for (Index v : baseVertices)
      baseOrigin += point3(mesh, v);
    for (size_t k = 0; k < 3; ++k)
    {
      baseOrigin(k) =
        boost::mpi::all_reduce(comm, baseOrigin(k), std::plus<Real>()) / baseCount;
    }
    Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
    for (Index v : baseVertices)
    {
      const Vec dx = point3(mesh, v) - baseOrigin;
      covariance += dx * dx.transpose();
    }
    for (size_t i = 0; i < 3; ++i)
    {
      for (size_t j = 0; j < 3; ++j)
      {
        covariance(i, j) =
          boost::mpi::all_reduce(comm, Real(covariance(i, j)), std::plus<Real>());
      }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> plane(covariance);
    if (plane.info() != Eigen::Success || plane.eigenvalues()(1) <= 0)
    {
      Alert::Exception() << "The basal rim does not define a plane" << Alert::Raise;
    }
    waveDirection = plane.eigenvectors().col(0);
    size_t vertexCount = 0;
    for (auto v = mesh.getVertex(); v; ++v)
    {
      if (shard.isOwned(0, v->getIndex()))
      {
        wallCenter += point3(mesh, v->getIndex());
        ++vertexCount;
      }
    }
    vertexCount = boost::mpi::all_reduce(comm, vertexCount, std::plus<size_t>());
    for (size_t k = 0; k < 3; ++k)
    {
      wallCenter(k) =
        boost::mpi::all_reduce(comm, wallCenter(k), std::plus<Real>()) / vertexCount;
    }
    if (waveDirection.dot(wallCenter - baseOrigin) < 0)
      waveDirection *= -1;
  }

  void ProblemData::generateTransmuralCoordinate()
  {
    constexpr Attribute Endocardium = 8, Epicardium = 9;
    const auto& shard = mesh.getShard();
    const auto& comm = mesh.getContext().getCommunicator();
    std::set<std::array<Index, 3>> inner, outer;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      std::array<Index, 3> a, b;
      for (size_t j = 0; j < 3; ++j)
      {
        a[j] = cell->getVertices()(j);
        b[j] = cell->getVertices()(3 + j);
      }
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      inner.insert(a);
      outer.insert(b);
    }
    bool reserved = false;
    for (auto face = mesh.getFace(); face; ++face)
      reserved = reserved || face->getAttribute() == Endocardium ||
        face->getAttribute() == Epicardium;
    if (boost::mpi::all_reduce(comm, reserved, std::logical_or<bool>()))
      Alert::Exception() << "Boundary labels 8 and 9 are reserved for coordinate solves"
                         << Alert::Raise;
    std::vector<std::pair<Index, Optional<Attribute>>> labels;
    for (auto face = mesh.getFace(); face; ++face)
    {
      if (!shard.isBoundary(face->getIndex()) ||
        face->getGeometry() != Polytope::Type::Triangle)
        continue;
      std::array<Index, 3> key;
      for (size_t j = 0; j < 3; ++j)
        key[j] = face->getVertices()(j);
      std::sort(key.begin(), key.end());
      if (!inner.contains(key) && !outer.contains(key))
        continue;
      labels.emplace_back(face->getIndex(), face->getAttribute());
      mesh.setAttribute({2, face->getIndex()}, inner.contains(key) ? Endocardium : Epicardium);
    }
    m_laplace = Integral(Grad(m_coordinate), Grad(m_coordinateTest)) +
      DirichletBC(m_coordinate, RealFunction(0)).on(Endocardium) +
      DirichletBC(m_coordinate, RealFunction(1)).on(Epicardium);
    m_coordinateSolver.setPrefix(std::string("atrium_transmural_"))
      .setType(KSPCG).setTolerances(1e-12, 1e-14, 1e5, 10000);
    m_laplace.solve(m_coordinateSolver);
    transmuralResidual = m_coordinateSolver.getError();
    for (const auto& item : labels)
      mesh.setAttribute({2, item.first}, item.second);
    if (!std::isfinite(transmuralResidual) || transmuralResidual > 1e-8)
      Alert::Exception() << "The transmural-coordinate Laplace solve failed: "
                         << transmuralResidual << Alert::Raise;
    transmuralCoordinate = m_coordinate.getSolution();
    transmuralCoordinate.setName("TransmuralCoordinate");
    m_coordinate.getSolution() = Real(0);
  }

  ProblemData::Frame ProblemData::frameAt(const Point& point) const
  {
    Frame frame;
    const auto transverse = Grad(transmuralCoordinate).getValue(point);
    frame.normal = Vec(transverse(0), transverse(1), transverse(2));
    if (!frame.normal.allFinite() || !(frame.normal.norm() > 0))
      Alert::Exception() << "The transmural coordinate has an invalid gradient"
                         << Alert::Raise;
    frame.normal.normalize();
    const auto gradient = Grad(fiberCoordinate).getValue(point);
    frame.longitudinal = Vec(gradient(0), gradient(1), gradient(2));
    constexpr Real tolerance = 1e-10;
    if (frame.longitudinal.norm() > tolerance /
      std::max(cfg.apexRadius, (cfg.appendageApex - baseOrigin).norm()))
      frame.longitudinal.normalize();
    else
      frame.longitudinal.setZero();
    frame.longitudinal -= frame.longitudinal.dot(frame.normal) * frame.normal;
    if (frame.longitudinal.norm() < tolerance)
    {
      frame.fallback = 1;
      const auto& x = point.getPhysicalCoordinates();
      frame.longitudinal = cfg.appendageApex - Vec(x(0), x(1), x(2));
      if (frame.longitudinal.norm() > 0)
        frame.longitudinal.normalize();
      frame.longitudinal -= frame.longitudinal.dot(frame.normal) * frame.normal;
    }
    if (frame.longitudinal.norm() < tolerance)
    {
      frame.fallback = 2;
      Eigen::Index axis;
      frame.normal.cwiseAbs().minCoeff(&axis);
      frame.longitudinal = Vec::Unit(axis) - frame.normal(axis) * frame.normal;
    }
    frame.longitudinal.normalize();
    frame.fiber = frame.longitudinal.cross(frame.normal).normalized();
    return frame;
  }

  void ProblemData::generateFibers()
  {
    const auto& comm = mesh.getContext().getCommunicator();
    const auto& shard = mesh.getShard();
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
    bool reservedLabel = false;
    for (auto face = mesh.getFace(); face; ++face)
    {
      if (!shard.isBoundary(face->getIndex()))
        continue;
      if (face->getAttribute() == ApexPatchLabel)
        reservedLabel = true;
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
    if (boost::mpi::all_reduce(comm, reservedLabel, std::logical_or<bool>()))
    {
      Alert::Exception() << "Boundary label 8 is reserved for the temporary apex patch"
                         << Alert::Raise;
    }
    for (const auto& item : restoredLabels)
    {
      if (shard.isOwned(2, item.first))
        ++apexPatchFaces;
    }
    apexPatchFaces = boost::mpi::all_reduce(comm, apexPatchFaces, std::plus<size_t>());
    if (apexPatchFaces == 0)
    {
      Alert::Exception() << "No wall triangles lie in the appendage-apex patch; check "
                            "--apex and --apex-radius"
                         << Alert::Raise;
    }
    bool overlaps = false;
    for (Index vertex : apexVertices)
      overlaps = overlaps || baseAllVertices.contains(vertex);
    if (boost::mpi::all_reduce(comm, overlaps, std::logical_or<bool>()))
      Alert::Exception() << "The appendage patch overlaps the base rim" << Alert::Raise;
    for (const auto& item : restoredLabels)
      mesh.setAttribute({2, item.first}, ApexPatchLabel);
    m_laplace = Integral(Grad(m_coordinate), Grad(m_coordinateTest)) +
      DirichletBC(m_coordinate, RealFunction(0)).on(cfg.baseLabel) +
      DirichletBC(m_coordinate, RealFunction(1)).on(ApexPatchLabel);
    m_coordinateSolver.setPrefix(std::string("atrium_fibers_"))
      .setType(KSPCG)
      .setTolerances(1e-12, 1e-14, 1e5, 10000);
    m_laplace.solve(m_coordinateSolver);
    coordinateResidual = m_coordinateSolver.getError();
    if (!std::isfinite(coordinateResidual) || coordinateResidual > 1e-8)
    {
      Alert::Exception() << "The fiber-coordinate Laplace solve failed: "
                         << coordinateResidual << Alert::Raise;
    }
    for (const auto& item : restoredLabels)
      mesh.setAttribute({2, item.first}, item.second);
    fiberCoordinate = m_coordinate.getSolution();
    fiberCoordinate.setName("FiberCoordinate");
    const auto& faceRule = QF::QuadratureFormula::get(4, Polytope::Type::Triangle);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      const auto frame = frameAt(point);
      normals[cell->getIndex()] = frame.normal;
      longitudinal[cell->getIndex()] = frame.longitudinal;
      directions[cell->getIndex()] = frame.fiber;
      fallback[cell->getIndex()] = frame.fallback;
      if (!shard.isOwned(3, cell->getIndex()))
        continue;
      fallbackCells += frame.fallback > 0;
      for (size_t layer = 0; layer < 2; ++layer)
      {
        for (size_t q = 0; q < faceRule.getSize(); ++q)
        {
          const auto& r = faceRule.getPoint(q);
          const Point sample(*cell, Math::SpatialVector<Real>{r(0), r(1), Real(layer)});
          const auto local = frameAt(sample);
          const auto& jacobian = sample.getJacobian();
          const Vec a(jacobian(0, 0), jacobian(1, 0), jacobian(2, 0));
          const Vec b(jacobian(0, 1), jacobian(1, 1), jacobian(2, 1));
          const Real error = std::abs(local.fiber.dot(a.cross(b).normalized()));
          auto& maximum = layer == 0 ? endocardialTangencyError : epicardialTangencyError;
          maximum = std::max(maximum, error);
        }
      }
    }
    endocardialTangencyError = boost::mpi::all_reduce(
      comm, endocardialTangencyError, boost::mpi::maximum<Real>());
    epicardialTangencyError = boost::mpi::all_reduce(
      comm, epicardialTangencyError, boost::mpi::maximum<Real>());
    if (std::max(endocardialTangencyError, epicardialTangencyError) > 1e-8)
      Alert::Exception() << "The fiber frame is not tangent to both wall surfaces"
                         << Alert::Raise;
    fallbackCells = boost::mpi::all_reduce(comm, fallbackCells, std::plus<size_t>());
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
    fiber = frameAt(point).fiber;
    return fiber;
  }
  Real ProblemData::electricalAt(const Point& point, Real time) const
  {
    return (cfg.activation + cfg.relaxation) * activationAt(point, time) - cfg.relaxation;
  }

  Optional<Real> ProblemData::minimumDeformedJacobian(
    const VectorField& displacement) const
  {
    Real minimum = std::numeric_limits<Real>::infinity();
    const auto& comm = mesh.getContext().getCommunicator();
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (!mesh.getShard().isOwned(3, cell->getIndex()))
        continue;
      std::array<Vec, 6> q;
      for (size_t j = 0; j < 6; ++j)
      {
        const Index vertex = cell->getVertices()(j);
        const auto x = mesh.getVertexCoordinates(vertex);
        const auto& dofs = displacement.getFiniteElementSpace().getDOFs(0, vertex);
        for (size_t k = 0; k < 3; ++k)
          q[j](k) = x(k) + displacement[dofs(k)];
      }
      const Real determinant = minimumJacobian(q);
      const Real area = (q[1] - q[0]).cross(q[2] - q[0]).norm();
      if (!std::isfinite(determinant) || !std::isfinite(area) || area <= 0)
        minimum = -std::numeric_limits<Real>::infinity();
      else
        minimum = std::min(minimum, determinant / area);
    }
    minimum = boost::mpi::all_reduce(comm, minimum, boost::mpi::minimum<Real>());
    if (!std::isfinite(minimum))
      return std::nullopt;
    return minimum;
  }

}
