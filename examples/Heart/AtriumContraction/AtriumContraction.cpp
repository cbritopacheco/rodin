/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @brief Quasistatic atrial wall contraction driven by a prescribed plane wave.
 *
 * A NeoHookean wall is augmented by ActiveContraction's static fiber branch.
 * The prescribed active extension follows a plane moving inward from basal label 6.
 * This is a demonstration with prescribed activation and target fiber strain;
 * it does not solve electrophysiology or evolve Hill-Maxwell internal states.
 * Coordinates and displacement are in mm, time in s, and stresses in kPa.
 * Simple circumferential fibers follow a harmonic base-to-appendage coordinate.
 * The inlet/outlet rim faces (labels 1,2,3,5,6) are fixed. Wall faces are free.
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <exception>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <Eigen/Eigenvalues>

#include <Rodin/Alert.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  using Vec = Eigen::Matrix<Real, 3, 1>;

  struct Options
  {
    std::string mesh =
      "resources/examples/Heart/Atrium/AtriumWall_Repaired_Wedges_1mm.mesh";
    std::string output = "output/atrium-contraction/Contraction";
    size_t steps = 24;
    size_t stopAfter = 0;
    Attribute baseLabel = 6;
    Vec appendageApex = Vec(46.39638367, 187.44307667, -205.03317767);
    Real apexRadius = 3;
    bool fibersOnly = false;
    Real contraction = 0.03;
    Real speed = 200;
    Real width = 6;
    Real young = 50;
    Real poisson = 0.3;
    Real activeStiffness = 20;
  };

  Options options(int argc, char** argv)
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
      !std::isfinite(out.apexRadius) || out.apexRadius <= 0)
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
}

int main(int argc, char** argv)
{
  try
  {
    const auto cfg = options(argc, argv);
    Mesh mesh;
    mesh.load(cfg.mesh, IO::FileFormat::MEDIT);
    if (mesh.getDimension() != 3 || mesh.getSpaceDimension() != 3)
    {
      Alert::Exception() << "Expected a 3D wall volume mesh" << Alert::Raise;
    }
    mesh.getConnectivity().compute(2, 3);
    std::vector<Vec> directions(mesh.getCellCount()), normals(mesh.getCellCount());
    std::vector<bool> innerVertex(mesh.getVertexCount(), false);
    std::set<Index> baseVertices, baseAllVertices;
    size_t badCells = 0, rimFaces = 0;
    Real minimumReferenceJacobian = std::numeric_limits<Real>::infinity();
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
    Vec baseOrigin = Vec::Zero(), wallCenter = Vec::Zero();
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
    Vec waveDirection = plane.eigenvectors().col(0);
    for (auto v = mesh.getVertex(); v; ++v)
      wallCenter += point3(mesh, v->getIndex());
    wallCenter /= mesh.getVertexCount();
    if (waveDirection.dot(wallCenter - baseOrigin) < 0)
      waveDirection *= -1;
    const auto waveCoordinateAt = [&](const Point& point) {
      const auto x = point.getPhysicalCoordinates();
      return waveDirection.dot(Vec(x(0), x(1), x(2)) - baseOrigin);
    };
    const auto parent = std::filesystem::path(cfg.output).parent_path();
    if (!parent.empty())
      std::filesystem::create_directories(parent);

    P1 vectorSpace(mesh, 3);
    P1 scalarSpace(mesh);
    P0<Math::SpatialVector<Real>> fiberSpace(mesh, 3);
    const VectorFunction zero(3, [](const Point&) {
      Math::SpatialVector<Real> result(3);
      result.setZero();
      return result;
    });
    GridFunction displacement(vectorSpace);
    GridFunction activation(scalarSpace);
    GridFunction extension(scalarSpace);
    displacement = zero;
    displacement.setName("Displacement");
    activation.setName("Activation");
    extension.setName("ActiveExtension");
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
    const Real coordinateResidual =
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
    GridFunction fiberCoordinate(scalarSpace);
    fiberCoordinate = coordinate.getSolution();
    fiberCoordinate.setName("FiberCoordinate");
    const auto coordinateGradient = Grad(fiberCoordinate);
    std::vector<Vec> longitudinal(mesh.getCellCount());
    std::vector<Real> fallback(mesh.getCellCount(), 0);
    size_t fallbackCells = 0;
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
    GridFunction fiberField(fiberSpace), normalField(fiberSpace),
      longitudinalField(fiberSpace);
    fiberField.setName("FiberDirection");
    fiberField = vectorField(directions);
    normalField.setName("WallNormal");
    normalField = vectorField(normals);
    longitudinalField.setName("LongitudinalDirection");
    longitudinalField = vectorField(longitudinal);
    P0 diagnosticSpace(mesh);
    GridFunction fallbackField(diagnosticSpace);
    fallbackField.setName("FiberFallback");
    fallbackField = [&](const Point& point) {
      return fallback.at(point.getPolytope().getIndex());
    };
    GridFunction waveCoordinate(scalarSpace);
    waveCoordinate.setName("WaveCoordinate");
    waveCoordinate = waveCoordinateAt;
    const Real minimumProjection = waveCoordinate.getData().minCoeff();
    const Real maximumProjection = waveCoordinate.getData().maxCoeff();
    constexpr Real PulsePadding = 3;
    const Real endTime = (maximumProjection + PulsePadding * cfg.width) / cfg.speed;
    const Real dt = endTime / cfg.steps;
    Real time = 0;
    const auto activate = [&](const Point& point) {
      const Real argument = (waveCoordinateAt(point) - cfg.speed * time) / cfg.width;
      return std::exp(-argument * argument);
    };
    const auto updateActivation = [&] {
      activation = [&](const Point& point) { return activate(point); };
      extension = [&](const Point& point) { return -cfg.contraction * activate(point); };
    };
    std::ofstream metadata(cfg.output + ".geometry.json");
    metadata << std::setprecision(17) << "{\n  \"base_label\": " << cfg.baseLabel
             << ",\n  \"base_origin_mm\": [" << baseOrigin(0) << ',' << baseOrigin(1)
             << ',' << baseOrigin(2) << "],\n  \"wave_direction\": [" << waveDirection(0)
             << ',' << waveDirection(1) << ',' << waveDirection(2)
             << "],\n  \"appendage_apex_mm\": [" << cfg.appendageApex(0) << ','
             << cfg.appendageApex(1) << ',' << cfg.appendageApex(2)
             << "],\n  \"apex_radius_mm\": " << cfg.apexRadius
             << ",\n  \"apex_patch_faces\": " << restoredLabels.size()
             << ",\n  \"fiber_fallback_cells\": " << fallbackCells
             << ",\n  \"fiber_coordinate_residual\": " << coordinateResidual
             << ",\n  \"minimum_wave_coordinate_mm\": " << minimumProjection
             << ",\n  \"maximum_wave_coordinate_mm\": " << maximumProjection
             << ",\n  \"speed_mm_per_s\": " << cfg.speed
             << ",\n  \"width_mm\": " << cfg.width << "\n}\n";
    metadata.close();
    std::cout << "Base label " << cfg.baseLabel << ": origin=" << baseOrigin.transpose()
              << "; inward normal=" << waveDirection.transpose()
              << "; appendage projection="
              << waveDirection.dot(cfg.appendageApex - baseOrigin)
              << " mm; harmonic residual=" << coordinateResidual
              << "; fiber fallback cells=" << fallbackCells << std::endl;
    const Real mu = cfg.young / (2 * (1 + cfg.poisson));
    const Real lambda =
      cfg.young * cfg.poisson / ((1 + cfg.poisson) * (1 - 2 * cfg.poisson));
    Solid::ActiveFiberLaw::Parameters parameters;
    parameters.stiffness = cfg.activeStiffness;
    Solid::ActiveContraction<Solid::NeoHookean, Solid::ActiveFiberLaw> law(
      Solid::NeoHookean(lambda, mu), Solid::ActiveFiberLaw(parameters));
    const auto input = [&](Solid::ConstitutivePoint& cp) {
      const Point& point = cp.getPoint()->get();
      Math::SpatialVector<Real> fiber(3);
      fiber = directions.at(point.getPolytope().getIndex());
      cp.set<Solid::Tags::FiberDirection>(fiber);
      cp.set<Solid::Tags::ActiveExtension>(-cfg.contraction * activate(point));
    };
    TrialFunction increment(vectorSpace);
    TestFunction test(vectorSpace);
    auto work = Solid::InternalVirtualWork(law, displacement).setInput(input);
    // A quadrature policy for this non-polynomial hyperelastic example.
    constexpr size_t QuadratureOrder = 4;
    work.setQuadratureOrder(QuadratureOrder);
    Problem mechanics(increment, test);
    mechanics = work(increment, test) +
      DirichletBC(increment, zero)
        .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
    Solver::SparseLU linearSolver(mechanics);
    Solver::NewtonSolver newton(linearSolver);
    newton.setMaxIterations(30).setAbsoluteTolerance(1e-8).setRelativeTolerance(1e-7);
    newton.setMonitor([](const auto& r) {
      std::cout << "  Newton " << r.iterations << " residual=" << r.finalResidual
                << std::endl;
    });
    IO::XDMF output(cfg.output);
    auto grid = output.grid();
    grid.setMesh(mesh);
    grid.add(displacement);
    grid.add(activation);
    grid.add(extension);
    grid.add(fiberField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(normalField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(
      longitudinalField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(fallbackField, IO::XDMF::Center::Cell, IO::XDMF::AttributePolicy::Static);
    grid.add(fiberCoordinate, IO::XDMF::Center::Node, IO::XDMF::AttributePolicy::Static);
    grid.add(waveCoordinate, IO::XDMF::Center::Node, IO::XDMF::AttributePolicy::Static);
    updateActivation();
    if (cfg.fibersOnly)
    {
      output.write(0).close();
      std::cout << "Fiber inspection only: no mechanical equilibrium solve" << std::endl;
      return 0;
    }
    std::ofstream log(cfg.output + ".csv");
    log << "time_s,max_displacement_mm,activation_max,newton_iterations,residual,min_"
           "scaled_jacobian_mm\n";
    log.flush();
    const size_t count = cfg.stopAfter ? std::min(cfg.steps, cfg.stopAfter) : cfg.steps;
    std::cout << std::setprecision(9) << "Wall: " << mesh.getCellCount() << " wedges, "
              << mesh.getVertexCount() << " vertices; fixed rim faces=" << rimFaces
              << "; minimum reference scaled Jacobian=" << minimumReferenceJacobian
              << "; dt=" << dt << "; end=" << endTime << std::endl;
    for (size_t step = 0; step <= count; ++step)
    {
      time = step * dt;
      updateActivation();
      std::cout << "Step " << step << '/' << cfg.steps << " time=" << time
                << " front_distance_from_base=" << cfg.speed * time << std::endl;
      newton.solve(displacement);
      if (!newton.getReport().converged)
      {
        Alert::Exception() << "Contraction Newton solve did not converge" << Alert::Raise;
      }
      Real maximum = 0;
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        const auto& dofs = vectorSpace.getDOFs(0, vertex->getIndex());
        Vec u;
        for (size_t j = 0; j < 3; ++j)
          u(j) = displacement.getData()(dofs(j));
        maximum = std::max(maximum, u.norm());
      }
      Real minimumCurrentJacobian = std::numeric_limits<Real>::infinity();
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        std::array<Vec, 6> q;
        for (size_t j = 0; j < 6; ++j)
        {
          const Index vertex = cell->getVertices()(j);
          const auto x = mesh.getVertexCoordinates(vertex);
          const auto& dofs = vectorSpace.getDOFs(0, vertex);
          for (size_t k = 0; k < 3; ++k)
            q[j](k) = x(k) + displacement.getData()(dofs(k));
        }
        const Real determinant = minimumJacobian(q);
        if (!(determinant > 0) || !std::isfinite(determinant))
        {
          Alert::Exception() << "The contraction produced a folded wedge" << Alert::Raise;
        }
        const Real area = (q[1] - q[0]).cross(q[2] - q[0]).norm();
        minimumCurrentJacobian = std::min(minimumCurrentJacobian, determinant / area);
      }
      output.write(time).flush();
      log << std::setprecision(17) << time << ',' << maximum << ','
          << activation.getData().maxCoeff() << ',' << newton.getReport().iterations
          << ',' << newton.getReport().finalResidual << ',' << minimumCurrentJacobian
          << '\n';
      log.flush();
      std::cout << "  max displacement=" << maximum << " mm" << std::endl;
    }
    output.close();
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  return 0;
}
