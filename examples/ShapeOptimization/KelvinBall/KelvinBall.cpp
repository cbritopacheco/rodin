/* Kelvin-ball shape optimisation in a 1/24 cubic chamber. */
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string_view>
#include <stdexcept>
#include <vector>

#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Alert/Info.h>
#include <Rodin/Alert/Notation.h>
#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MMG.h>
#include <Rodin/Variational.h>

#include "RotatedNitsche.h"
#include "SewedOutput.h"
#include "Metrics.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  constexpr Attribute Obstacle = 2;
  constexpr Attribute Fluid = 3;
  constexpr Attribute Outer = 7;
  constexpr Attribute Gamma = 13;
  constexpr Attribute SigmaPlus = 31;
  constexpr Attribute SigmaMinus = 32;
  constexpr Attribute SigmaXYPlus = 33;
  constexpr Attribute SigmaXYMinus = 35;

  constexpr Real mu = 1.0;
  constexpr Real remeshGradation = 2.0;
  constexpr Real periodicPhysicalTolerance = 0.01;
  constexpr Real periodicReferenceTolerance = 0.25;
  constexpr Real defaultNitschePenalty = 320;
  constexpr Real linearResidualTolerance = 1e-8;
  constexpr size_t chamberMultiplicity = 24;
#ifdef RODIN_USE_OPENMP
  constexpr const char* assemblyBackend = "OpenMP";
#else
  constexpr const char* assemblyBackend = "Sequential";
#endif

  Alert::Text<Alert::CyanT> stageHeading(const std::string& text)
  {
    Alert::Text<Alert::CyanT> heading(Alert::Cyan, text);
    return heading.setBold();
  }

  Alert::Text<Alert::YellowT> substageHeading(const std::string& text)
  {
    Alert::Text<Alert::YellowT> heading(Alert::Yellow, text);
    return heading.setBold();
  }

  std::string diagnosticLabel(const std::string& text)
  {
    constexpr size_t width = 36;
    return text + std::string(text.size() < width ? width - text.size() : 1, ' ');
  }

  struct GradientDiagnostics
  {
      Real residual;
      Real jump;
  };

  struct MeshDiagnostics
  {
      size_t vertices = 0;
      size_t cells = 0;
      size_t obstacleCells = 0;
      size_t fluidCells = 0;
      size_t interfaceTriangles = 0;
  };

  struct ReconstructionDiagnostics
  {
      Real minimumSize = 0;
      Real maximumSize = 0;
      Real hausdorffTolerance = 0;
      size_t requiredBoundaryTriangles = 0;
      size_t cellsBefore = 0;
      size_t cellsAfter = 0;
  };

  void printUsage(const char* executable)
  {
    Alert::Info() << "Usage" << Alert::NewLine << "  " << executable << " [options]"
                  << Alert::NewLine << Alert::Notation("--n=<points>")
                  << "              Background points per edge (default: 13)."
                  << Alert::NewLine << Alert::Notation("--h=<size>")
                  << "                Requested background size; alternative to --n."
                  << Alert::NewLine << Alert::Notation("--iterations=<count>")
                  << "       Number of evaluated designs (default: 1)." << Alert::NewLine
                  << Alert::Notation("--penalty=<value>")
                  << "          Rotational Nitsche penalty (default: 320)."
                  << Alert::NewLine << Alert::Notation("--stabilization=<value>")
                  << "    P1--P1 pressure-stabilization factor (default: 0.05)."
                  << Alert::NewLine << Alert::Notation("--regularization=<value>")
                  << "   H1 smoothing length in multiples of h (default: 4)."
                  << Alert::NewLine << Alert::Notation("--step=<value>")
                  << "             Advection time step in multiples of h (default: 0.1)."
                  << Alert::NewLine << Alert::Notation("--geometry-only")
                  << "             Stop after initial MMG reconstruction."
                  << Alert::NewLine << Alert::Notation("--state-only")
                  << "                Stop after the first Stokes evaluation."
                  << Alert::NewLine << Alert::Notation("--save-mesh")
                  << "                 Save KelvinBallInitial.mesh." << Alert::NewLine
                  << Alert::Notation("--help")
                  << "                      Show this message." << Alert::Raise;
  }

  void announce(const char* stage)
  {
    Alert::Info() << stageHeading(stage) << Alert::Raise;
  }

  using LocalMesh = Geometry::Mesh<Context::Local>;
  using VelocitySpace = P1<Math::SpatialVector<Real>, LocalMesh>;
  using PressureSpace = P1<Real, LocalMesh>;

  VelocitySpace makeVelocitySpace(const LocalMesh& mesh)
  {
    return VelocitySpace(mesh, 3);
  }

  Math::SpatialMatrix<Real> rotation(std::initializer_list<Real> values)
  {
    Math::SpatialMatrix<Real> R(3, 3);
    auto value = values.begin();
    for (size_t i = 0; i < 3; ++i)
      for (size_t j = 0; j < 3; ++j)
        R(i, j) = *value++;
    return R;
  }

  const std::array<KelvinBall::RotationPair, 2> pairs{
    {{SigmaXYMinus, SigmaXYPlus, rotation({0, 1, 0, 1, 0, 0, 0, 0, -1})},
      {SigmaMinus, SigmaPlus, rotation({1, 0, 0, 0, 0, -1, 0, 1, 0})}}};
  const FlatSet<Attribute> masterCuts{SigmaPlus, SigmaXYPlus};

  Math::SpatialPoint centroid(const LocalMesh& mesh, const Polytope& face)
  {
    Math::SpatialPoint c(3);
    c.setZero();
    for (const Index vertex : face.getVertices())
      c += mesh.getVertexCoordinates(vertex);
    return c / static_cast<Real>(face.getVertices().size());
  }

  void splitSelfPairedCut(LocalMesh& mesh)
  {
    const size_t faceDimension = mesh.getDimension() - 1;
    for (auto face = mesh.getPolytope(faceDimension); face; ++face)
    {
      if (face->getAttribute() == SigmaXYPlus && centroid(mesh, *face).z() < 0)
        mesh.setAttribute({faceDimension, face->getIndex()}, SigmaXYMinus);
    }
  }

  Real planeResidual(Attribute attribute, const Math::SpatialPoint& x)
  {
    switch (attribute)
    {
      case Outer:
        return std::abs(x(0) - 2.0);
      case SigmaPlus:
        return std::abs(x(1) - x(2));
      case SigmaMinus:
        return std::abs(x(1) + x(2));
      case SigmaXYPlus:
      case SigmaXYMinus:
        return std::abs(x(0) - x(1));
      default:
        return 0;
    }
  }

  void checkFixedGeometry(const LocalMesh& mesh)
  {
    const FlatSet<Attribute> fixed{
      Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
    std::map<Attribute, size_t> counts;
    Real residual = 0;
    for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
    {
      const auto attribute = face->getAttribute();
      if (!attribute || !fixed.contains(*attribute))
        continue;
      ++counts[*attribute];
      for (const Index vertex : face->getVertices())
      {
        const auto& x = mesh.getVertexCoordinates(vertex);
        residual = std::max(residual, planeResidual(*attribute, x));
      }
    }
    for (const Attribute attribute : fixed)
    {
      if (counts[attribute] == 0)
        throw std::runtime_error("A fixed chamber boundary has no triangles.");
    }
    Alert::Info() << substageHeading("Mesh validation") << Alert::NewLine
                  << diagnosticLabel("Vertices:")
                  << Alert::Notation::Number(mesh.getVertexCount()) << Alert::NewLine
                  << diagnosticLabel("Cells:")
                  << Alert::Notation::Number(mesh.getCellCount()) << Alert::NewLine
                  << diagnosticLabel("Fixed-boundary planarity residual:")
                  << Alert::Notation::Number(residual) << Alert::Raise;
    if (residual > 1e-6)
      throw std::runtime_error(
        "Mesh reconstruction changed the fixed boundary geometry.");
  }

  MeshDiagnostics getMeshDiagnostics(const LocalMesh& mesh)
  {
    MeshDiagnostics diagnostics;
    diagnostics.vertices = mesh.getVertexCount();
    diagnostics.cells = mesh.getCellCount();
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (cell->getAttribute() == Obstacle)
        ++diagnostics.obstacleCells;
      else if (cell->getAttribute() == Fluid)
        ++diagnostics.fluidCells;
      else
        throw std::runtime_error("A cell has an unexpected material attribute.");
    }
    for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
    {
      if (face->getAttribute() == Gamma)
        ++diagnostics.interfaceTriangles;
    }
    if (diagnostics.obstacleCells == 0 || diagnostics.fluidCells == 0 ||
      diagnostics.interfaceTriangles == 0)
      throw std::runtime_error("The fitted body/fluid partition is incomplete.");
    return diagnostics;
  }

  void checkMaterials(const LocalMesh& mesh)
  {
    const auto diagnostics = getMeshDiagnostics(mesh);
    Alert::Info() << substageHeading("Material partition") << Alert::NewLine
                  << diagnosticLabel("Obstacle cells:")
                  << Alert::Notation::Number(diagnostics.obstacleCells) << Alert::NewLine
                  << diagnosticLabel("Fluid cells:")
                  << Alert::Notation::Number(diagnostics.fluidCells) << Alert::NewLine
                  << diagnosticLabel("Interface triangles:")
                  << Alert::Notation::Number(diagnostics.interfaceTriangles)
                  << Alert::Raise;
  }

  void reportSphereGeometry(const LocalMesh& mesh)
  {
    std::set<Index> vertices;
    for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
      if (face->getAttribute() == Gamma)
        vertices.insert(face->getVertices().begin(), face->getVertices().end());
    Real squared = 0;
    Real supremum = 0;
    for (const Index vertex : vertices)
    {
      const Real error = std::abs(mesh.getVertexCoordinates(vertex).norm() - 1.0);
      squared += error * error;
      supremum = std::max(supremum, error);
    }
    const Real rms = vertices.empty()
      ? std::numeric_limits<Real>::quiet_NaN()
      : std::sqrt(squared / static_cast<Real>(vertices.size()));
    Alert::Info() << substageHeading("Sphere geometry") << Alert::NewLine
                  << diagnosticLabel("Chamber volume:")
                  << Alert::Notation::Number(mesh.getVolume(Obstacle)) << Alert::NewLine
                  << diagnosticLabel("Vertex RMS error:") << Alert::Notation::Number(rms)
                  << Alert::NewLine << diagnosticLabel("Vertex maximum error:")
                  << Alert::Notation::Number(supremum) << Alert::Raise;
  }

  Math::SpatialVector<Real> boundingBoxSize(const LocalMesh& mesh)
  {
    Math::SpatialPoint minimum(3);
    Math::SpatialPoint maximum(3);
    minimum.setConstant(std::numeric_limits<Real>::infinity());
    maximum.setConstant(-std::numeric_limits<Real>::infinity());
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const auto& point = mesh.getVertexCoordinates(vertex);
      for (size_t component = 0; component < 3; ++component)
      {
        minimum(component) = std::min(minimum(component), point(component));
        maximum(component) = std::max(maximum(component), point(component));
      }
    }
    return maximum - minimum;
  }

  Math::SpatialVector<Real> sewnBoundingBoxSize(const LocalMesh& mesh)
  {
    Math::SpatialPoint minimum(3);
    Math::SpatialPoint maximum(3);
    minimum.setConstant(std::numeric_limits<Real>::infinity());
    maximum.setConstant(-std::numeric_limits<Real>::infinity());
    const auto rotations = KelvinBall::properOctahedralRotations();
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      for (const auto& rotation : rotations)
      {
        const Math::SpatialPoint point = rotation * mesh.getVertexCoordinates(vertex);
        for (size_t component = 0; component < 3; ++component)
        {
          minimum(component) = std::min(minimum(component), point(component));
          maximum(component) = std::max(maximum(component), point(component));
        }
      }
    }
    return maximum - minimum;
  }

  LocalMesh makeUniformChamber(size_t points)
  {
    if (points < 5)
      throw std::runtime_error("The chamber uniform grid requires at least five points.");
    LocalMesh cube =
      LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {points, points, points});
    const Real h = Real(2) / static_cast<Real>(points - 1);
    cube.scale(h);

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

    LocalMesh::Builder builder;
    builder.initialize(3).nodes(coordinates.size());
    for (const auto& x : coordinates)
      builder.vertex(x);
    for (auto& [vertices, attribute] : cells)
    {
      Index index;
      builder.polytope(Polytope::Type::Tetrahedron, std::move(vertices), index);
      builder.attribute({3, index}, attribute);
    }
    LocalMesh chamber = builder.finalize();
    chamber.getConnectivity().compute(2, 3);
    for (auto face = chamber.getBoundary(); face; ++face)
    {
      const auto c = centroid(chamber, *face);
      Attribute attribute = 0;
      if (std::abs(c(0) - Real(2)) < tolerance)
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

  size_t protectFixedGeometry(MMG::Mesh& mesh)
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

  template <class LevelSet>
  ReconstructionDiagnostics discretizeLevelSetMMG(
    MMG::Mesh& mesh, const LevelSet& levelSet, Real h, bool initial)
  {
    const size_t previousCells = mesh.getCellCount();
    const size_t requiredTriangles = protectFixedGeometry(mesh);
    const Real hmin = 0.8 * h;
    const Real hmax = 1.25 * h;
    const Real hausdorff = 0.1 * h * h;

    MMG::LevelSetDiscretizer discretizer;
    discretizer.split(Fluid, {Obstacle, Fluid})
      .setHMin(hmin)
      .setHMax(hmax)
      .setHausdorff(hausdorff)
      .setGradation(remeshGradation)
      .setBaseReferences(
        initial ? FlatSet<Attribute>{Fluid} : FlatSet<Attribute>{Obstacle, Fluid})
      .setBoundaryReference(Gamma)
      .setAngleDetection(true);
    if (!initial)
      discretizer.split(Obstacle, {Obstacle, Fluid});
    mesh = discretizer.discretize(levelSet);

    splitSelfPairedCut(mesh);
    Alert::Info() << substageHeading("MMG reconstruction") << Alert::NewLine
                  << diagnosticLabel("Minimum size:") << Alert::Notation::Number(hmin)
                  << Alert::NewLine << diagnosticLabel("Maximum size:")
                  << Alert::Notation::Number(hmax) << Alert::NewLine
                  << diagnosticLabel("Hausdorff tolerance:")
                  << Alert::Notation::Number(hausdorff) << Alert::NewLine
                  << diagnosticLabel("Required boundary triangles:")
                  << Alert::Notation::Number(requiredTriangles) << Alert::NewLine
                  << diagnosticLabel("Cell count:")
                  << Alert::Notation::Number(previousCells) << " -> "
                  << Alert::Notation::Number(mesh.getCellCount()) << Alert::Raise;
    checkFixedGeometry(mesh);
    checkMaterials(mesh);
    return {hmin, hmax, hausdorff, requiredTriangles, previousCells, mesh.getCellCount()};
  }

  ReconstructionDiagnostics initializeSphereMMG(MMG::Mesh& mesh, Real h)
  {
    P1 levelSetSpace(mesh);
    GridFunction sphere(levelSetSpace);
    sphere = RealFunction([](const Geometry::Point& point) {
      return point.getPhysicalCoordinates().norm() - Real(1);
    });
    return discretizeLevelSetMMG(mesh, sphere, h, true);
  }

  template <class ScalarGridFunction, class Locator>
  Real stitchScalarTrace(ScalarGridFunction& field, const Locator& locator)
  {
    const auto copy = field;
    const auto& space = field.getFiniteElementSpace();
    const auto& mesh = space.getMesh();
    Math::Vector<Real> sum(space.getSize());
    sum.setZero();
    std::vector<size_t> count(space.getSize(), 0);
    Real correction = 0;
    for (const KelvinBall::RotationPair& pair : pairs)
    {
      std::set<Index> vertices;
      for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
      {
        if (face->getAttribute() == pair.slave)
          vertices.insert(face->getVertices().begin(), face->getVertices().end());
      }
      for (const Index vertex : vertices)
      {
        const auto mapped =
          locator.locate(pair.master, pair.rotation * mesh.getVertexCoordinates(vertex));
        if (!mapped)
          throw std::runtime_error("A scalar trace node has no rotated master.");
        const Index dof = space.getDOFs(0, vertex)(0);
        sum(dof) += copy.getValue(*mapped);
        ++count[dof];
      }
    }
    for (Eigen::Index dof = 0; dof < sum.size(); ++dof)
    {
      if (count[dof] == 0)
        continue;
      const Real value = sum(dof) / count[dof];
      correction = std::max(correction, std::abs(field.getData()(dof) - value));
      field.getData()(dof) = value;
    }
    Alert::Info() << "Scalar trace stitching correction: "
                  << Alert::Notation::Number(correction) << Alert::Raise;
    return correction;
  }

  template <class Space, class Density, class Locator, class Output>
  GradientDiagnostics identifyGradient(const Space& shapeSpace, const Density& density,
    const Locator& locator, Output& gradient, Real regularizationLength,
    Real nitschePenalty, const char* name)
  {
    TrialFunction gradientTrial(shapeSpace);
    TestFunction testField(shapeSpace);
    auto interfaceNormal = FaceNormal(shapeSpace.getMesh());
    interfaceNormal.traceOf(Fluid);
    Problem identification(gradientTrial, testField);
    identification =
      Integral(regularizationLength * regularizationLength * Jacobian(gradientTrial),
        Jacobian(testField)) +
      Integral(gradientTrial, testField) -
      FaceIntegral(density, Dot(interfaceNormal, testField)).over(Gamma) +
      DirichletBC(gradientTrial, VectorFunction{0, 0, 0}).on(Outer);
    identification.assemble();
    auto& system = identification.getLinearSystem();
    KelvinBall::addRotatedVectorNitsche(shapeSpace, locator, pairs, system,
      regularizationLength * regularizationLength, nitschePenalty,
      FlatSet<Attribute>{Outer});
    KelvinBall::solveDirect(identification);
    const Real linearResidual =
      (system.getOperator() * system.getSolution() - system.getVector()).norm() /
      std::max(system.getVector().norm(), Real(1));
    if (!std::isfinite(linearResidual) || linearResidual > linearResidualTolerance)
      throw std::runtime_error("The gradient identification problem did not converge.");
    gradient = gradientTrial.getSolution();
    const Real jump = KelvinBall::rotatedVectorJump(locator, pairs, gradient);
    Alert::Info() << substageHeading(std::string(name) + " gradient") << Alert::NewLine
                  << diagnosticLabel("Linear residual:")
                  << Alert::Notation::Number(linearResidual) << Alert::NewLine
                  << diagnosticLabel("Rotated jump:") << Alert::Notation::Number(jump)
                  << Alert::Raise;
    return {linearResidual, jump};
  }

}

int main(int argc, char** argv)
{
  size_t points = 13;
  size_t maxIterations = 1;
  Real requestedH = 0;
  bool pointsSpecified = false;
  bool hSpecified = false;
  bool saveMeshDiagnostic = false;
  bool geometryOnly = false;
  bool stateOnly = false;
  Real nitschePenalty = defaultNitschePenalty;
  Real stabilizationFactor = 0.05;
  Real regularizationFactor = 4.0;
  Real stepFactor = 0.1;
  for (int argument = 1; argument < argc; ++argument)
  {
    const std::string_view mode(argv[argument]);
    if (mode.rfind("--n=", 0) == 0)
    {
      points = std::stoul(std::string(mode.substr(4)));
      pointsSpecified = true;
    }
    else if (mode.rfind("--h=", 0) == 0)
    {
      requestedH = std::stod(std::string(mode.substr(4)));
      hSpecified = true;
    }
    else if (mode.rfind("--iterations=", 0) == 0)
      maxIterations = std::stoul(std::string(mode.substr(13)));
    else if (mode.rfind("--penalty=", 0) == 0)
      nitschePenalty = std::stod(std::string(mode.substr(10)));
    else if (mode.rfind("--stabilization=", 0) == 0)
      stabilizationFactor = std::stod(std::string(mode.substr(16)));
    else if (mode.rfind("--regularization=", 0) == 0)
      regularizationFactor = std::stod(std::string(mode.substr(17)));
    else if (mode.rfind("--step=", 0) == 0)
      stepFactor = std::stod(std::string(mode.substr(7)));
    else if (mode == "--save-mesh")
      saveMeshDiagnostic = true;
    else if (mode == "--geometry-only")
      geometryOnly = true;
    else if (mode == "--state-only")
      stateOnly = true;
    else if (mode == "--help")
    {
      printUsage(argv[0]);
      return 0;
    }
    else
      throw std::runtime_error("Unknown KelvinBall option: " + std::string(mode));
  }
  if (pointsSpecified && hSpecified)
    throw std::runtime_error("Use either --n or --h, not both.");
  if (hSpecified)
  {
    if (!(requestedH > 0))
      throw std::runtime_error("The requested mesh size must be positive.");
    points = static_cast<size_t>(std::ceil(Real(2) / requestedH)) + 1;
  }
  if (points < 5)
    throw std::runtime_error("The chamber uniform grid requires --n >= 5.");
  if (maxIterations == 0)
    throw std::runtime_error("The iteration count must be positive.");
  if (!(nitschePenalty > 0))
    throw std::runtime_error("The Nitsche penalty must be positive.");
  if (stabilizationFactor < 0)
    throw std::runtime_error("The pressure stabilization must be nonnegative.");
  if (!(regularizationFactor > 0))
    throw std::runtime_error("The regularization factor must be positive.");
  if (!(stepFactor > 0))
    throw std::runtime_error("The advection step factor must be positive.");
  if (geometryOnly && stateOnly)
    throw std::runtime_error("Use either --geometry-only or --state-only, not both.");
  const Real h = Real(2) / static_cast<Real>(points - 1);
  const Real hilbertLength = regularizationFactor * h;
  const Real dt = stepFactor * h;
  Alert::Info() << substageHeading("Configuration") << Alert::NewLine
                << diagnosticLabel("Grid points:") << Alert::Notation::Number(points)
                << Alert::NewLine << diagnosticLabel("Effective h:")
                << Alert::Notation::Number(h) << Alert::NewLine
                << diagnosticLabel("Evaluated designs:")
                << Alert::Notation::Number(maxIterations) << Alert::NewLine
                << diagnosticLabel("Nitsche penalty:")
                << Alert::Notation::Number(nitschePenalty) << Alert::NewLine
                << diagnosticLabel("Stabilization factor:")
                << Alert::Notation::Number(stabilizationFactor) << Alert::NewLine
                << diagnosticLabel("Assembly backend:") << assemblyBackend
                << Alert::NewLine << diagnosticLabel("Direct solver:")
                << KelvinBall::DirectSolverName << Alert::NewLine
                << diagnosticLabel("Regularization length:")
                << Alert::Notation::Number(hilbertLength) << " = "
                << Alert::Notation::Number(regularizationFactor) << " h" << Alert::NewLine
                << diagnosticLabel("Advection step:") << Alert::Notation::Number(dt)
                << " = " << Alert::Notation::Number(stepFactor) << " h" << Alert::Raise;
  announce("Stage 1: Discretizing the initial sphere with MMG.");
  MMG::Mesh mesh(std::move(makeUniformChamber(points)));
  ReconstructionDiagnostics reconstruction = initializeSphereMMG(mesh, h);
  if (saveMeshDiagnostic)
    mesh.save("KelvinBallInitial.mesh", IO::FileFormat::MEDIT);
  reportSphereGeometry(mesh);
  if (geometryOnly)
    return 0;
  const Real targetVolume = mesh.getVolume(Obstacle);

  std::ofstream history("kelvin-ball.csv");
  history.precision(17);
  history << "iteration,h,dt,regularization_length,nitsche_penalty,"
             "stabilization_factor,assembly_backend,direct_solver,"
             "vertices,cells,obstacle_cells,fluid_cells,"
             "interface_triangles,remesh_hmin,remesh_hmax,remesh_hausdorff,"
             "required_boundary_triangles,remesh_cells_before,remesh_cells_after,"
             "k,c,q,rho,coupling_symmetry,nitsche_jump,volume,target_volume,"
             "volume_error,chamber_bbox_x,chamber_bbox_y,chamber_bbox_z,"
             "sewn_bbox_x,sewn_bbox_y,sewn_bbox_z,actual_delta_rho,"
             "predicted_delta_rho,actual_delta_volume,predicted_delta_volume,"
             "incoming_scalar_trace_jump,incoming_stitch_correction,"
             "null_space_multiplier,xi_rho_inf_norm,theta_inf_norm,d_rho_theta,"
             "d_volume_theta,required_d_volume_theta,rho_gradient_residual,"
             "rho_gradient_jump,volume_gradient_residual,"
             "volume_gradient_jump\n";
  IO::XDMF xdmf("KelvinBall");
  auto chamber = xdmf.grid("Chamber");
  chamber.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  auto fluidState = xdmf.grid("Fluid");
  IO::XDMF sewedXdmf("KelvinBallSewed");
  auto sewedDesignOutput = sewedXdmf.grid("Design");
  auto sewedFluidOutput = sewedXdmf.grid("Fluid");
  Optional<Real> previousRho;
  Optional<Real> predictedRhoChange;
  Optional<Real> previousVolume;
  Optional<Real> predictedVolumeChange;
  Optional<Real> incomingScalarTraceJump;
  Optional<Real> incomingStitchCorrection;

  for (size_t iteration = 0; iteration < maxIterations; ++iteration)
  {
    Alert::Info() << " ------------------------------------------------------------"
                  << Alert::NewLine << " Iteration "
                  << Alert::Notation::Number(iteration + 1) << " of "
                  << Alert::Notation::Number(maxIterations) << Alert::NewLine
                  << " ------------------------------------------------------------"
                  << Alert::Raise;
    announce("Stage 2: Preparing the chamber fluid mesh.");
    auto& connectivity = mesh.getConnectivity();
    connectivity.discover(3, 2);
    connectivity.discover(3, 1);
    connectivity.restrict(1, 0);
    connectivity.restrict(2, 0);
    connectivity.restrict(2, 3);
    connectivity.discover(0, 0);
    SubMesh fluid = mesh.trim(Obstacle);
    splitSelfPairedCut(fluid);
    auto& fluidConnectivity = fluid.getConnectivity();
    fluidConnectivity.discover(2, 1);
    fluidConnectivity.restrict(1, 0);
    fluidConnectivity.restrict(2, 3);

    VelocitySpace Vh = makeVelocitySpace(fluid);
    PressureSpace Qh(fluid);
    KelvinBall::AttributeFaceLocator fluidLocator(
      fluid, masterCuts, periodicPhysicalTolerance, periodicReferenceTolerance);
    GridFunction uT0(Vh), uT1(Vh), uT2(Vh), uR0(Vh), uR1(Vh), uR2(Vh);
    GridFunction pT0(Qh), pT1(Qh), pT2(Qh), pR0(Qh), pR1(Qh), pR2(Qh);

    announce("Stage 3: Solving the translational and rotational Stokes states.");
    const KelvinBall::Parameters metricParameters{h, nitschePenalty, stabilizationFactor};
    const KelvinBall::Values metrics = KelvinBall::evaluateChamber(Vh, Qh, fluidLocator,
      uT0, uT1, uT2, uR0, uR1, uR2, pT0, pT1, pT2, pR0, pR1, pR2, metricParameters);
    const Real k = metrics.k;
    const Real c = metrics.c;
    const Real q = metrics.q;
    const Real rho = metrics.rho;
    const Real couplingSymmetry = metrics.couplingSymmetry;
    const Real volume = mesh.getVolume(Obstacle);
    const Real nitscheJump = metrics.nitscheJump;
    const auto chamberBoundingBox = boundingBoxSize(mesh);
    const auto sewnBoundingBox = sewnBoundingBoxSize(mesh);
    const auto meshDiagnostics = getMeshDiagnostics(mesh);
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const Real actualRhoChange = previousRho ? rho - *previousRho : nan;
    const Real incomingPredictedRhoChange =
      predictedRhoChange ? *predictedRhoChange : nan;
    const Real actualVolumeChange = previousVolume ? volume - *previousVolume : nan;
    const Real incomingPredictedVolumeChange =
      predictedVolumeChange ? *predictedVolumeChange : nan;
    const Real scalarTraceJump = incomingScalarTraceJump ? *incomingScalarTraceJump : nan;
    const Real stitchCorrection =
      incomingStitchCorrection ? *incomingStitchCorrection : nan;
    const auto writeHistory = [&](Real nullSpaceMultiplier, Real xiRhoInfinityNorm,
                                Real thetaInfinityNorm, Real dRhoTheta,
                                Real dVolumeTheta, Real requiredDVolumeTheta,
                                const GradientDiagnostics& rhoGradientDiagnostics,
                                const GradientDiagnostics& volumeGradientDiagnostics) {
      history << iteration << ',' << h << ',' << dt << ',' << hilbertLength << ','
              << nitschePenalty << ',' << stabilizationFactor << ',' << assemblyBackend
              << ',' << KelvinBall::DirectSolverName << ',' << meshDiagnostics.vertices
              << ',' << meshDiagnostics.cells << ',' << meshDiagnostics.obstacleCells
              << ',' << meshDiagnostics.fluidCells << ','
              << meshDiagnostics.interfaceTriangles << ',' << reconstruction.minimumSize
              << ',' << reconstruction.maximumSize << ','
              << reconstruction.hausdorffTolerance << ','
              << reconstruction.requiredBoundaryTriangles << ','
              << reconstruction.cellsBefore << ',' << reconstruction.cellsAfter << ','
              << k << ',' << c << ',' << q << ',' << rho << ',' << couplingSymmetry << ','
              << nitscheJump << ',' << volume << ',' << targetVolume << ','
              << volume - targetVolume << ',' << chamberBoundingBox(0) << ','
              << chamberBoundingBox(1) << ',' << chamberBoundingBox(2) << ','
              << sewnBoundingBox(0) << ',' << sewnBoundingBox(1) << ','
              << sewnBoundingBox(2) << ',' << actualRhoChange << ','
              << incomingPredictedRhoChange << ',' << actualVolumeChange << ','
              << incomingPredictedVolumeChange << ',' << scalarTraceJump << ','
              << stitchCorrection << ',' << nullSpaceMultiplier << ','
              << xiRhoInfinityNorm << ',' << thetaInfinityNorm << ',' << dRhoTheta << ','
              << dVolumeTheta << ',' << requiredDVolumeTheta << ','
              << rhoGradientDiagnostics.residual << ',' << rhoGradientDiagnostics.jump
              << ',' << volumeGradientDiagnostics.residual << ','
              << volumeGradientDiagnostics.jump << '\n';
      history.flush();
    };
    announce("Stage 4: Evaluating the resistance objective.");
    if (previousRho)
    {
      Alert::Info() << substageHeading("Previous update") << Alert::NewLine
                    << diagnosticLabel("Actual delta rho:")
                    << Alert::Notation::Number(actualRhoChange) << Alert::NewLine
                    << diagnosticLabel("Predicted delta rho:")
                    << Alert::Notation::Number(incomingPredictedRhoChange)
                    << Alert::NewLine << diagnosticLabel("Actual delta volume:")
                    << Alert::Notation::Number(actualVolumeChange) << Alert::NewLine
                    << diagnosticLabel("Predicted delta volume:")
                    << Alert::Notation::Number(incomingPredictedVolumeChange)
                    << Alert::Raise;
    }
    Alert::Info() << substageHeading("Resistance metrics") << Alert::NewLine
                  << diagnosticLabel("Translational resistance k:")
                  << Alert::Notation::Number(k) << Alert::NewLine
                  << diagnosticLabel("Coupling coefficient c:")
                  << Alert::Notation::Number(c) << Alert::NewLine
                  << diagnosticLabel("Rotational resistance q:")
                  << Alert::Notation::Number(q) << Alert::NewLine
                  << diagnosticLabel("Objective rho:") << Alert::Notation::Number(rho)
                  << Alert::NewLine << diagnosticLabel("Volume:")
                  << Alert::Notation::Number(volume) << Alert::NewLine
                  << diagnosticLabel("Chamber bounding box:")
                  << Alert::Notation::Print(chamberBoundingBox.transpose())
                  << Alert::NewLine << diagnosticLabel("Sewn mesh bounding box:")
                  << Alert::Notation::Print(sewnBoundingBox.transpose()) << Alert::NewLine
                  << diagnosticLabel("Nitsche jump:")
                  << Alert::Notation::Number(nitscheJump) << Alert::NewLine
                  << diagnosticLabel("Coupling symmetry residual:")
                  << Alert::Notation::Number(couplingSymmetry) << Alert::Raise;
    if (stateOnly)
    {
      const GradientDiagnostics unavailable{nan, nan};
      writeHistory(nan, nan, nan, nan, nan, nan, unavailable, unavailable);
      return 0;
    }

    announce("Stage 5: Computing the shape derivative.");
    auto gradUT0 = Jacobian(uT0);
    gradUT0.traceOf(Fluid);
    auto gradUT1 = Jacobian(uT1);
    gradUT1.traceOf(Fluid);
    auto gradUT2 = Jacobian(uT2);
    gradUT2.traceOf(Fluid);
    auto gradUR0 = Jacobian(uR0);
    gradUR0.traceOf(Fluid);
    auto gradUR1 = Jacobian(uR1);
    gradUR1.traceOf(Fluid);
    auto gradUR2 = Jacobian(uR2);
    gradUR2.traceOf(Fluid);
    const auto strainUT0 = 0.5 * (gradUT0 + gradUT0.T());
    const auto strainUT1 = 0.5 * (gradUT1 + gradUT1.T());
    const auto strainUT2 = 0.5 * (gradUT2 + gradUT2.T());
    const auto strainUR0 = 0.5 * (gradUR0 + gradUR0.T());
    const auto strainUR1 = 0.5 * (gradUR1 + gradUR1.T());
    const auto strainUR2 = 0.5 * (gradUR2 + gradUR2.T());
    const Real signC = c < 0 ? -1.0 : 1.0;
    const auto kDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUT0, strainUT0) + Dot(strainUT1, strainUT1) +
        Dot(strainUT2, strainUT2));
    const auto cDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUT0, strainUR0) + Dot(strainUT1, strainUR1) +
        Dot(strainUT2, strainUR2));
    const auto qDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUR0, strainUR0) + Dot(strainUR1, strainUR1) +
        Dot(strainUR2, strainUR2));
    const auto rhoDensity = signC / std::sqrt(k * q) * cDensity -
      0.5 * rho * (kDensity / k + qDensity / q);

    P1 levelSetSpace(mesh);
    P1 shapeSpace(mesh, 3);
    KelvinBall::AttributeFaceLocator shapeCutLocator(mesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus},
      periodicPhysicalTolerance, periodicReferenceTolerance);
    GridFunction rhoGradient(shapeSpace), volumeGradient(shapeSpace);
    announce("Stage 6: Regularizing and constraining the shape direction.");
    const GradientDiagnostics rhoGradientDiagnostics = identifyGradient(shapeSpace,
      -rhoDensity, shapeCutLocator, rhoGradient, hilbertLength, nitschePenalty, "Rho");
    const GradientDiagnostics volumeGradientDiagnostics = identifyGradient(shapeSpace,
      RealFunction{-1}, shapeCutLocator, volumeGradient, hilbertLength,
      nitschePenalty, "Volume");
    TestFunction testField(shapeSpace);
    auto interfaceNormal = FaceNormal(mesh);
    interfaceNormal.traceOf(Fluid);
    LinearForm dRho(testField), dVolume(testField);
    dRho = FaceIntegral(-rhoDensity, Dot(interfaceNormal, testField)).over(Gamma);
    dVolume = FaceIntegral(RealFunction{-1}, Dot(interfaceNormal, testField)).over(Gamma);
    dRho.assemble();
    dVolume.assemble();
    const Real nullSpaceMultiplier = -dVolume(rhoGradient) / dVolume(volumeGradient);
    GridFunction xiRho(shapeSpace);
    xiRho = rhoGradient + nullSpaceMultiplier * volumeGradient;
    GridFunction xiRhoNorm(levelSetSpace);
    xiRhoNorm = Frobenius(xiRho);
    const Real xiRhoInfinityNorm = std::max(xiRhoNorm.max(), Real(1e-30));
    GridFunction theta(shapeSpace);
    theta = xiRho;
    theta /= xiRhoInfinityNorm;
    xiRhoNorm = Frobenius(theta);
    const Real thetaInfinityNorm = xiRhoNorm.max();
    const Real dRhoTheta = dRho(theta);
    const Real dVolumeTheta = dVolume(theta);
    const Real requiredDVolumeTheta = 0;
    Alert::Info() << substageHeading("Constrained direction") << Alert::NewLine
                  << diagnosticLabel("Null-space multiplier:")
                  << Alert::Notation::Number(nullSpaceMultiplier) << Alert::NewLine
                  << diagnosticLabel("Xi rho infinity norm:")
                  << Alert::Notation::Number(xiRhoInfinityNorm) << Alert::NewLine
                  << diagnosticLabel("Theta infinity norm:")
                  << Alert::Notation::Number(thetaInfinityNorm)
                  << Alert::NewLine << diagnosticLabel("D rho [theta]:")
                  << Alert::Notation::Number(dRhoTheta) << Alert::NewLine
                  << diagnosticLabel("D volume [theta]:")
                  << Alert::Notation::Number(dVolumeTheta) << Alert::NewLine
                  << diagnosticLabel("Required D volume [theta]:")
                  << Alert::Notation::Number(requiredDVolumeTheta) << Alert::Raise;
    writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
      dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
      volumeGradientDiagnostics);
    previousRho = rho;
    predictedRhoChange = dt * dRhoTheta;
    previousVolume = volume;
    predictedVolumeChange = dt * dVolumeTheta;
    GridFunction distance(levelSetSpace);
    Distance::Eikonal(distance).setInterior(Obstacle).setInterface(Gamma).solve().sign();

    announce("Stage 7: Writing the chamber and sewn fields.");
    chamber.clear();
    chamber.add("Distance", distance, IO::XDMF::Center::Node);
    chamber.add("Theta", theta, IO::XDMF::Center::Node);
    fluidState.clear();
    fluidState.setMesh(fluid, IO::XDMF::MeshPolicy::Transient);
    fluidState.add("Translation_0", uT0, IO::XDMF::Center::Node);
    fluidState.add("Translation_1", uT1, IO::XDMF::Center::Node);
    fluidState.add("Translation_2", uT2, IO::XDMF::Center::Node);
    fluidState.add("Rotation_0", uR0, IO::XDMF::Center::Node);
    fluidState.add("Rotation_1", uR1, IO::XDMF::Center::Node);
    fluidState.add("Rotation_2", uR2, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Translation_0", pT0, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Translation_1", pT1, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Translation_2", pT2, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Rotation_0", pR0, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Rotation_1", pR1, IO::XDMF::Center::Node);
    fluidState.add("Pressure_Rotation_2", pR2, IO::XDMF::Center::Node);

    auto sewedDesign = KelvinBall::sew(mesh, FlatSet<Attribute>{Gamma, Outer});
    P1 sewedDesignScalar(sewedDesign.mesh);
    P1 sewedDesignVector(sewedDesign.mesh, 3);
    GridFunction sewedDistance(sewedDesignScalar);
    GridFunction sewedVelocity(sewedDesignVector);
    KelvinBall::sewScalar(sewedDistance, distance, sewedDesign);
    KelvinBall::sewVector(sewedVelocity, theta, sewedDesign);
    sewedDesignOutput.clear();
    sewedDesignOutput.setMesh(sewedDesign.mesh, IO::XDMF::MeshPolicy::Transient);
    sewedDesignOutput.add("Distance", sewedDistance, IO::XDMF::Center::Node);
    sewedDesignOutput.add("Theta", sewedVelocity, IO::XDMF::Center::Node);

    auto sewedFluid = KelvinBall::sew(fluid, FlatSet<Attribute>{Gamma, Outer});
    VelocitySpace sewedVelocitySpace = makeVelocitySpace(sewedFluid.mesh);
    PressureSpace sewedPressureSpace(sewedFluid.mesh);
    GridFunction sewedUT0(sewedVelocitySpace), sewedUT1(sewedVelocitySpace),
      sewedUT2(sewedVelocitySpace), sewedUR0(sewedVelocitySpace),
      sewedUR1(sewedVelocitySpace), sewedUR2(sewedVelocitySpace);
    GridFunction sewedPT0(sewedPressureSpace), sewedPT1(sewedPressureSpace),
      sewedPT2(sewedPressureSpace), sewedPR0(sewedPressureSpace),
      sewedPR1(sewedPressureSpace), sewedPR2(sewedPressureSpace);
    const auto translations = std::array{&uT0, &uT1, &uT2};
    const auto rotations = std::array{&uR0, &uR1, &uR2};
    const auto translationPressures = std::array{&pT0, &pT1, &pT2};
    const auto rotationPressures = std::array{&pR0, &pR1, &pR2};
    KelvinBall::sewVectorLoad(sewedUT0, translations, 0, sewedFluid);
    KelvinBall::sewVectorLoad(sewedUT1, translations, 1, sewedFluid);
    KelvinBall::sewVectorLoad(sewedUT2, translations, 2, sewedFluid);
    KelvinBall::sewVectorLoad(sewedUR0, rotations, 0, sewedFluid);
    KelvinBall::sewVectorLoad(sewedUR1, rotations, 1, sewedFluid);
    KelvinBall::sewVectorLoad(sewedUR2, rotations, 2, sewedFluid);
    KelvinBall::sewScalarLoad(sewedPT0, translationPressures, 0, sewedFluid);
    KelvinBall::sewScalarLoad(sewedPT1, translationPressures, 1, sewedFluid);
    KelvinBall::sewScalarLoad(sewedPT2, translationPressures, 2, sewedFluid);
    KelvinBall::sewScalarLoad(sewedPR0, rotationPressures, 0, sewedFluid);
    KelvinBall::sewScalarLoad(sewedPR1, rotationPressures, 1, sewedFluid);
    KelvinBall::sewScalarLoad(sewedPR2, rotationPressures, 2, sewedFluid);
    sewedFluidOutput.clear();
    sewedFluidOutput.setMesh(sewedFluid.mesh, IO::XDMF::MeshPolicy::Transient);
    sewedFluidOutput.add("Translation_0", sewedUT0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Translation_1", sewedUT1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Translation_2", sewedUT2, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Rotation_0", sewedUR0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Rotation_1", sewedUR1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Rotation_2", sewedUR2, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Translation_0", sewedPT0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Translation_1", sewedPT1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Translation_2", sewedPT2, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Rotation_0", sewedPR0, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Rotation_1", sewedPR1, IO::XDMF::Center::Node);
    sewedFluidOutput.add("Pressure_Rotation_2", sewedPR2, IO::XDMF::Center::Node);

    xdmf.write(static_cast<Real>(iteration)).flush();
    sewedXdmf.write(static_cast<Real>(iteration)).flush();

    if (iteration + 1 == maxIterations)
      continue;
    announce("Stage 8: Advecting the level set.");
    TrialFunction advected(levelSetSpace);
    TestFunction test(levelSetSpace);
    KelvinBall::RotatedBoundaryPolicy periodicBoundary(
      -dt, mesh, shapeCutLocator, pairs);
    Problem transport(advected, test);
    transport = Integral(advected, test) -
      Integral(
        Flow(-dt, distance, theta, Math::RungeKutta::RK4{}, periodicBoundary), test);
    transport.assemble();
    KelvinBall::addRotatedScalarPenalty(
      levelSetSpace, shapeCutLocator, pairs, transport.getLinearSystem(), nitschePenalty);
    Solver::CG(transport).solve();
    incomingScalarTraceJump =
      KelvinBall::rotatedScalarJump(shapeCutLocator, pairs, advected.getSolution());
    Alert::Info() << "Weak scalar trace jump: "
                  << Alert::Notation::Number(*incomingScalarTraceJump) << Alert::Raise;
    incomingStitchCorrection = stitchScalarTrace(advected.getSolution(), shapeCutLocator);
    announce("Stage 9: Reconstructing the advected interface with MMG.");
    reconstruction = discretizeLevelSetMMG(mesh, advected.getSolution(), h, false);
  }

  Alert::Success() << "Wrote KelvinBall.xdmf, KelvinBallSewed.xdmf, and kelvin-ball.csv"
                   << Alert::Raise;
  return 0;
}
