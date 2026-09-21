/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "KelvinBallOptimization.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <stdexcept>
#include <vector>

#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Adaptation.h>
#include <Rodin/Alert/Info.h>
#include <Rodin/Location.h>
#include <Rodin/Alert/Notation.h>
#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MMG.h>
#include <Rodin/Variational.h>

#include "Configuration.h"
#include "Metrics.h"
#include "RotatedCharacteristicContinuation.h"
#include "RotatedNitscheIntegrator.h"
#include "SewedOutput.h"
#include "Sphere.h"
#include "../../WNGIRExampleParameters.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace KelvinBall
{
  class KelvinBallOptimization::Implementation
  {
    public:
      Implementation(int argc, char** argv)
        : m_argc(argc),
          m_argv(argv)
      {}

      int run();

    private:
      static constexpr Attribute Obstacle = KelvinBall::Obstacle;
      static constexpr Attribute Fluid = KelvinBall::Fluid;
      static constexpr Attribute Outer = KelvinBall::Outer;
      static constexpr Attribute Gamma = KelvinBall::Gamma;
      static constexpr Attribute SigmaPlus = KelvinBall::SigmaPlus;
      static constexpr Attribute SigmaMinus = KelvinBall::SigmaMinus;
      static constexpr Attribute SigmaXYPlus = KelvinBall::SigmaXYPlus;
      static constexpr Attribute SigmaXYMinus = KelvinBall::SigmaXYMinus;

      static constexpr Real mu = KelvinBall::Mu;
      static constexpr Real remeshGradation = 2.0;
      static constexpr Real rotatedTracePhysicalTolerance = 0.01;
      static constexpr Real rotatedTraceReferenceTolerance = 0.25;
      static constexpr Real linearResidualTolerance = KelvinBall::LinearResidualTolerance;
      static constexpr size_t chamberMultiplicity = KelvinBall::ChamberMultiplicity;
#ifdef RODIN_USE_OPENMP
      static constexpr const char* assemblyBackend = "OpenMP";
#else
      static constexpr const char* assemblyBackend = "Sequential";
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
          Real minimumQuality = std::numeric_limits<Real>::infinity();
          Real meanQuality = 0;
          Real maximumQuality = 0;
          Real meanElementSize = 0;
      };

      using ReconstructionDiagnostics = KelvinBall::ReconstructionDiagnostics;

      struct MMGReconstruction
      {
          MMG::Mesh mesh;
          ReconstructionDiagnostics diagnostics;
      };

      template <class Displacement>
      void moveMesh(KelvinBall::Mesh& moved, const KelvinBall::Mesh& mesh,
        const Displacement& displacement)
      {
        const auto& space = displacement.getFiniteElementSpace();
        const auto& coefficients = displacement.getData();
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          auto coordinates = mesh.getVertexCoordinates(vertex);
          const auto& dofs = space.getDOFs(0, vertex);
          for (size_t component = 0; component < mesh.getSpaceDimension(); ++component)
            coordinates(component) += coefficients(dofs[component]);
          moved.setVertexCoordinates(vertex, coordinates);
        }
      }

      template <class LevelSet>
      MMG::Mesh classifyLevelSetForWNGIR(
        const MMG::Mesh& background, const LevelSet& levelSet)
      {
        MMG::Mesh classified(background);
        const auto& space = levelSet.getFiniteElementSpace();
        const auto& coefficients = levelSet.getData();
        for (auto cell = classified.getCell(); cell; ++cell)
        {
          Real meanLevelSet = 0;
          for (const Index vertex : cell->getVertices())
          {
            const auto& dofs = space.getDOFs(0, vertex);
            meanLevelSet += coefficients(dofs[0]);
          }
          meanLevelSet /= static_cast<Real>(cell->getVertices().size());
          classified.setAttribute({classified.getDimension(), cell->getIndex()},
            meanLevelSet < 0 ? Obstacle : Fluid);
        }

        classified.getConnectivity().compute(2, 3);
        for (auto face = classified.getFace(); face; ++face)
        {
          const auto& incident =
            classified.getConnectivity().getIncidence({2, 3}, face->getIndex());
          if (incident.size() != 2)
            continue;
          const auto first = classified.getCell(incident[0])->getAttribute();
          const auto second = classified.getCell(incident[1])->getAttribute();
          if (first && second && *first != *second)
            classified.setAttribute({2, face->getIndex()}, Gamma);
        }
        return classified;
      }

      template <class LevelSet>
      MMGReconstruction fitLevelSetWNGIR(const KelvinBall::Mesh& mesh,
        const LevelSet& levelSet, Real h, Real outerRadius, int argc, char** argv)
      {
        P1<Math::SpatialVector<Real>, KelvinBall::Mesh> gradientSpace(mesh, 3);
        GridFunction projectedGradient(gradientSpace);
        TrialFunction gradientTrial(gradientSpace);
        TestFunction gradientTest(gradientSpace);
        Problem gradientProjection(gradientTrial, gradientTest);
        gradientProjection = Integral(gradientTrial, gradientTest)
                           - Integral(Grad(levelSet), gradientTest);
        gradientProjection.assemble();
        Solver::CG(gradientProjection).solve();
        projectedGradient.getData() = gradientTrial.getSolution().getData();

        P1<Math::SpatialVector<Real>, KelvinBall::Mesh> displacementSpace(mesh, 3);
        TrialFunction displacementTrial(displacementSpace);
        TestFunction displacementTest(displacementSpace);
        Rodin::Examples::WNGIRExampleDefaults defaults;
        defaults.kappaBulk = Real(8e-4);
        // A fit from the classified staircase converges in seven to nine steps.
        defaults.maxIterations = 12;
        auto parameters =
          Rodin::Examples::makeWNGIRParameters(argc, argv, h, Gamma, defaults);
        if (!Rodin::Examples::findOption(
              argc, argv, "wngir-primal-barrier-iterations", nullptr))
          parameters.primalBarrierIterations = 30;
        if (!Rodin::Examples::findOption(argc, argv, "wngir-cg-max-iters", nullptr))
          parameters.cgMaxIterations = 10000;
        // The outer sphere is curved, so its nodes are pinned: sliding in a
        // facet plane would walk them off the sphere. The cut planes bound the
        // wedge but are not walls, and the rim of the interface lies entirely
        // on them, so they slide within themselves instead of being frozen.
        parameters.fixedBoundaryAttributes = {Outer};
        parameters.slipBoundaryAttributes =
          {SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
        // Pinning the outer sphere removes every rigid mode, so no floor is
        // needed; the option still overrides this default.
        if (!Rodin::Examples::findOption(argc, argv, "wngir-rigid-stabilisation", nullptr))
          parameters.rigidStabilisationLevel = 0;
        Adaptation::WNGIR fitting(displacementTrial, displacementTest);
        fitting.setParameters(parameters);

        std::vector<Index> interface;
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          if (face->getAttribute() == Attribute{Gamma})
            interface.push_back(face->getIndex());

        RealFunction target(
          [&](const Geometry::Point& point) { return levelSet.getValue(point); });
        Adaptation::AnalyticVectorFunction targetGradient(
          [&](const Geometry::Point& point) {
            return projectedGradient.getValue(point);
          }, 3);

        const auto report = fitting.solve(mesh, interface, target, targetGradient);
        Alert::Info() << substageHeading("WNGIR reconstruction") << Alert::NewLine
                      << diagnosticLabel("Iterations:")
                      << Alert::Notation::Number(report.iterations) << Alert::NewLine


                      << diagnosticLabel("Skeleton normal jump RMS:")
                      << Alert::Notation::Number(report.normalJumpRMS) << Alert::NewLine
                      << diagnosticLabel("Exit reason:") << report.exitReason
                      << Alert::NewLine << diagnosticLabel("Active RMS:")
                      << Alert::Notation::Number(report.activeRMS) << Alert::NewLine
                      << diagnosticLabel("Active RMS / level-set mesh scale:")
                      << Alert::Notation::Number(
                           report.levelSetGradientScale > 0
                             ? report.activeRMS /
                                 (parameters.h * report.levelSetGradientScale)
                             : 0)
                      << Alert::NewLine
                      << diagnosticLabel("RMS stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauRms)
                      << Alert::NewLine
                      << diagnosticLabel("Scaled RMS stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauRmsH)
                      << Alert::NewLine
                      << diagnosticLabel("Active supremum:")
                      << Alert::Notation::Number(report.activeSup) << Alert::NewLine
                      << diagnosticLabel("Supremum stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauInf)
                      << Alert::NewLine
                      << diagnosticLabel("Scaled supremum stopping tolerance:")
                      << Alert::Notation::Number(report.effectiveTauInfH)
                      << Alert::NewLine
                      << diagnosticLabel("Minimum Jacobian:")
                      << Alert::Notation::Number(report.minJ) << Alert::NewLine
                      << diagnosticLabel("Maximum relative distortion:")
                      << Alert::Notation::Number(report.maxQRel) << Alert::Raise;

        KelvinBall::Mesh moved(mesh);
        moveMesh(moved, mesh, displacementTrial.getSolution());
        checkFixedGeometry(moved, outerRadius);
        checkMaterials(moved);
        const MeshDiagnostics diagnostics = getMeshDiagnostics(moved);
        return {MMG::Mesh(std::move(moved)),
          {h, h, 0, parameters.fixedBoundaryAttributes.size(), diagnostics.cells,
            diagnostics.cells}};
      }

      void printUsage(const char* executable)
      {
        Alert::Info()
          << "Usage" << Alert::NewLine << "  " << executable << " [options]"
          << Alert::NewLine << Alert::Notation("--n=<points>")
          << "              Background points per edge (default: 13)." << Alert::NewLine
          << Alert::Notation("--h=<size>")
          << "                Requested background size; alternative to --n."
          << Alert::NewLine << Alert::Notation("--outer-radius=<value>")
          << "     Chamber outer radius (default: 2)." << Alert::NewLine
          << Alert::Notation("--iterations=<count>")
          << "       Number of evaluated designs (default: 1)." << Alert::NewLine
          << Alert::Notation("--penalty=<value>")
          << "          Rotational Nitsche penalty (default: 320)." << Alert::NewLine
          << Alert::Notation("--stabilization=<value>")
          << "    P1--P1 pressure-stabilization factor (default: 0.05)." << Alert::NewLine
          << Alert::Notation("--regularization=<value>")
          << "   H1 smoothing length in multiples of h (default: 4)." << Alert::NewLine
          << Alert::Notation("--step=<value>")
          << "             Advection time step in multiples of h (default: 0.1)."
          << Alert::NewLine << Alert::Notation("--advection-quadrature=<order>")
          << " Quadrature order for the transported distance (default: 8)."
          << Alert::NewLine << Alert::Notation("--reconstruction=<method>")
          << "  Interface reconstruction: mmg or wngir (default: mmg)."
          << Alert::NewLine << Alert::Notation("--background-hmin=<value>")
          << "  WNGIR background minimum size, in h (default: 0.1)."
          << Alert::NewLine << Alert::Notation("--background-hmax=<value>")
          << "  WNGIR background maximum size, in h (default: 1)."
          << Alert::NewLine << Alert::Notation("--background-hausdorff=<value>")
          << " WNGIR background Hausdorff tolerance, in h (default: 0.05)."
          << Alert::NewLine << Alert::Notation("--background-gradation=<value>")
          << " WNGIR background gradation (default: 2)."
          << Alert::NewLine << Alert::Notation("--wngir-*=<value>")
          << "        WNGIR fitting parameters (--wngir-steps defaults to 12;"
          << Alert::NewLine
          << "                              also --trace, --j-safe, --j-ls, --j-min)."
          << Alert::NewLine << Alert::Notation("--geometry-only")
          << "             Stop after initial MMG reconstruction." << Alert::NewLine
          << Alert::Notation("--state-only")
          << "                Stop after the first Stokes evaluation." << Alert::NewLine
          << Alert::Notation("--save-mesh")
          << "                 Save KelvinBallInitial.mesh." << Alert::NewLine
          << Alert::Notation("--help") << "                      Show this message."
          << Alert::Raise;
      }

      void announce(const char* stage)
      {
        Alert::Info() << stageHeading(stage) << Alert::Raise;
      }

      using Clock = std::chrono::steady_clock;

      Real elapsedSeconds(const Clock::time_point& start) const
      {
        return std::chrono::duration<Real>(Clock::now() - start).count();
      }

      void reportStageTiming(size_t stage, Real seconds)
      {
        Alert::Info() << substageHeading("Stage " + std::to_string(stage) + " timing")
                      << Alert::NewLine << diagnosticLabel("Elapsed time:")
                      << Alert::Notation::Number(seconds) << " s" << Alert::Raise;
      }

      using LocalMesh = Geometry::Mesh<Context::Local>;
      using VelocitySpace = P1<Math::SpatialVector<Real>, LocalMesh>;
      using PressureSpace = P1<Real, LocalMesh>;

      VelocitySpace makeVelocitySpace(const LocalMesh& mesh)
      {
        return VelocitySpace(mesh, 3);
      }

      Real planeResidual(
        Attribute attribute, const Math::SpatialPoint& x, Real outerRadius)
      {
        switch (attribute)
        {
          case Outer:
            return std::abs(x(0) - outerRadius);
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

      void checkFixedGeometry(const LocalMesh& mesh, Real outerRadius)
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
            residual = std::max(residual, planeResidual(*attribute, x, outerRadius));
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

      MeshDiagnostics getMeshDiagnostics(
        const LocalMesh& mesh, bool requireCompletePartition = true)
      {
        MeshDiagnostics diagnostics;
        diagnostics.vertices = mesh.getVertexCount();
        diagnostics.cells = mesh.getCellCount();
        Real qualitySum = 0;
        Real elementSizeSum = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          if (cell->getAttribute() == Obstacle)
            ++diagnostics.obstacleCells;
          else if (cell->getAttribute() == Fluid)
            ++diagnostics.fluidCells;
          else
            throw std::runtime_error("A cell has an unexpected material attribute.");

          const auto& vertices = cell->getVertices();
          if (vertices.size() != 4)
            throw std::runtime_error("Kelvin-ball mesh quality requires tetrahedra.");
          Real squaredEdgeLengthSum = 0;
          Real edgeLengthSum = 0;
          for (size_t i = 0; i < vertices.size(); ++i)
          {
            for (size_t j = i + 1; j < vertices.size(); ++j)
            {
              const Real edgeLength = (mesh.getVertexCoordinates(vertices[i]) -
                mesh.getVertexCoordinates(vertices[j]))
                                          .norm();
              edgeLengthSum += edgeLength;
              squaredEdgeLengthSum += edgeLength * edgeLength;
            }
          }
          const Real quality =
            12 * std::pow(3 * cell->getMeasure(), Real(2) / 3) / squaredEdgeLengthSum;
          diagnostics.minimumQuality = std::min(diagnostics.minimumQuality, quality);
          diagnostics.maximumQuality = std::max(diagnostics.maximumQuality, quality);
          qualitySum += quality;
          elementSizeSum += edgeLengthSum / 6;
        }
        diagnostics.meanQuality = qualitySum / static_cast<Real>(diagnostics.cells);
        diagnostics.meanElementSize =
          elementSizeSum / static_cast<Real>(diagnostics.cells);
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          if (face->getAttribute() == Gamma)
            ++diagnostics.interfaceTriangles;
        }
        if (requireCompletePartition &&
          (diagnostics.obstacleCells == 0 || diagnostics.fluidCells == 0 ||
            diagnostics.interfaceTriangles == 0))
          throw std::runtime_error("The fitted body/fluid partition is incomplete.");
        return diagnostics;
      }

      void checkMaterials(const LocalMesh& mesh)
      {
        const auto diagnostics = getMeshDiagnostics(mesh);
        Alert::Info() << substageHeading("Material partition") << Alert::NewLine
                      << diagnosticLabel("Obstacle cells:")
                      << Alert::Notation::Number(diagnostics.obstacleCells)
                      << Alert::NewLine << diagnosticLabel("Fluid cells:")
                      << Alert::Notation::Number(diagnostics.fluidCells) << Alert::NewLine
                      << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(diagnostics.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Maximum tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.maximumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean element size:")
                      << Alert::Notation::Number(diagnostics.meanElementSize)
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
                      << Alert::Notation::Number(mesh.getVolume(Obstacle))
                      << Alert::NewLine << diagnosticLabel("Vertex RMS error:")
                      << Alert::Notation::Number(rms) << Alert::NewLine
                      << diagnosticLabel("Vertex maximum error:")
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
        const auto& rotations = KelvinBall::SewedOutput::getCubeRotations();
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
      MMGReconstruction discretizeLevelSetMMG(
        MMG::Mesh& mesh, const LevelSet& levelSet, Real h)
      {
        const size_t previousCells = mesh.getCellCount();
        const size_t requiredTriangles = protectFixedGeometry(mesh);
        const Real hmin = 0.1 * h;
        const Real hmax = 10 * h;
        const Real hausdorff = 0.1 * h * h;
        const MeshDiagnostics inputDiagnostics = getMeshDiagnostics(mesh);

        Alert::Info() << substageHeading("MMG input") << Alert::NewLine
                      << diagnosticLabel("Level-set minimum:")
                      << Alert::Notation::Number(levelSet.min()) << Alert::NewLine
                      << diagnosticLabel("Level-set maximum:")
                      << Alert::Notation::Number(levelSet.max()) << Alert::NewLine
                      << diagnosticLabel("Vertices:")
                      << Alert::Notation::Number(inputDiagnostics.vertices)
                      << Alert::NewLine << diagnosticLabel("Cells:")
                      << Alert::Notation::Number(inputDiagnostics.cells) << Alert::NewLine
                      << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(inputDiagnostics.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(inputDiagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(inputDiagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Mean element size:")
                      << Alert::Notation::Number(inputDiagnostics.meanElementSize)
                      << Alert::Raise;

        MMG::LevelSetDiscretizer discretizer;
        discretizer.split(Fluid, {Obstacle, Fluid})
          .setHMin(hmin)
          .setHMax(hmax)
          .setHausdorff(hausdorff)
          .setGradation(remeshGradation)
          .setBaseReferences(FlatSet<Attribute>{Obstacle, Fluid})
          .setBoundaryReference(Gamma)
          .setAngleDetection(false);
        discretizer.split(Obstacle, {Obstacle, Fluid});
        MMG::Mesh reconstructed = discretizer.discretize(levelSet);

        splitSelfPairedCut(reconstructed);
        const MeshDiagnostics reconstructionDiagnostics =
          getMeshDiagnostics(reconstructed, false);
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
                      << Alert::Notation::Number(reconstructionDiagnostics.cells)
                      << Alert::NewLine << diagnosticLabel("Obstacle cells:")
                      << Alert::Notation::Number(reconstructionDiagnostics.obstacleCells)
                      << Alert::NewLine << diagnosticLabel("Fluid cells:")
                      << Alert::Notation::Number(reconstructionDiagnostics.fluidCells)
                      << Alert::NewLine << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(inputDiagnostics.interfaceTriangles)
                      << " -> "
                      << Alert::Notation::Number(
                           reconstructionDiagnostics.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(inputDiagnostics.minimumQuality)
                      << " -> "
                      << Alert::Notation::Number(reconstructionDiagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(inputDiagnostics.meanQuality) << " -> "
                      << Alert::Notation::Number(reconstructionDiagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Mean element size:")
                      << Alert::Notation::Number(inputDiagnostics.meanElementSize) << " -> "
                      << Alert::Notation::Number(reconstructionDiagnostics.meanElementSize)
                      << Alert::Raise;

        const size_t requiredTrianglesBeforeOptimization =
          protectFixedGeometry(reconstructed);
        MMG::Optimizer()
          .setHMin(hmin)
          .setHMax(hmax)
          .setHausdorff(hausdorff)
          .setGradation(remeshGradation)
          .setAngleDetection(false)
          .optimize(reconstructed);
        splitSelfPairedCut(reconstructed);
        const size_t requiredTrianglesAfterOptimization =
          protectFixedGeometry(reconstructed);
        const MeshDiagnostics outputDiagnostics = getMeshDiagnostics(reconstructed, false);
        Alert::Info() << substageHeading("MMG optimization") << Alert::NewLine
                      << diagnosticLabel("Required boundary triangles:")
                      << Alert::Notation::Number(requiredTrianglesBeforeOptimization)
                      << " -> "
                      << Alert::Notation::Number(requiredTrianglesAfterOptimization)
                      << Alert::NewLine << diagnosticLabel("Cell count:")
                      << Alert::Notation::Number(reconstructionDiagnostics.cells) << " -> "
                      << Alert::Notation::Number(outputDiagnostics.cells)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(reconstructionDiagnostics.minimumQuality)
                      << " -> "
                      << Alert::Notation::Number(outputDiagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(reconstructionDiagnostics.meanQuality) << " -> "
                      << Alert::Notation::Number(outputDiagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Mean element size:")
                      << Alert::Notation::Number(reconstructionDiagnostics.meanElementSize)
                      << " -> " << Alert::Notation::Number(outputDiagnostics.meanElementSize)
                      << Alert::Raise;
        return {std::move(reconstructed),
          {hmin, hmax, hausdorff, requiredTrianglesAfterOptimization, previousCells,
            outputDiagnostics.cells}};
      }

      template <class Space, class Density, class Output>
      GradientDiagnostics identifyGradient(const Space& shapeSpace,
        const Density& density, const KelvinBall::RotatedNitscheIntegrator& coupling,
        Output& gradient, Real regularizationLength, Real nitschePenalty,
        const char* name)
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
        coupling.assembleVector(shapeSpace, system,
          regularizationLength * regularizationLength, nitschePenalty,
          FlatSet<Attribute>{Outer});
        KelvinBall::solveDirect(identification);
        const Real linearResidual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          std::max(system.getVector().norm(), Real(1));
        if (!std::isfinite(linearResidual) || linearResidual > linearResidualTolerance)
          throw std::runtime_error(
            "The gradient identification problem did not converge.");
        gradient = gradientTrial.getSolution();
        const Real jump = coupling.vectorJump(gradient);
        Alert::Info() << substageHeading(std::string(name) + " gradient")
                      << Alert::NewLine << diagnosticLabel("Linear residual:")
                      << Alert::Notation::Number(linearResidual) << Alert::NewLine
                      << diagnosticLabel("Rotated jump:") << Alert::Notation::Number(jump)
                      << Alert::Raise;
        return {linearResidual, jump};
      }

      int m_argc;
      char** m_argv;
  };
}

KelvinBall::KelvinBallOptimization::KelvinBallOptimization(int argc, char** argv)
  : m_argc(argc),
    m_argv(argv)
{}

int KelvinBall::KelvinBallOptimization::Implementation::run()
{
  const int argc = m_argc;
  char** argv = m_argv;
  Configuration configuration;
  size_t maxIterations = 1;
  bool saveMeshDiagnostic = false;
  bool geometryOnly = false;
  bool stateOnly = false;
  Real regularizationFactor = 4.0;
  Real stepFactor = 0.1;
  size_t advectionQuadratureOrder = 8;
  std::string reconstructionMethod = "mmg";
  for (int argument = 1; argument < argc; ++argument)
  {
    const std::string_view mode(argv[argument]);
    if (configuration.parse(mode))
      continue;
    if (mode.rfind("--iterations=", 0) == 0)
      maxIterations = std::stoul(std::string(mode.substr(13)));
    else if (mode.rfind("--regularization=", 0) == 0)
      regularizationFactor = std::stod(std::string(mode.substr(17)));
    else if (mode.rfind("--step=", 0) == 0)
      stepFactor = std::stod(std::string(mode.substr(7)));
    else if (mode.rfind("--advection-quadrature=", 0) == 0)
      advectionQuadratureOrder = std::stoul(std::string(mode.substr(23)));
    else if (mode.rfind("--reconstruction=", 0) == 0)
      reconstructionMethod = std::string(mode.substr(17));
    else if (mode.rfind("--wngir-", 0) == 0 || mode.rfind("--quad-order=", 0) == 0 ||
      mode == "--trace" || mode.rfind("--trace=", 0) == 0 ||
      mode.rfind("--j-safe=", 0) == 0 || mode.rfind("--j-ls=", 0) == 0 ||
      mode.rfind("--j-min=", 0) == 0)
      continue;
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
  configuration.finalize();
  const size_t points = configuration.points;
  const Real outerRadius = configuration.outerRadius;
  const Real nitschePenalty = configuration.nitschePenalty;
  const Real stabilizationFactor = configuration.stabilizationFactor;
  if (maxIterations == 0)
    throw std::runtime_error("The iteration count must be positive.");
  if (!(regularizationFactor > 0))
    throw std::runtime_error("The regularization factor must be positive.");
  if (!(stepFactor > 0))
    throw std::runtime_error("The advection step factor must be positive.");
  if (advectionQuadratureOrder == 0)
    throw std::runtime_error("The advection quadrature order must be positive.");
  if (reconstructionMethod != "mmg" && reconstructionMethod != "wngir")
    throw std::runtime_error("The reconstruction method must be mmg or wngir.");
  if (geometryOnly && stateOnly)
    throw std::runtime_error("Use either --geometry-only or --state-only, not both.");
  const Real h = configuration.getH();
  const Real hilbertLength = regularizationFactor * h;
  const Real dt = stepFactor * h;
  Alert::Info configurationInfo;
  configurationInfo << substageHeading("Configuration") << Alert::NewLine
                << diagnosticLabel("Grid points:") << Alert::Notation::Number(points)
                << Alert::NewLine << diagnosticLabel("Outer radius:")
                << Alert::Notation::Number(outerRadius) << Alert::NewLine
                << diagnosticLabel("Effective h:") << Alert::Notation::Number(h)
                << Alert::NewLine << diagnosticLabel("Evaluated designs:")
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
                << " = " << Alert::Notation::Number(stepFactor) << " h"
                << Alert::NewLine << diagnosticLabel("Advection quadrature order:")
                << Alert::Notation::Number(advectionQuadratureOrder) << Alert::NewLine
                << diagnosticLabel("Reconstruction method:") << reconstructionMethod;
  if (reconstructionMethod == "wngir")
  {
    configurationInfo << Alert::NewLine << diagnosticLabel("Background minimum size:")
                      << Alert::Notation::Number(configuration.backgroundHMin) << " h"
                      << Alert::NewLine << diagnosticLabel("Background maximum size:")
                      << Alert::Notation::Number(configuration.backgroundHMax) << " h"
                      << Alert::NewLine << diagnosticLabel("Background Hausdorff:")
                      << Alert::Notation::Number(configuration.backgroundHausdorff) << " h"
                      << Alert::NewLine << diagnosticLabel("Background gradation:")
                      << Alert::Notation::Number(configuration.backgroundGradation);
  }
  configurationInfo << Alert::Raise;
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const auto stage1Start = Clock::now();
  announce(reconstructionMethod == "wngir"
      ? "Stage 1: Preparing the background mesh and fitting the initial sphere with WNGIR."
      : "Stage 1: Discretizing the initial sphere with MMG.");
  Sphere sphere(configuration);
  SphereDiscretization initial = reconstructionMethod == "wngir"
    ? sphere.prepareWNGIRBackground()
    : sphere.discretize();
  ReconstructionDiagnostics reconstruction = initial.diagnostics;
  Optional<MMG::Mesh> wngirBackground;
  MMG::Mesh mesh;
  if (reconstructionMethod == "wngir")
  {
    wngirBackground.emplace(std::move(initial.mesh));
    P1 sphereSpace(*wngirBackground);
    GridFunction sphereLevelSet(sphereSpace);
    sphereLevelSet = RealFunction([](const Geometry::Point& point) {
      return point.getPhysicalCoordinates().norm() - Real(1);
    });
    MMG::Mesh classified =
      classifyLevelSetForWNGIR(*wngirBackground, sphereLevelSet);
    P1 classifiedSphereSpace(classified);
    GridFunction classifiedSphereLevelSet(classifiedSphereSpace);
    classifiedSphereLevelSet.getData() = sphereLevelSet.getData();
    MMGReconstruction fitted =
      fitLevelSetWNGIR(
        classified, classifiedSphereLevelSet, h, outerRadius, argc, argv);
    mesh = std::move(fitted.mesh);
  }
  else
  {
    mesh = std::move(initial.mesh);
  }
  checkFixedGeometry(mesh, outerRadius);
  checkMaterials(mesh);
  if (saveMeshDiagnostic)
    mesh.save("KelvinBallInitial.mesh", IO::FileFormat::MEDIT);
  reportSphereGeometry(mesh);
  const Real stage1Seconds = elapsedSeconds(stage1Start);
  reportStageTiming(1, stage1Seconds);
  if (geometryOnly)
    return 0;
  const Real targetVolume = mesh.getVolume(Obstacle);

  std::ofstream history("kelvin-ball.csv");
  history.precision(17);
  history << "iteration,outer_radius,h,dt,advection_quadrature_order,"
             "regularization_length,nitsche_penalty,"
             "stabilization_factor,assembly_backend,direct_solver,"
             "vertices,cells,obstacle_cells,fluid_cells,"
             "interface_triangles,mesh_quality_min,mesh_quality_mean,mesh_quality_max,"
             "mean_element_size,"
             "remesh_hmin,remesh_hmax,remesh_hausdorff,"
             "required_boundary_triangles,remesh_cells_before,remesh_cells_after,"
             "k,c,q,rho,coupling_symmetry,nitsche_jump,volume,target_volume,"
             "volume_error,chamber_bbox_x,chamber_bbox_y,chamber_bbox_z,"
             "sewn_bbox_x,sewn_bbox_y,sewn_bbox_z,actual_delta_rho,"
             "predicted_delta_rho,actual_delta_volume,predicted_delta_volume,"
             "null_space_multiplier,xi_rho_inf_norm,theta_inf_norm,d_rho_theta,"
             "d_volume_theta,required_d_volume_theta,rho_gradient_residual,"
             "rho_gradient_jump,volume_gradient_residual,"
             "volume_gradient_jump,stage_1_seconds,stage_2_seconds,"
             "stage_3_seconds,stage_4_seconds,stage_5_seconds,"
             "stage_6_seconds,stage_7_seconds,stage_8_seconds,"
             "stage_9_seconds\n";
  IO::XDMF xdmf("KelvinBall");
  auto chamber = xdmf.grid("Chamber");
  chamber.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
  auto fluidState = xdmf.grid("Fluid");
  IO::XDMF sewedXdmf("KelvinBallSewed");
  auto sewedDesignOutput = sewedXdmf.grid("Design");
  auto sewedFluidOutput = sewedXdmf.grid("Fluid");
  const std::string reconstructionName = reconstructionMethod == "wngir"
    ? "KelvinBallWNGIR"
    : "KelvinBallMMG";
  IO::XDMF reconstructionXdmf(reconstructionName);
  auto reconstructionOutput = reconstructionXdmf.grid("Reconstructed");
  Optional<Real> previousRho;
  Optional<Real> predictedRhoChange;
  Optional<Real> previousVolume;
  Optional<Real> predictedVolumeChange;
  Optional<MMG::Mesh> nextMesh;

  for (size_t iteration = 0; iteration < maxIterations; ++iteration)
  {
    std::array<Real, 9> stageSeconds;
    stageSeconds.fill(nan);
    if (iteration == 0)
      stageSeconds[0] = stage1Seconds;
    if (nextMesh)
    {
      mesh = std::move(*nextMesh);
      nextMesh.reset();
    }
    Alert::Info() << " ------------------------------------------------------------"
                  << Alert::NewLine << " Iteration "
                  << Alert::Notation::Number(iteration + 1) << " of "
                  << Alert::Notation::Number(maxIterations) << Alert::NewLine
                  << " ------------------------------------------------------------"
                  << Alert::Raise;
    const auto stage2Start = Clock::now();
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
    KelvinBall::RotatedNitscheIntegrator fluidCoupling(
      fluid, MasterCuts, rotatedTracePhysicalTolerance, rotatedTraceReferenceTolerance);
    GridFunction uT0(Vh), uT1(Vh), uT2(Vh), uR0(Vh), uR1(Vh), uR2(Vh);
    GridFunction pT0(Qh), pT1(Qh), pT2(Qh), pR0(Qh), pR1(Qh), pR2(Qh);
    stageSeconds[1] = elapsedSeconds(stage2Start);
    reportStageTiming(2, stageSeconds[1]);

    const auto stage3Start = Clock::now();
    announce("Stage 3: Solving the translational and rotational Stokes states.");
    const KelvinBall::Parameters metricParameters{h, nitschePenalty, stabilizationFactor};
    const KelvinBall::Metrics resistanceMetrics(metricParameters);
    const KelvinBall::Values metrics = resistanceMetrics.evaluateChamber(
      Vh, Qh, fluidCoupling, uT0, uT1, uT2, uR0, uR1, uR2, pT0, pT1, pT2, pR0, pR1, pR2);
    stageSeconds[2] = elapsedSeconds(stage3Start);
    reportStageTiming(3, stageSeconds[2]);
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
    const Real actualRhoChange = previousRho ? rho - *previousRho : nan;
    const Real incomingPredictedRhoChange =
      predictedRhoChange ? *predictedRhoChange : nan;
    const Real actualVolumeChange = previousVolume ? volume - *previousVolume : nan;
    const Real incomingPredictedVolumeChange =
      predictedVolumeChange ? *predictedVolumeChange : nan;
    const auto writeHistory = [&](Real nullSpaceMultiplier, Real xiRhoInfinityNorm,
                                Real thetaInfinityNorm, Real dRhoTheta, Real dVolumeTheta,
                                Real requiredDVolumeTheta,
                                const GradientDiagnostics& rhoGradientDiagnostics,
                                const GradientDiagnostics& volumeGradientDiagnostics) {
      history << iteration << ',' << outerRadius << ',' << h << ',' << dt << ','
              << advectionQuadratureOrder << ',' << hilbertLength << ','
              << nitschePenalty << ',' << stabilizationFactor << ',' << assemblyBackend
              << ',' << KelvinBall::DirectSolverName << ','
              << meshDiagnostics.vertices << ',' << meshDiagnostics.cells << ','
              << meshDiagnostics.obstacleCells << ',' << meshDiagnostics.fluidCells << ','
              << meshDiagnostics.interfaceTriangles << ','
              << meshDiagnostics.minimumQuality << ',' << meshDiagnostics.meanQuality
              << ',' << meshDiagnostics.maximumQuality << ','
              << meshDiagnostics.meanElementSize << ','
              << reconstruction.minimumSize << ',' << reconstruction.maximumSize << ','
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
              << incomingPredictedVolumeChange << ',' << nullSpaceMultiplier << ','
              << xiRhoInfinityNorm << ',' << thetaInfinityNorm << ',' << dRhoTheta << ','
              << dVolumeTheta << ',' << requiredDVolumeTheta << ','
              << rhoGradientDiagnostics.residual << ',' << rhoGradientDiagnostics.jump
              << ',' << volumeGradientDiagnostics.residual << ','
              << volumeGradientDiagnostics.jump;
      for (const Real seconds : stageSeconds)
        history << ',' << seconds;
      history << '\n';
      history.flush();
    };
    const auto stage4Start = Clock::now();
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
    stageSeconds[3] = elapsedSeconds(stage4Start);
    reportStageTiming(4, stageSeconds[3]);
    if (stateOnly)
    {
      const GradientDiagnostics unavailable{nan, nan};
      writeHistory(nan, nan, nan, nan, nan, nan, unavailable, unavailable);
      return 0;
    }

    const auto stage5Start = Clock::now();
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
      (Dot(strainUT0, strainUT0) + Dot(strainUT1, strainUT1) + Dot(strainUT2, strainUT2));
    const auto cDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUT0, strainUR0) + Dot(strainUT1, strainUR1) + Dot(strainUT2, strainUR2));
    const auto qDensity = chamberMultiplicity * 2.0 * mu / 3.0 *
      (Dot(strainUR0, strainUR0) + Dot(strainUR1, strainUR1) + Dot(strainUR2, strainUR2));
    const auto rhoDensity =
      signC / std::sqrt(k * q) * cDensity - 0.5 * rho * (kDensity / k + qDensity / q);

    P1 levelSetSpace(mesh);
    P1 shapeSpace(mesh, 3);
    KelvinBall::RotatedNitscheIntegrator shapeCoupling(mesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus},
      rotatedTracePhysicalTolerance, rotatedTraceReferenceTolerance);
    GridFunction rhoGradient(shapeSpace), volumeGradient(shapeSpace);
    stageSeconds[4] = elapsedSeconds(stage5Start);
    reportStageTiming(5, stageSeconds[4]);
    const auto stage6Start = Clock::now();
    announce("Stage 6: Regularizing and constraining the shape direction.");
    const GradientDiagnostics rhoGradientDiagnostics = identifyGradient(shapeSpace,
      -rhoDensity, shapeCoupling, rhoGradient, hilbertLength, nitschePenalty, "Rho");
    const GradientDiagnostics volumeGradientDiagnostics =
      identifyGradient(shapeSpace, RealFunction{-1}, shapeCoupling, volumeGradient,
        hilbertLength, nitschePenalty, "Volume");
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
                  << Alert::Notation::Number(thetaInfinityNorm) << Alert::NewLine
                  << diagnosticLabel("D rho [theta]:")
                  << Alert::Notation::Number(dRhoTheta) << Alert::NewLine
                  << diagnosticLabel("D volume [theta]:")
                  << Alert::Notation::Number(dVolumeTheta) << Alert::NewLine
                  << diagnosticLabel("Required D volume [theta]:")
                  << Alert::Notation::Number(requiredDVolumeTheta) << Alert::Raise;
    previousRho = rho;
    predictedRhoChange = dt * dRhoTheta;
    previousVolume = volume;
    predictedVolumeChange = dt * dVolumeTheta;
    GridFunction eikonalDistance(levelSetSpace);
    Distance::Eikonal(eikonalDistance)
      .setInterior(Obstacle)
      .setInterface(Gamma)
      .solve()
      .sign();
    TrialFunction periodicDistance(levelSetSpace);
    TestFunction distanceTest(levelSetSpace);
    Problem distanceProjection(periodicDistance, distanceTest);
    auto eikonalDistanceLoad = Integral(eikonalDistance, distanceTest);
    eikonalDistanceLoad.setOrder(2);
    distanceProjection =
      Integral(periodicDistance, distanceTest) - eikonalDistanceLoad;
    distanceProjection.assemble();
    shapeCoupling.assembleScalarTracePenalty(
      levelSetSpace, distanceProjection.getLinearSystem(), nitschePenalty);
    Solver::CG(distanceProjection).solve();
    GridFunction distance(levelSetSpace);
    distance = periodicDistance.getSolution();
    const auto& distanceSystem = distanceProjection.getLinearSystem();
    const Real distanceProjectionResidual =
      (distanceSystem.getOperator() * distanceSystem.getSolution() -
        distanceSystem.getVector())
          .norm() /
      std::max(distanceSystem.getVector().norm(), Real(1));
    const Real distanceProjectionCorrection =
      (distance.getData() - eikonalDistance.getData()).lpNorm<Eigen::Infinity>();
    Alert::Info() << substageHeading("Distance trace projection") << Alert::NewLine
                  << diagnosticLabel("Linear residual:")
                  << Alert::Notation::Number(distanceProjectionResidual)
                  << Alert::NewLine << diagnosticLabel("Eikonal rotated jump:")
                  << Alert::Notation::Number(shapeCoupling.scalarJump(eikonalDistance))
                  << Alert::NewLine << diagnosticLabel("Projected rotated jump:")
                  << Alert::Notation::Number(shapeCoupling.scalarJump(distance))
                  << Alert::NewLine << diagnosticLabel("Infinity correction:")
                  << Alert::Notation::Number(distanceProjectionCorrection)
                  << Alert::Raise;

    stageSeconds[5] = elapsedSeconds(stage6Start);
    reportStageTiming(6, stageSeconds[5]);
    const auto stage7Start = Clock::now();
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

    KelvinBall::SewedOutput sewedDesign(mesh, FlatSet<Attribute>{Gamma, Outer});
    P1 sewedDesignScalar(sewedDesign.getMesh());
    P1 sewedDesignVector(sewedDesign.getMesh(), 3);
    GridFunction sewedDistance(sewedDesignScalar);
    GridFunction sewedVelocity(sewedDesignVector);
    sewedDesign.setScalar(sewedDistance, distance);
    sewedDesign.setVector(sewedVelocity, theta);
    sewedDesignOutput.clear();
    sewedDesignOutput.setMesh(sewedDesign.getMesh(), IO::XDMF::MeshPolicy::Transient);
    sewedDesignOutput.add("Distance", sewedDistance, IO::XDMF::Center::Node);
    sewedDesignOutput.add("Theta", sewedVelocity, IO::XDMF::Center::Node);

    KelvinBall::SewedOutput sewedFluid(fluid, FlatSet<Attribute>{Gamma, Outer});
    VelocitySpace sewedVelocitySpace = makeVelocitySpace(sewedFluid.getMesh());
    PressureSpace sewedPressureSpace(sewedFluid.getMesh());
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
    sewedFluid.setVectorLoad(sewedUT0, translations, 0);
    sewedFluid.setVectorLoad(sewedUT1, translations, 1);
    sewedFluid.setVectorLoad(sewedUT2, translations, 2);
    sewedFluid.setVectorLoad(sewedUR0, rotations, 0);
    sewedFluid.setVectorLoad(sewedUR1, rotations, 1);
    sewedFluid.setVectorLoad(sewedUR2, rotations, 2);
    sewedFluid.setScalarLoad(sewedPT0, translationPressures, 0);
    sewedFluid.setScalarLoad(sewedPT1, translationPressures, 1);
    sewedFluid.setScalarLoad(sewedPT2, translationPressures, 2);
    sewedFluid.setScalarLoad(sewedPR0, rotationPressures, 0);
    sewedFluid.setScalarLoad(sewedPR1, rotationPressures, 1);
    sewedFluid.setScalarLoad(sewedPR2, rotationPressures, 2);
    sewedFluidOutput.clear();
    sewedFluidOutput.setMesh(sewedFluid.getMesh(), IO::XDMF::MeshPolicy::Transient);
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

    if (iteration + 1 == maxIterations)
    {
      xdmf.write(static_cast<Real>(iteration)).flush();
      sewedXdmf.write(static_cast<Real>(iteration)).flush();
      stageSeconds[6] = elapsedSeconds(stage7Start);
      reportStageTiming(7, stageSeconds[6]);
      writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
        dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
        volumeGradientDiagnostics);
      continue;
    }

    stageSeconds[6] = elapsedSeconds(stage7Start);
    reportStageTiming(7, stageSeconds[6]);
    const auto stage8Start = Clock::now();
    announce("Stage 8: Advecting the level set.");
    MMG::Mesh& advectionMesh = wngirBackground ? *wngirBackground : mesh;
    auto& advectionConnectivity = advectionMesh.getConnectivity();
    advectionConnectivity.discover(3, 2);
    advectionConnectivity.discover(3, 1);
    advectionConnectivity.restrict(1, 0);
    advectionConnectivity.restrict(2, 0);
    advectionConnectivity.restrict(2, 3);
    advectionConnectivity.discover(0, 0);
    P1 advectionLevelSetSpace(advectionMesh);
    P1 advectionShapeSpace(advectionMesh, 3);
    GridFunction advectionDistance(advectionLevelSetSpace);
    GridFunction advectionDirection(advectionShapeSpace);
    if (wngirBackground)
    {
      // The fitted mesh is the background with its vertices moved: the same
      // numbering at different positions. Copying nodal values would carry
      // each one from x + u(x) back to x and undo the fit near the interface,
      // so both fields are evaluated at the positions of the background.
      const Location::AABB<MMG::Mesh> fittedLocator(mesh);
      std::size_t unlocated = 0;
      advectionDistance = RealFunction([&](const Geometry::Point& point) {
        const auto located = fittedLocator.locate(3, point.getPhysicalCoordinates());
        if (!located)
        {
          ++unlocated;
          return Real(0);
        }
        return distance.getValue(*located);
      });
      advectionDirection = VectorFunction(static_cast<size_t>(3),
        [&](const Geometry::Point& point) {
          Math::SpatialVector<Real> value(3);
          value.setZero();
          const auto located =
            fittedLocator.locate(3, point.getPhysicalCoordinates());
          if (!located)
          {
            ++unlocated;
            return value;
          }
          const auto sample = theta.getValue(*located);
          for (Eigen::Index component = 0; component < 3; ++component)
            value(component) = sample(component);
          return value;
        });
      if (unlocated > 0)
        throw std::runtime_error(
          "Background points fell outside the fitted mesh during the transfer.");
      Alert::Info() << substageHeading("Fitted-to-background transfer")
                    << Alert::NewLine << diagnosticLabel("Nodal-copy error (distance):")
                    << Alert::Notation::Number(
                         (advectionDistance.getData() - distance.getData())
                           .lpNorm<Eigen::Infinity>())
                    << Alert::NewLine << diagnosticLabel("Nodal-copy error (direction):")
                    << Alert::Notation::Number(
                         (advectionDirection.getData() - theta.getData())
                           .lpNorm<Eigen::Infinity>())
                    << Alert::Raise;
    }
    else
    {
      advectionDistance.getData() = distance.getData();
      advectionDirection.getData() = theta.getData();
    }
    KelvinBall::RotatedNitscheIntegrator advectionCoupling(advectionMesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus},
      rotatedTracePhysicalTolerance, rotatedTraceReferenceTolerance);
    TrialFunction advected(advectionLevelSetSpace);
    TestFunction test(advectionLevelSetSpace);
    KelvinBall::RotatedCharacteristicContinuation rotationalContinuation(
      -dt, advectionMesh, advectionCoupling.getLocator(), RotationPairs);
    Problem transport(advected, test);
    auto transportedDistance =
      Integral(Flow(-dt, advectionDistance, advectionDirection, Math::RungeKutta::RK4{},
                 rotationalContinuation),
        test);
    transportedDistance.setOrder(advectionQuadratureOrder);
    transport = Integral(advected, test) - transportedDistance;
    transport.assemble();
    advectionCoupling.assembleScalarTracePenalty(advectionLevelSetSpace,
      transport.getLinearSystem(), nitschePenalty, advectionDistance);
    Solver::CG(transport).solve();
    const auto& advectedDistance = advected.getSolution();
    const auto& transportSystem = transport.getLinearSystem();
    const Real transportResidual =
      (transportSystem.getOperator() * transportSystem.getSolution() -
        transportSystem.getVector())
          .norm() /
      std::max(transportSystem.getVector().norm(), Real(1));
    const Real advectionIncrement =
      (advectedDistance.getData() - advectionDistance.getData()).lpNorm<Eigen::Infinity>();
    const Real distanceJump = advectionCoupling.scalarJump(advectionDistance);
    const Real advectedJump = advectionCoupling.scalarJump(advectedDistance);
    Alert::Info() << substageHeading("Advected distance") << Alert::NewLine
                  << diagnosticLabel("Advection step:")
                  << Alert::Notation::Number(dt) << Alert::NewLine
                  << diagnosticLabel("Linear residual:")
                  << Alert::Notation::Number(transportResidual) << Alert::NewLine
                  << diagnosticLabel("Minimum:")
                  << Alert::Notation::Number(advectedDistance.min()) << Alert::NewLine
                  << diagnosticLabel("Maximum:")
                  << Alert::Notation::Number(advectedDistance.max()) << Alert::NewLine
                  << diagnosticLabel("Increment infinity norm:")
                  << Alert::Notation::Number(advectionIncrement) << Alert::NewLine
                  << diagnosticLabel("Distance rotated jump:")
                  << Alert::Notation::Number(distanceJump) << Alert::NewLine
                  << diagnosticLabel("Advected rotated jump:")
                  << Alert::Notation::Number(advectedJump) << Alert::Raise;
    GridFunction advectedOutput(levelSetSpace);
    advectedOutput.getData() = advectedDistance.getData();
    chamber.add("Advected", advectedOutput, IO::XDMF::Center::Node);
    GridFunction sewedAdvected(sewedDesignScalar);
    sewedDesign.setScalar(sewedAdvected, advectedOutput);
    sewedDesignOutput.add("Advected", sewedAdvected, IO::XDMF::Center::Node);

    xdmf.write(static_cast<Real>(iteration)).flush();
    sewedXdmf.write(static_cast<Real>(iteration)).flush();
    stageSeconds[7] = elapsedSeconds(stage8Start);
    reportStageTiming(8, stageSeconds[7]);

    const auto stage9Start = Clock::now();
    announce(reconstructionMethod == "wngir"
        ? "Stage 9: Fitting the advected interface with WNGIR."
        : "Stage 9: Reconstructing the advected interface with MMG.");
    MMGReconstruction result = [&]() {
      if (reconstructionMethod == "wngir")
      {
        MMG::Mesh classified =
          classifyLevelSetForWNGIR(*wngirBackground, advectedDistance);
        P1 classifiedLevelSetSpace(classified);
        GridFunction classifiedLevelSet(classifiedLevelSetSpace);
        classifiedLevelSet.getData() = advectedDistance.getData();
        return fitLevelSetWNGIR(
          classified, classifiedLevelSet, h, outerRadius, argc, argv);
      }
      return discretizeLevelSetMMG(mesh, advectedDistance, h);
    }();
    reconstruction = result.diagnostics;
    reconstructionOutput.clear();
    reconstructionOutput.setMesh(result.mesh, IO::XDMF::MeshPolicy::Transient);
    reconstructionXdmf.write(static_cast<Real>(iteration + 1)).flush();
    checkFixedGeometry(result.mesh, outerRadius);
    checkMaterials(result.mesh);
    nextMesh.emplace(std::move(result.mesh));
    stageSeconds[8] = elapsedSeconds(stage9Start);
    reportStageTiming(9, stageSeconds[8]);
    writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
      dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
      volumeGradientDiagnostics);
  }

  Alert::Success()
    << "Wrote KelvinBall.xdmf, KelvinBallSewed.xdmf, " << reconstructionName
    << ".xdmf, and "
       "kelvin-ball.csv"
    << Alert::Raise;
  return 0;
}

int KelvinBall::KelvinBallOptimization::run()
{
  return Implementation(m_argc, m_argv).run();
}
