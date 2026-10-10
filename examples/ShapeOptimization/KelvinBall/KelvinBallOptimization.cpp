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
#include <sstream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Adaptation.h>
#include <Rodin/Alert/Info.h>
#include <Rodin/Location.h>
#include <Rodin/Alert/Notation.h>
#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Alert/Warning.h>
#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MMG.h>
#include <Rodin/Variational.h>

#include "Configuration.h"
#include "Metrics.h"
#include "RotatedCharacteristicContinuation.h"
#include "TransportProjection.h"
#include "RotatedNitscheIntegrator.h"
#include "SewedOutput.h"
#include "Thickness.h"
#include "HilbertIdentification.h"
#include "Sphere.h"
#include "../../Adaptation/SWIFT/Options.h"

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
      struct Options
      {
          Options()
            : motionForce(3)
          {
            motionForce.setZero();
            motionForce(0) = 1;
          }

          Configuration configuration;
          size_t maxIterations = 1;
          bool saveMeshDiagnostic = false;
          bool geometryOnly = false;
          bool stateOnly = false;
          bool thicknessDiagnostics = false;
          Real regularizationFactor = 4;
          Real normalRegularizationFactor = 1;
          Real stepFactor = 0.1;
          Real levelSetPenalty = 1;
          Real minimumThickness = 0;
          Math::SpatialVector<Real> motionForce;
          size_t advectionQuadratureOrder = 8;
          std::string reconstructionMethod = "mmg";

          bool parse(int argc, char** argv)
          {
            for (int argument = 1; argument < argc; ++argument)
            {
              const std::string_view mode(argv[argument]);
              if (configuration.parse(mode))
                continue;
              if (mode.rfind("--iterations=", 0) == 0)
                maxIterations = std::stoul(std::string(mode.substr(13)));
              else if (mode.rfind("--regularization=", 0) == 0)
                regularizationFactor = std::stod(std::string(mode.substr(17)));
              else if (mode.rfind("--normal-regularization=", 0) == 0)
                normalRegularizationFactor = std::stod(std::string(mode.substr(24)));
              else if (mode.rfind("--step=", 0) == 0)
                stepFactor = std::stod(std::string(mode.substr(7)));
              else if (mode.rfind("--level-set-penalty=", 0) == 0)
                levelSetPenalty = std::stod(std::string(mode.substr(20)));
              else if (mode.rfind("--thickness-min=", 0) == 0)
                minimumThickness = std::stod(std::string(mode.substr(16)));
              else if (mode.rfind("--motion-force=", 0) == 0)
              {
                std::istringstream components{std::string(mode.substr(15))};
                std::string component;
                for (Eigen::Index i = 0; i < 3; ++i)
                {
                  if (!std::getline(components, component, ','))
                    throw std::runtime_error(
                      "--motion-force needs three components fx,fy,fz.");
                  motionForce(i) = std::stod(component);
                }
              }
              else if (mode.rfind("--advection-quadrature=", 0) == 0)
                advectionQuadratureOrder = std::stoul(std::string(mode.substr(23)));
              else if (mode.rfind("--reconstruction=", 0) == 0)
                reconstructionMethod = std::string(mode.substr(17));
              else if (mode.rfind("--swift-", 0) == 0)
                continue;
              else if (mode == "--save-mesh")
                saveMeshDiagnostic = true;
              else if (mode == "--geometry-only")
                geometryOnly = true;
              else if (mode == "--thickness-diagnostics")
                thicknessDiagnostics = true;
              else if (mode == "--state-only")
                stateOnly = true;
              else if (mode == "--help")
              {
                printUsage(argv[0]);
                return false;
              }
              else
                throw std::runtime_error(
                  "Unknown KelvinBall option: " + std::string(mode));
            }
            configuration.finalize();
            if (maxIterations == 0)
              throw std::runtime_error("The iteration count must be positive.");
            if (!(regularizationFactor > 0))
              throw std::runtime_error("The regularization factor must be positive.");
            if (!(normalRegularizationFactor > 0))
              throw std::runtime_error(
                "The normal regularization factor must be positive.");
            if (!(stepFactor > 0))
              throw std::runtime_error("The advection step factor must be positive.");
            if (!(levelSetPenalty > 0))
              throw std::runtime_error("The level-set trace penalty must be positive.");
            if (!std::isfinite(minimumThickness) || minimumThickness < 0)
              throw std::runtime_error("--thickness-min must be nonnegative.");
            if (!std::isfinite(motionForce.norm()) || !(motionForce.norm() > 0))
              throw std::runtime_error("The motion force must be finite and nonzero.");
            if (advectionQuadratureOrder == 0)
              throw std::runtime_error(
                "The advection quadrature order must be positive.");
            if (reconstructionMethod != "mmg" && reconstructionMethod != "swift")
              throw std::runtime_error("The reconstruction method must be mmg or swift.");
            if (configuration.mmgSnap > 0 && reconstructionMethod != "mmg")
              throw std::runtime_error(
                "--mmg-snap applies only to --reconstruction=mmg.");
            if (geometryOnly && stateOnly)
              throw std::runtime_error(
                "Use either --geometry-only or --state-only, not both.");
            return true;
          }
      };

      Options m_options;

      static constexpr Real mu = KelvinBall::Mu;
      static constexpr Real remeshGradation = 2.0;
      static constexpr Real rotatedTracePhysicalTolerance = 0.01;
      static constexpr Real rotatedTraceReferenceTolerance = 0.25;
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

      /// Formats a vector as (x, y, z), with signed zeros printed as 0.
      std::string vectorText(const Math::SpatialVector<Real>& v)
      {
        std::ostringstream text;
        text << '(';
        for (Eigen::Index i = 0; i < v.size(); ++i)
          text << (i ? ", " : "") << v(i) + Real(0);
        text << ')';
        return text.str();
      }

      struct InterfaceSeamJump
      {
          Real maximum = 0;
          size_t samples = 0;
          size_t unmatched = 0;
      };

      struct InterfaceComponents
      {
          Real normalRMS = 0;
          Real tangentialRMS = 0;
          Real tangentialFraction = 0;
      };

      template <class Field>
      InterfaceComponents interfaceComponents(const Mesh& mesh, const Field& field) const
      {
        static constexpr std::array<std::array<Real, 3>, 3> quadrature{
          {{Real(2) / 3, Real(1) / 6, Real(1) / 6},
            {Real(1) / 6, Real(2) / 3, Real(1) / 6},
            {Real(1) / 6, Real(1) / 6, Real(2) / 3}}};
        auto fluidNormal = FaceNormal(mesh);
        fluidNormal.traceOf(Fluid);
        Real areaTotal = 0, normalSquared = 0, tangentialSquared = 0;
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          if (face->getAttribute() != Gamma)
            continue;
          const auto& vertices = face->getVertices();
          const auto a = mesh.getVertexCoordinates(vertices[0]);
          const auto b = mesh.getVertexCoordinates(vertices[1]);
          const auto c = mesh.getVertexCoordinates(vertices[2]);
          const Real area = (b - a).cross(c - a).norm() / Real(2);
          areaTotal += area;
          for (const auto& barycentric : quadrature)
          {
            const Geometry::Point point(
              *face, barycentric[0] * a + barycentric[1] * b + barycentric[2] * c);
            const auto velocity = field.getValue(point);
            const auto normal = -fluidNormal.getValue(point);
            const Real component = velocity.dot(normal);
            normalSquared += area / Real(3) * component * component;
            tangentialSquared += area / Real(3) *
              std::max(Real(0), velocity.squaredNorm() - component * component);
          }
        }
        InterfaceComponents result;
        if (areaTotal > 0)
        {
          result.normalRMS = std::sqrt(normalSquared / areaTotal);
          result.tangentialRMS = std::sqrt(tangentialSquared / areaTotal);
        }
        if (normalSquared + tangentialSquared > 0)
          result.tangentialFraction =
            std::sqrt(tangentialSquared / (normalSquared + tangentialSquared));
        return result;
      }

      /// Compare rotated traces only where the design interface meets a cut.
      template <class Field>
      InterfaceSeamJump interfaceSeamJump(const Mesh& mesh, const Field& field,
        const KelvinBall::RotatedNitscheIntegrator::Locator& locator) const
      {
        std::vector<bool> interfaceVertices(mesh.getVertexCount(), false);
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          if (face->getAttribute() == Gamma)
            for (const Index vertex : face->getVertices())
              interfaceVertices[vertex] = true;

        InterfaceSeamJump jump;
        for (const auto& pair : RotationPairs)
        {
          std::set<Index> sampled;
          for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
          {
            if (face->getAttribute() != pair.slave)
              continue;
            for (const Index vertex : face->getVertices())
            {
              if (!interfaceVertices[vertex] || !sampled.insert(vertex).second)
                continue;
              const auto mapped = locator.locate(
                pair.master, pair.rotation * mesh.getVertexCoordinates(vertex));
              if (!mapped)
              {
                ++jump.unmatched;
                continue;
              }
              const auto dofs = field.getFiniteElementSpace().getDOFs(0, vertex);
              Math::SpatialVector<Real> source(3);
              for (size_t component = 0; component < 3; ++component)
                source(component) = field.getData()(dofs(component));
              jump.maximum = std::max(
                jump.maximum, (field.getValue(*mapped) - pair.rotation * source).norm());
              ++jump.samples;
            }
          }
        }
        return jump;
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

      using GradientDiagnostics = HilbertIdentification::Diagnostics;

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

      using MMGReconstruction = SphereDiscretization;

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
      MMG::Mesh classifyLevelSetForSWIFT(
        const MMG::Mesh& background, const LevelSet& levelSet)
      {
        MMG::Mesh classified(background);
        const auto& space = levelSet.getFiniteElementSpace();
        const auto& coefficients = levelSet.getData();
        std::vector<Real> moments(classified.getCellCount());
        std::vector<Real> sizes(classified.getCellCount());
        constexpr Real PhaseWidth = Real(1.25);
        for (auto cell = classified.getCell(); cell; ++cell)
        {
          const Index index = cell->getIndex();
          Real meanLevelSet = 0;
          for (const Index vertex : cell->getVertices())
          {
            const auto& dofs = space.getDOFs(0, vertex);
            const Real value = coefficients(dofs[0]);
            meanLevelSet += value;
          }
          meanLevelSet /= static_cast<Real>(cell->getVertices().size());
          sizes[index] = cellSize(*cell);
          moments[index] = std::tanh(meanLevelSet / (PhaseWidth * sizes[index]));
        }

        classified.getConnectivity().compute(2, 3);
        for (auto face = classified.getFace(); face; ++face)
        {
          const auto& incident =
            classified.getConnectivity().getIncidence({2, 3}, face->getIndex());
          if (incident.size() != 2)
            continue;
          classified.setAttribute({2, face->getIndex()}, {});
        }
        const auto options = getReconstructionOptions(m_argc, m_argv);
        MinSTCut classifier(static_cast<const KelvinBall::Mesh&>(classified));
        decltype(classifier)::Parameters parameters;
        parameters.fidelity = options.classification.fidelity;
        parameters.smoothing = [&](const Polytope& face) {
          const auto& cells =
            classified.getConnectivity().getIncidence({2, 3}, face.getIndex());
          return options.classification.smoothing *
            std::min(sizes[cells[0]], sizes[cells[1]]);
        };
        classifier.setParameters(parameters);
        const auto partition = classifier.classify(
          [&](const Polytope& cell) { return moments[cell.getIndex()]; });
        Alert::Info() << substageHeading("MinSTCut classification") << Alert::NewLine
                      << diagnosticLabel("Fidelity:")
                      << Alert::Notation::Number(options.classification.fidelity)
                      << Alert::NewLine << diagnosticLabel("Smoothing / local cell size:")
                      << Alert::Notation::Number(options.classification.smoothing)
                      << Alert::NewLine << diagnosticLabel("Solid cells:")
                      << Alert::Notation::Number(partition.inside.size())
                      << Alert::NewLine << diagnosticLabel("Fluid cells:")
                      << Alert::Notation::Number(partition.outside.size())
                      << Alert::NewLine << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(partition.cut.size()) << Alert::Raise;
        for (const Index cell : partition.inside)
          classified.setAttribute({3, cell}, Obstacle);
        for (const Index cell : partition.outside)
          classified.setAttribute({3, cell}, Fluid);
        for (const Index face : partition.cut)
          classified.setAttribute({2, face}, Gamma);
        return classified;
      }

      Rodin::Examples::ReconstructionOptions getReconstructionOptions(
        int argc, char** argv) const
      {
        // Retain KelvinBall's existing dimensionless perimeter weight.
        std::vector<std::string> arguments{argv[0], "--classification-smoothing=0.04"};
        for (int i = 1; i < argc; ++i)
        {
          const std::string argument(argv[i]);
          if (!argument.starts_with("--swift-"))
            continue;
          const std::string option = "--" + argument.substr(8);
          if (!option.starts_with("--classification-") &&
            !option.starts_with("--model-") && !option.starts_with("--globalization-") &&
            !option.starts_with("--linear-") && !option.starts_with("--convergence-") &&
            !option.starts_with("--quadrature-") && !option.starts_with("--sampling-") &&
            option != "--trace" && !option.starts_with("--trace=") &&
            !option.starts_with("--trace-quality-witness"))
            throw std::runtime_error("Unknown SWIFT fitting option: " + argument);
          arguments.push_back(option);
        }
        std::vector<char*> pointers;
        for (auto& argument : arguments)
          pointers.push_back(argument.data());
        return Rodin::Examples::ReconstructionOptions(
          static_cast<int>(pointers.size()), pointers.data());
      }

      Adaptation::SWIFT::Parameters getFittingParameters(
        Real referenceSpacing, int argc, char** argv) const
      {
        auto parameters = getReconstructionOptions(argc, argv).parameters;
        FlatSet<std::string> specified;
        for (int i = 1; i < argc; ++i)
        {
          const std::string argument(argv[i]);
          if (argument.starts_with("--swift-"))
            specified.insert(argument.substr(8, argument.find('=') - 8));
        }
        parameters.model.h = referenceSpacing;
        parameters.interfaceAttribute = Gamma;
        if (!specified.contains("convergence-tolerance-geometric"))
          parameters.convergence.tolerance.geometric =
            Real(0.1) * referenceSpacing * referenceSpacing;
        if (!specified.contains("convergence-tolerance-step"))
          parameters.convergence.tolerance.step =
            Real(1e-3) * referenceSpacing * referenceSpacing;
        if (!specified.contains("convergence-tolerance-step-over-h"))
          parameters.convergence.tolerance.stepOverH = Real(1e-3) * referenceSpacing;
        if (!specified.contains("trace"))
          parameters.trace = true;
        return parameters;
      }

      template <class LevelSet>
      MMGReconstruction fitLevelSetSWIFT(const KelvinBall::Mesh& mesh,
        const LevelSet& levelSet, Real backgroundH, Real referenceSpacing,
        Real outerRadius, int argc, char** argv)
      {
        size_t interfaceCount = 0;
        Real interfaceSizeSum = 0;
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          if (face->getAttribute() == Attribute{Gamma})
          {
            ++interfaceCount;
            interfaceSizeSum +=
              std::sqrt(Real(4) * face->getMeasure() / std::sqrt(Real(3)));
          }
        }
        const Real meanInterfaceSize = interfaceCount == 0
          ? backgroundH
          : interfaceSizeSum / static_cast<Real>(interfaceCount);

        // Classification changes labels only. Copy the P1 coefficients onto
        // that mesh so the locator used by SWIFT resolves both the target and
        // its exact element-wise gradient on the same geometry.
        P1<Real, KelvinBall::Mesh> targetSpace(mesh);
        GridFunction targetLevelSet(targetSpace);
        targetLevelSet.getData() = levelSet.getData();
        auto targetGradient = Grad(targetLevelSet);
        targetGradient.traceOf(Fluid);

        P1<Math::SpatialVector<Real>, KelvinBall::Mesh> displacementSpace(mesh, 3);
        TrialFunction displacementTrial(displacementSpace);
        TestFunction displacementTest(displacementSpace);
        auto parameters = getFittingParameters(referenceSpacing, argc, argv);
        // The outer sphere is curved, so its nodes are pinned: sliding in a
        // facet plane would walk them off the sphere. The cut planes bound the
        // wedge but are not walls, and the rim of the interface lies entirely
        // on them, so they slide within themselves instead of being frozen.
        parameters.fixedBoundaryAttributes = {Outer};
        parameters.slipBoundaryAttributes = {
          SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
        Adaptation::SWIFT::Problem fitting(displacementTrial, displacementTest);
        fitting.setParameters(parameters);

        RealFunction target(
          [&](const Geometry::Point& point) { return targetLevelSet.getValue(point); });
        const auto report = fitting.solve(target, targetGradient);
        Alert::Info()
          << substageHeading("SWIFT reconstruction") << Alert::NewLine
          << diagnosticLabel("Background mean edge h:")
          << Alert::Notation::Number(backgroundH) << Alert::NewLine
          << diagnosticLabel("Fitting reference spacing h0:")
          << Alert::Notation::Number(parameters.model.h) << Alert::NewLine
          << diagnosticLabel("Mean interface triangle size:")
          << Alert::Notation::Number(meanInterfaceSize) << Alert::NewLine
          << diagnosticLabel("Outer iteration cap:")
          << Alert::Notation::Number(parameters.convergence.iterations.outer)
          << Alert::NewLine << diagnosticLabel("Inner correction cap:")
          << Alert::Notation::Number(parameters.convergence.iterations.inner)
          << Alert::NewLine << diagnosticLabel("Inner cap policy:")
          << "Stop on an uncertified inner residual" << Alert::NewLine
          << diagnosticLabel("Inner relative residual tolerance:")
          << Alert::Notation::Number(parameters.convergence.tolerance.innerRelative)
          << Alert::NewLine << diagnosticLabel("Barrier model:")
          << "Adaptive affine quadratic quality hinges" << Alert::NewLine
          << diagnosticLabel("Soft quality guard fraction:")
          << Alert::Notation::Number(parameters.model.qualityGuard) << Alert::NewLine
          << diagnosticLabel("Fitting / deviatoric / divergence:")
          << Alert::Notation::Number(parameters.model.fit) << " / "
          << Alert::Notation::Number(parameters.model.distribution.deviatoric) << " / "
          << Alert::Notation::Number(parameters.model.distribution.divergence)
          << Alert::NewLine << diagnosticLabel("Inactive hinge skips:")
          << Alert::Notation::Number(report.inactiveHingeSkips) << Alert::NewLine
          << diagnosticLabel("Direct symbolic analyses:")
          << Alert::Notation::Number(report.directAnalyses) << Alert::NewLine
          << diagnosticLabel("Direct numeric factorizations:")
          << Alert::Notation::Number(report.directFactorizations) << Alert::NewLine
          << diagnosticLabel("Inner Newton steps:")
          << "Full with fixed-inner merit backtracking" << Alert::NewLine
          << diagnosticLabel("Linear backend:")
          << (parameters.linear.solver == Adaptation::SWIFT::Parameters::LinearSolver::CG
                 ? "CG"
                 : (parameters.linear.solver ==
                         Adaptation::SWIFT::Parameters::LinearSolver::MUMPS
                       ? "MUMPS"
                       : "SparseLU"))
          << Alert::NewLine << diagnosticLabel("Linear relative tolerance:")
          << Alert::Notation::Number(parameters.convergence.tolerance.linearRelative)
          << Alert::NewLine << diagnosticLabel("CG iterations per solve:")
          << Alert::Notation::Number(parameters.convergence.iterations.linear)
          << Alert::NewLine << diagnosticLabel("Outer iterations:")
          << Alert::Notation::Number(report.iterations) << Alert::NewLine
          << diagnosticLabel("Inner corrections (last):")
          << Alert::Notation::Number(report.lastInnerIterations) << Alert::NewLine
          << diagnosticLabel("Inner relative correction:")
          << Alert::Notation::Number(report.innerRelativeCorrection) << Alert::NewLine
          << diagnosticLabel("Inner relative residual:")
          << Alert::Notation::Number(report.innerRelativeResidual) << Alert::NewLine
          << diagnosticLabel("Inner converged (last):")
          << (report.innerConverged ? "Yes" : "No") << Alert::NewLine
          << diagnosticLabel("Full inner Newton steps:")
          << Alert::Notation::Number(report.fullInnerSteps) << Alert::NewLine
          << diagnosticLabel("Minimum inner step factor:")
          << Alert::Notation::Number(report.minInnerAlpha) << Alert::NewLine
          << diagnosticLabel("Maximum linear iterations:")
          << Alert::Notation::Number(report.maxLinearIterations) << Alert::NewLine
          << diagnosticLabel("Skeleton normal jump RMS:")
          << Alert::Notation::Number(report.normalJumpRMS) << Alert::NewLine
          << diagnosticLabel("Exit reason:") << report.getReasonString() << Alert::NewLine
          << diagnosticLabel("Geometric RMS distance:")
          << Alert::Notation::Number(report.geometricRMS) << Alert::NewLine
          << diagnosticLabel("Geometric D infinity:")
          << Alert::Notation::Number(report.geometricSup) << Alert::NewLine
          << diagnosticLabel("Geometric D infinity target:")
          << Alert::Notation::Number(parameters.convergence.tolerance.geometric)
          << Alert::NewLine << diagnosticLabel("Geometric target reached:")
          << (report.geometricSup <= parameters.convergence.tolerance.geometric ? "Yes"
                                                                                : "No")
          << Alert::NewLine << diagnosticLabel("Welsch residual scale:")
          << Alert::Notation::Number(report.sigma) << Alert::NewLine
          << diagnosticLabel("Geometric RMS / reference spacing:")
          << Alert::Notation::Number(report.geometricRMS / parameters.model.h)
          << Alert::NewLine << diagnosticLabel("Minimum Jacobian:")
          << Alert::Notation::Number(report.minJ) << Alert::NewLine
          << diagnosticLabel("Maximum relative distortion:")
          << Alert::Notation::Number(report.maxQRel) << Alert::NewLine
          << diagnosticLabel("Setup time:") << Alert::Notation::Number(report.tSetup)
          << " s" << Alert::NewLine << diagnosticLabel("Step assembly time:")
          << Alert::Notation::Number(report.tAssembly) << " s" << Alert::NewLine
          << diagnosticLabel("Linear solve time:")
          << Alert::Notation::Number(report.tSolve) << " s" << Alert::NewLine
          << diagnosticLabel("Outer line-search time:")
          << Alert::Notation::Number(report.tLineSearch) << " s" << Alert::NewLine
          << diagnosticLabel("Linear solves:")
          << Alert::Notation::Number(report.linearSolveCount) << Alert::NewLine
          << diagnosticLabel("Total linear iterations:")
          << Alert::Notation::Number(report.linearIterations) << Alert::Raise;

        KelvinBall::Mesh moved(mesh);
        moveMesh(moved, mesh, displacementTrial.getSolution());
        checkFixedGeometry(moved, outerRadius);
        checkMaterials(moved);
        const MeshDiagnostics diagnostics = getMeshDiagnostics(moved);
        return {MMG::Mesh(std::move(moved)),
          {backgroundH, backgroundH, 0, parameters.fixedBoundaryAttributes.size(),
            diagnostics.cells, diagnostics.cells}};
      }

      void adaptSWIFT(MMGReconstruction& fitted, const Sphere& sphere,
        const Configuration& configuration, Real requestedWelschScale)
      {
        const auto before = getMeshDiagnostics(fitted.mesh);
        try
        {
          sphere.adapt(fitted.mesh, requestedWelschScale);
        }
        catch (const std::exception& error)
        {
          fitted.diagnostics.cellsBefore = before.cells;
          fitted.diagnostics.cellsAfter = before.cells;
          Alert::Warning() << "Skipped MMG adaptation after SWIFT reconstruction: "
                           << error.what() << Alert::NewLine
                           << "Retaining the SWIFT-fitted mesh with "
                           << Alert::Notation::Number(before.cells)
                           << " cells and continuing optimization." << Alert::Raise;
          return;
        }
        // MMG preserves the material partition, but may change the internal
        // face references. Gamma is the boundary between the two materials.
        fitted.mesh.getConnectivity().compute(2, 3);
        for (auto face = fitted.mesh.getFace(); face; ++face)
        {
          const auto& cells =
            fitted.mesh.getConnectivity().getIncidence({2, 3}, face->getIndex());
          if (cells.size() == 2 &&
            fitted.mesh.getPolytope(3, cells[0])->getAttribute() !=
              fitted.mesh.getPolytope(3, cells[1])->getAttribute())
            fitted.mesh.setAttribute({2, face->getIndex()}, Gamma);
        }
        const auto after = getMeshDiagnostics(fitted.mesh);
        fitted.diagnostics.minimumSize = configuration.hmin;
        fitted.diagnostics.maximumSize = configuration.hmax;
        fitted.diagnostics.hausdorffTolerance =
          Real(0.1) * configuration.getGridSpacing() * configuration.getGridSpacing();
        fitted.diagnostics.cellsBefore = before.cells;
        fitted.diagnostics.cellsAfter = after.cells;
        fitted.diagnostics.requiredBoundaryTriangles =
          sphere.protectFixedGeometry(fitted.mesh, false);
        Alert::Info() << substageHeading("MMG adaptation after SWIFT reconstruction")
                      << Alert::NewLine << diagnosticLabel("Cell count:")
                      << Alert::Notation::Number(before.cells) << " -> "
                      << Alert::Notation::Number(after.cells) << Alert::NewLine
                      << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(before.interfaceTriangles) << " -> "
                      << Alert::Notation::Number(after.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(before.minimumQuality) << " -> "
                      << Alert::Notation::Number(after.minimumQuality) << Alert::NewLine
                      << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(before.meanQuality) << " -> "
                      << Alert::Notation::Number(after.meanQuality) << Alert::NewLine
                      << diagnosticLabel("Mean tetra edge length h:")
                      << Alert::Notation::Number(before.meanElementSize) << " -> "
                      << Alert::Notation::Number(after.meanElementSize) << Alert::Raise;
      }

      void reportConfiguration(Real requestedWelschScale)
      {
        const size_t points = m_options.configuration.points;
        const Real outerRadius = m_options.configuration.outerRadius;
        const Real nitschePenalty = m_options.configuration.nitschePenalty;
        const Real stabilizationFactor = m_options.configuration.stabilizationFactor;
        const Real initialGridSpacing = m_options.configuration.getGridSpacing();
        Alert::Info configurationInfo;
        configurationInfo << substageHeading("Configuration") << Alert::NewLine
                          << diagnosticLabel("Grid points:")
                          << Alert::Notation::Number(points) << Alert::NewLine
                          << diagnosticLabel("Outer radius:")
                          << Alert::Notation::Number(outerRadius) << Alert::NewLine
                          << diagnosticLabel("Reference spacing h0:")
                          << Alert::Notation::Number(initialGridSpacing) << Alert::NewLine
                          << diagnosticLabel("Requested minimum size:")
                          << Alert::Notation::Number(m_options.configuration.hmin)
                          << Alert::NewLine << diagnosticLabel("Requested maximum size:")
                          << Alert::Notation::Number(m_options.configuration.hmax)
                          << Alert::NewLine << diagnosticLabel("Evaluated designs:")
                          << Alert::Notation::Number(m_options.maxIterations)
                          << Alert::NewLine << diagnosticLabel("Nitsche penalty:")
                          << Alert::Notation::Number(nitschePenalty) << Alert::NewLine
                          << diagnosticLabel("Stabilization factor:")
                          << Alert::Notation::Number(stabilizationFactor)
                          << Alert::NewLine << diagnosticLabel("Assembly backend:")
                          << assemblyBackend << Alert::NewLine
                          << diagnosticLabel("Direct solver:")
                          << KelvinBall::DirectSolverName << Alert::NewLine
                          << diagnosticLabel("Regularization length:")
                          << Alert::Notation::Number(m_options.regularizationFactor)
                          << " h0" << Alert::NewLine
                          << diagnosticLabel("Thickness model:")
                          << "Paired first-exit surrogate with aggregate correction"
                          << Alert::NewLine << diagnosticLabel("Normal smoothing length:")
                          << Alert::Notation::Number(m_options.normalRegularizationFactor)
                          << " h0" << Alert::NewLine << diagnosticLabel("Advection step:")
                          << Alert::Notation::Number(m_options.stepFactor) << " h0"
                          << Alert::NewLine << diagnosticLabel("Level-set trace penalty:")
                          << Alert::Notation::Number(m_options.levelSetPenalty)
                          << Alert::NewLine
                          << diagnosticLabel("Advection quadrature order:")
                          << Alert::Notation::Number(m_options.advectionQuadratureOrder)
                          << Alert::NewLine << diagnosticLabel("Transport velocity:")
                          << "Full shape direction" << Alert::NewLine
                          << diagnosticLabel("Reconstruction method:")
                          << m_options.reconstructionMethod;
        if (m_options.reconstructionMethod == "swift" && !m_options.configuration.adapt)
        {
          configurationInfo << Alert::NewLine << diagnosticLabel("Background Hausdorff:")
                            << Alert::Notation::Number(
                                 m_options.configuration.backgroundHausdorff)
                            << " h0" << Alert::NewLine
                            << diagnosticLabel("Background gradation:")
                            << Alert::Notation::Number(
                                 m_options.configuration.backgroundGradation);
        }
        if (m_options.configuration.adapt)
        {
          configurationInfo << Alert::NewLine
                            << diagnosticLabel("Adaptation hmin (interface):")
                            << Alert::Notation::Number(m_options.configuration.hmin)
                            << Alert::NewLine
                            << diagnosticLabel("Adaptation hmax (far field):")
                            << Alert::Notation::Number(m_options.configuration.hmax);
          if (m_options.reconstructionMethod == "swift")
            configurationInfo << Alert::NewLine
                              << diagnosticLabel("Background Hausdorff:")
                              << Alert::Notation::Number(
                                   m_options.configuration.backgroundHausdorff)
                              << " h0";
          if (m_options.reconstructionMethod == "swift")
            configurationInfo << Alert::NewLine
                              << diagnosticLabel("Welsch size-map scale:")
                              << (requestedWelschScale > 0
                                     ? "User-specified"
                                     : "3 h0 (fixed reference spacing)");
          configurationInfo << Alert::NewLine << diagnosticLabel("Adaptation gradation:")
                            << Alert::Notation::Number(
                                 m_options.configuration.adaptGradation);
        }
        if (m_options.configuration.mmgSnap > 0)
          configurationInfo << Alert::NewLine
                            << diagnosticLabel("Level-set snapping fraction:")
                            << Alert::Notation::Number(m_options.configuration.mmgSnap);
        configurationInfo << Alert::Raise;
      }

      struct InitialDesign
      {
          MMG::Mesh mesh;
          Optional<MMG::Mesh> background;
          Real backgroundH;
          ReconstructionDiagnostics diagnostics;
      };

      InitialDesign initialize(const Sphere& sphere, Real requestedWelschScale)
      {
        const Real initialGridSpacing = m_options.configuration.getGridSpacing();
        const Real outerRadius = m_options.configuration.outerRadius;
        const int argc = m_argc;
        char** argv = m_argv;
        SphereDiscretization initial = m_options.reconstructionMethod == "swift"
          ? sphere.prepareSWIFTBackground(requestedWelschScale)
          : sphere.discretize(false, requestedWelschScale);
        ReconstructionDiagnostics reconstruction = initial.diagnostics;
        Optional<MMG::Mesh> swiftBackground;
        Real backgroundH = std::numeric_limits<Real>::quiet_NaN();
        MMG::Mesh mesh;
        if (m_options.reconstructionMethod == "swift")
        {
          swiftBackground.emplace(std::move(initial.mesh));
          const MeshDiagnostics backgroundDiagnostics =
            getMeshDiagnostics(*swiftBackground, false);
          backgroundH = backgroundDiagnostics.meanElementSize;
          Alert::Info backgroundInfo;
          backgroundInfo
            << substageHeading(m_options.configuration.adapt
                   ? "MMG background adaptation"
                   : "MMG background optimization")
            << Alert::NewLine << diagnosticLabel("Cell count:")
            << Alert::Notation::Number(reconstruction.cellsBefore) << " -> "
            << Alert::Notation::Number(backgroundDiagnostics.cells) << Alert::NewLine
            << diagnosticLabel("Required boundary triangles:")
            << Alert::Notation::Number(reconstruction.requiredBoundaryTriangles)
            << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
            << Alert::Notation::Number(backgroundDiagnostics.minimumQuality)
            << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
            << Alert::Notation::Number(backgroundDiagnostics.meanQuality)
            << Alert::NewLine << diagnosticLabel("Mean tetra edge length h:")
            << Alert::Notation::Number(backgroundDiagnostics.meanElementSize);
          if (m_options.configuration.adapt)
            backgroundInfo << Alert::NewLine << diagnosticLabel("Input chamber h:")
                           << Alert::Notation::Number(
                                reconstruction.backgroundMeanElementSize)
                           << Alert::NewLine << diagnosticLabel("Welsch size-map scale:")
                           << Alert::Notation::Number(reconstruction.welschScale);
          backgroundInfo << Alert::Raise;
          P1 sphereSpace(*swiftBackground);
          GridFunction sphereLevelSet(sphereSpace);
          sphereLevelSet = RealFunction([](const Geometry::Point& point) {
            return point.getPhysicalCoordinates().norm() - Real(1);
          });
          MMG::Mesh classified =
            classifyLevelSetForSWIFT(*swiftBackground, sphereLevelSet);
          P1 classifiedSphereSpace(classified);
          GridFunction classifiedSphereLevelSet(classifiedSphereSpace);
          classifiedSphereLevelSet.getData() = sphereLevelSet.getData();
          MMGReconstruction fitted =
            fitLevelSetSWIFT(classified, classifiedSphereLevelSet, backgroundH,
              initialGridSpacing, outerRadius, argc, argv);
          if (m_options.configuration.adapt)
          {
            adaptSWIFT(fitted, sphere, m_options.configuration, requestedWelschScale);
            reconstruction = fitted.diagnostics;
          }
          mesh = std::move(fitted.mesh);
        }
        else
        {
          mesh = std::move(initial.mesh);
        }
        if (swiftBackground && m_options.configuration.adapt)
        {
          swiftBackground.emplace(mesh);
          backgroundH = meanElementSize(*swiftBackground);
        }
        return {std::move(mesh), std::move(swiftBackground), backgroundH, reconstruction};
      }

      static void printUsage(const char* executable)
      {
        Alert::Info()
          << "Usage" << Alert::NewLine << "  " << executable << " [options]"
          << Alert::NewLine << Alert::Notation("--n=<points>")
          << "              Background points per edge (default: 13)." << Alert::NewLine
          << Alert::Notation("--h=<size>")
          << "                Initial grid spacing; alternative to --n." << Alert::NewLine
          << Alert::Notation("--hmin-factor=<value>")
          << "       Minimum MMG size / reference spacing h0 (default: 0.1)."
          << Alert::NewLine << Alert::Notation("--hmax-factor=<value>")
          << "       Maximum MMG size / reference spacing h0 (default: 10)."
          << Alert::NewLine << Alert::Notation("--outer-radius=<value>")
          << "     Chamber outer radius (default: 2)." << Alert::NewLine
          << Alert::Notation("--iterations=<count>")
          << "       Number of evaluated designs (default: 1)." << Alert::NewLine
          << Alert::Notation("--penalty=<value>")
          << "          Rotational Nitsche penalty (default: 320)." << Alert::NewLine
          << Alert::Notation("--stabilization=<value>")
          << "    P1--P1 pressure-stabilization factor (default: 0.05)." << Alert::NewLine
          << Alert::Notation("--regularization=<value>")
          << "   H1 smoothing length in multiples of h0 (default: 4)." << Alert::NewLine
          << Alert::Notation("--thickness-diagnostics")
          << "       Export diagnostic smoothed normals and curvature (default: off)."
          << Alert::NewLine << Alert::Notation("--normal-regularization=<value>")
          << " Thickness-normal smoothing length in h0 (default: 1)." << Alert::NewLine
          << Alert::Notation("--step=<value>")
          << "             Advection time step in multiples of h0 (default: 0.1)."
          << Alert::NewLine << Alert::Notation("--level-set-penalty=<value>")
          << " Rotated trace penalty of the level set (default: 1)." << Alert::NewLine
          << Alert::Notation("--thickness-min=<value>")
          << "      Absolute minimum body thickness (default: 0, off)." << Alert::NewLine
          << Alert::Notation("--motion-force=<fx,fy,fz>")
          << "  Force for the Motion field at each iterate (default: 1,0,0)."
          << Alert::NewLine << Alert::Notation("--advection-quadrature=<order>")
          << " Quadrature order for the transported distance (default: 8)."
          << Alert::NewLine << Alert::Notation("--reconstruction=<method>")
          << "  Interface reconstruction: mmg or swift (default: mmg)." << Alert::NewLine
          << Alert::Notation("--background-hausdorff=<value>")
          << " SWIFT background Hausdorff tolerance, in h0 (default: 0.05)."
          << Alert::NewLine << Alert::Notation("--background-gradation=<value>")
          << " SWIFT background gradation (default: 2)." << Alert::NewLine
          << Alert::Notation("--mmg-adapt")
          << "                Adapt near the interface after each reconstruction:"
          << Alert::NewLine
          << "                              MMG cut or SWIFT fit, including the initial "
             "design."
          << Alert::NewLine
          << "                              Adaptation uses the fixed size factors times "
             "h0."
          << Alert::NewLine << Alert::Notation("--mmg-adapt-gradation=<value>")
          << "  Adaptation gradation (default: 1.3)." << Alert::NewLine
          << Alert::Notation("--mmg-snap=<value>")
          << "         MMG path: snap edge crossings closer than this fraction"
          << Alert::NewLine
          << "                              of the edge to a vertex (default: 0, off)."
          << Alert::NewLine << Alert::Notation("--mmg-retries=<count>")
          << "      MMG path: retries of a failed reconstruction, each at half"
          << Alert::NewLine
          << "                              the previous MMG scale (default: 2)."
          << Alert::NewLine << Alert::Notation("--swift-*=<value>")
          << "        SWIFT fitting parameters (defaults: 30 outer / 15 inner;"
          << Alert::NewLine
          << "                              iteration trace on; --swift-trace=0 disables;"
          << Alert::NewLine
          << "                              grouped model and convergence controls)."
          << Alert::NewLine << Alert::Notation("--swift-classification-fidelity=<value>")
          << "  MinSTCut fidelity (default: 1)." << Alert::NewLine
          << Alert::Notation("--swift-classification-smoothing=<value>")
          << "  MinSTCut dimensionless smoothing (default: 0.04);" << Alert::NewLine
          << "                              facet weight = value times smaller incident "
             "cell size."
          << Alert::NewLine << Alert::Notation("--swift-model-fit=<value>")
          << "  Fitting curvature weight (default: 1)." << Alert::NewLine
          << Alert::Notation("--swift-model-distribution-deviatoric=<value>")
          << "  Centered deviatoric weight (default: 0.0001)." << Alert::NewLine
          << Alert::Notation("--swift-model-distribution-divergence=<value>")
          << "  Centered divergence weight (default: 0.01)." << Alert::NewLine
          << Alert::Notation("--swift-linear-solver=<name>")
          << "  SWIFT linear backend: mumps, sparse-lu, or cg (default: MUMPS when "
             "built)."
          << Alert::NewLine
          << Alert::Notation("--swift-globalization-directional-newton[=0|1]")
          << "  Scale the frozen model with directional Newton (default: 1)."
          << Alert::NewLine << Alert::Notation("--swift-model-hinge=<value>")
          << "  Dimensionless quadratic-hinge weight (default: 10)." << Alert::NewLine
          << Alert::Notation("--swift-model-quality-guard=<fraction>")
          << "  Soft guard fraction for that penalty (default: 0.1)." << Alert::NewLine
          << Alert::Notation("--geometry-only")
          << "             Stop after initial reconstruction." << Alert::NewLine
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
          for (size_t i = 0; i < vertices.size(); ++i)
          {
            for (size_t j = i + 1; j < vertices.size(); ++j)
            {
              const Real edgeLength = (mesh.getVertexCoordinates(vertices[i]) -
                mesh.getVertexCoordinates(vertices[j]))
                                        .norm();
              squaredEdgeLengthSum += edgeLength * edgeLength;
            }
          }
          const Real quality =
            12 * std::pow(3 * cell->getMeasure(), Real(2) / 3) / squaredEdgeLengthSum;
          diagnostics.minimumQuality = std::min(diagnostics.minimumQuality, quality);
          diagnostics.maximumQuality = std::max(diagnostics.maximumQuality, quality);
          qualitySum += quality;
        }
        diagnostics.meanQuality = qualitySum / static_cast<Real>(diagnostics.cells);
        diagnostics.meanElementSize = meanElementSize(mesh);
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
                      << diagnosticLabel("Solid cells:")
                      << Alert::Notation::Number(diagnostics.obstacleCells)
                      << Alert::NewLine << diagnosticLabel("Fluid cells:")
                      << Alert::Notation::Number(diagnostics.fluidCells) << Alert::NewLine
                      << diagnosticLabel("Total cells:")
                      << Alert::Notation::Number(diagnostics.cells) << Alert::NewLine
                      << diagnosticLabel("Interface triangles:")
                      << Alert::Notation::Number(diagnostics.interfaceTriangles)
                      << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.minimumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.meanQuality)
                      << Alert::NewLine << diagnosticLabel("Maximum tetrahedron quality:")
                      << Alert::Notation::Number(diagnostics.maximumQuality)
                      << Alert::NewLine << diagnosticLabel("Mean tetra edge length h:")
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

      template <class LevelSet>
      MMGReconstruction discretizeLevelSetMMG(MMG::Mesh& mesh, const LevelSet& levelSet,
        Real targetSize, Real hmin, Real hmax, const Sphere& sphere, bool adapt,
        Real snap, Real requestedWelschScale)
      {
        const size_t previousCells = mesh.getCellCount();
        const Real hausdorff = 0.1 * targetSize * targetSize;
        // A retried reconstruction receives the mesh the failed attempt already
        // reduced to a single material, so the input need not be partitioned.
        const MeshDiagnostics inputDiagnostics = getMeshDiagnostics(mesh, false);

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
                      << Alert::NewLine << diagnosticLabel("Mean tetra edge length h:")
                      << Alert::Notation::Number(inputDiagnostics.meanElementSize)
                      << Alert::Raise;

        // The level set alone defines the new partition. MMG keeps the faces
        // between differently labelled cells as internal boundaries through
        // the cut, where they end up inside one material; the optimization
        // pass (opnbdy) then preserves them, so every past interface would
        // accumulate in the mesh. The cut therefore starts, as the initial
        // one does, from a single material and only the fixed boundary.
        const FlatSet<Attribute> fixed{
          Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
        for (auto cell = mesh.getCell(); cell; ++cell)
          mesh.setAttribute({mesh.getDimension(), cell->getIndex()}, Fluid);
        for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (attribute && !fixed.contains(*attribute))
            mesh.setAttribute({mesh.getDimension() - 1, face->getIndex()}, {});
        }
        const size_t requiredTriangles = sphere.protectFixedGeometry(mesh, false);

        MMG::LevelSetDiscretizer discretizer;
        // Under RMC, retain the body component attached to the chamber cuts;
        // Fluid is a cell label, not a boundary-face base reference.
        discretizer.split(Fluid, {Obstacle, Fluid})
          .setHMin(hmin)
          .setHMax(hmax)
          .setHausdorff(hausdorff)
          .setGradation(remeshGradation)
          .setBaseReferences(
            FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus})
          .setBoundaryReference(Gamma)
          .setRMC(1e-5)
          .setAngleDetection(false);
        // A crossed edge (i, j) is cut at t = phi_i / (phi_i - phi_j). A cut
        // with t near 0 or 1 puts the new vertex next to an existing one and
        // leaves a sliver. Snapping the near endpoint to zero instead makes the
        // interface pass through that vertex, so every remaining cut satisfies
        // snap <= t <= 1 - snap. For a distance-like level set the interface
        // moves by at most snap times the edge length.
        std::decay_t<LevelSet> sanitized(levelSet.getFiniteElementSpace());
        sanitized.getData() = levelSet.getData();
        const auto crossingFraction = [&](Index i, Index j) {
          const Real a = sanitized[i], b = sanitized[j];
          return (a < 0) != (b < 0) && a != 0 && b != 0
            ? a / (a - b)
            : std::numeric_limits<Real>::quiet_NaN();
        };
        const auto scanCrossings = [&]() {
          size_t edges = 0, nearVertex = 0;
          Real minimum = 0.5;
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& vertices = cell->getVertices();
            for (size_t a = 0; a < vertices.size(); ++a)
              for (size_t b = a + 1; b < vertices.size(); ++b)
              {
                const Real t = crossingFraction(vertices[a], vertices[b]);
                if (std::isnan(t))
                  continue;
                ++edges;
                const Real distanceToVertex = std::min(t, 1 - t);
                minimum = std::min(minimum, distanceToVertex);
                if (distanceToVertex < Real(1e-3))
                  ++nearVertex;
              }
          }
          return std::tuple{edges, nearVertex, minimum};
        };
        const auto [cutEdges, nearVertexCuts, minimumCrossing] = scanCrossings();
        Alert::Info crossingInfo;
        crossingInfo << substageHeading("Level-set crossings") << Alert::NewLine
                     << diagnosticLabel("Crossed cell edges:")
                     << Alert::Notation::Number(cutEdges) << Alert::NewLine
                     << diagnosticLabel("Crossings within 1e-3 of a vertex:")
                     << Alert::Notation::Number(nearVertexCuts) << Alert::NewLine
                     << diagnosticLabel("Minimum crossing fraction:")
                     << Alert::Notation::Number(minimumCrossing);
        size_t snappedCount = 0;
        if (snap > 0)
        {
          std::vector<Real> original(
            sanitized.getData().begin(), sanitized.getData().end());
          std::vector<char> snapped(mesh.getVertexCount(), 0);
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& vertices = cell->getVertices();
            for (size_t a = 0; a < vertices.size(); ++a)
              for (size_t b = a + 1; b < vertices.size(); ++b)
              {
                const Index i = vertices[a], j = vertices[b];
                const Real t = original[i] * original[j] < 0
                  ? original[i] / (original[i] - original[j])
                  : Real(0.5);
                if (t < snap)
                  snapped[i] = 1;
                else if (t > 1 - snap)
                  snapped[j] = 1;
              }
          }
          for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
            if (snapped[vertex])
              sanitized[vertex] = 0;
          // MMG requires every tetrahedron to keep a vertex off the level set;
          // the vertex farthest from it keeps its value.
          size_t restored = 0;
          for (auto cell = mesh.getCell(); cell; ++cell)
          {
            const auto& vertices = cell->getVertices();
            Index farthest = vertices[0];
            bool allZero = true;
            for (const Index vertex : vertices)
            {
              allZero = allZero && sanitized[vertex] == 0;
              if (std::abs(original[vertex]) > std::abs(original[farthest]))
                farthest = vertex;
            }
            if (allZero)
            {
              sanitized[farthest] = original[farthest];
              ++restored;
            }
          }
          Real displacement = 0;
          size_t snappedVertices = 0;
          for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
          {
            if (sanitized[vertex] == 0 && original[vertex] != 0)
            {
              ++snappedVertices;
              ++snappedCount;
              displacement = std::max(displacement, std::abs(original[vertex]));
            }
          }
          const auto [snappedEdges, snappedNearVertex, snappedMinimum] = scanCrossings();
          crossingInfo << Alert::NewLine << diagnosticLabel("Snapping fraction:")
                       << Alert::Notation::Number(snap) << Alert::NewLine
                       << diagnosticLabel("Snapped vertices:")
                       << Alert::Notation::Number(snappedVertices) << Alert::NewLine
                       << diagnosticLabel("Restored in all-zero cells:")
                       << Alert::Notation::Number(restored) << Alert::NewLine
                       << diagnosticLabel("Maximum snapped level set:")
                       << Alert::Notation::Number(displacement) << " = "
                       << Alert::Notation::Number(displacement / targetSize)
                       << " target sizes" << Alert::NewLine
                       << diagnosticLabel("Minimum crossing after snapping:")
                       << Alert::Notation::Number(snappedMinimum);
        }
        crossingInfo << Alert::Raise;
        MMG::Mesh reconstructed = discretizer.discretize(sanitized);
        sphere.projectFixedGeometry(reconstructed);

        splitSelfPairedCut(reconstructed);
        const MeshDiagnostics reconstructionDiagnostics =
          getMeshDiagnostics(reconstructed, false);
        Alert::Info()
          << substageHeading("MMG reconstruction") << Alert::NewLine
          << diagnosticLabel("Minimum size:") << Alert::Notation::Number(hmin)
          << Alert::NewLine << diagnosticLabel("Maximum size:")
          << Alert::Notation::Number(hmax) << Alert::NewLine
          << diagnosticLabel("Hausdorff tolerance:") << Alert::Notation::Number(hausdorff)
          << Alert::NewLine << diagnosticLabel("Required boundary triangles:")
          << Alert::Notation::Number(requiredTriangles) << Alert::NewLine
          << diagnosticLabel("Cell count:") << Alert::Notation::Number(previousCells)
          << " -> " << Alert::Notation::Number(reconstructionDiagnostics.cells)
          << Alert::NewLine << diagnosticLabel("Solid cells:")
          << Alert::Notation::Number(reconstructionDiagnostics.obstacleCells)
          << Alert::NewLine << diagnosticLabel("Fluid cells:")
          << Alert::Notation::Number(reconstructionDiagnostics.fluidCells)
          << Alert::NewLine << diagnosticLabel("Total cells:")
          << Alert::Notation::Number(reconstructionDiagnostics.cells) << Alert::NewLine
          << diagnosticLabel("Interface triangles:")
          << Alert::Notation::Number(inputDiagnostics.interfaceTriangles) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.interfaceTriangles)
          << Alert::NewLine << diagnosticLabel("Minimum tetrahedron quality:")
          << Alert::Notation::Number(inputDiagnostics.minimumQuality) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.minimumQuality)
          << Alert::NewLine << diagnosticLabel("Mean tetrahedron quality:")
          << Alert::Notation::Number(inputDiagnostics.meanQuality) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.meanQuality)
          << Alert::NewLine << diagnosticLabel("Mean tetra edge length h:")
          << Alert::Notation::Number(inputDiagnostics.meanElementSize) << " -> "
          << Alert::Notation::Number(reconstructionDiagnostics.meanElementSize)
          << Alert::Raise;

        const size_t requiredTrianglesBeforeOptimization =
          sphere.protectFixedGeometry(reconstructed, false);
        if (adapt)
        {
          sphere.adapt(reconstructed, requestedWelschScale);
        }
        else
        {
          MMG::Optimizer()
            .setHMin(hmin)
            .setHMax(hmax)
            .setHausdorff(hausdorff)
            .setGradation(remeshGradation)
            .setAngleDetection(false)
            .optimize(reconstructed);
          sphere.projectFixedGeometry(reconstructed);
          splitSelfPairedCut(reconstructed);
        }
        const size_t requiredTrianglesAfterOptimization =
          sphere.protectFixedGeometry(reconstructed, false);
        const MeshDiagnostics outputDiagnostics =
          getMeshDiagnostics(reconstructed, false);
        Alert::Info()
          << substageHeading(adapt ? "MMG adaptation" : "MMG optimization")
          << Alert::NewLine << diagnosticLabel("Required boundary triangles:")
          << Alert::Notation::Number(requiredTrianglesBeforeOptimization) << " -> "
          << Alert::Notation::Number(requiredTrianglesAfterOptimization) << Alert::NewLine
          << diagnosticLabel("Cell count:")
          << Alert::Notation::Number(reconstructionDiagnostics.cells) << " -> "
          << Alert::Notation::Number(outputDiagnostics.cells) << Alert::NewLine
          << diagnosticLabel("Minimum tetrahedron quality:")
          << Alert::Notation::Number(reconstructionDiagnostics.minimumQuality) << " -> "
          << Alert::Notation::Number(outputDiagnostics.minimumQuality) << Alert::NewLine
          << diagnosticLabel("Mean tetrahedron quality:")
          << Alert::Notation::Number(reconstructionDiagnostics.meanQuality) << " -> "
          << Alert::Notation::Number(outputDiagnostics.meanQuality) << Alert::NewLine
          << diagnosticLabel("Mean tetra edge length h:")
          << Alert::Notation::Number(reconstructionDiagnostics.meanElementSize) << " -> "
          << Alert::Notation::Number(outputDiagnostics.meanElementSize) << Alert::Raise;
        ReconstructionDiagnostics diagnostics{hmin, hmax, hausdorff,
          requiredTrianglesAfterOptimization, previousCells, outputDiagnostics.cells};
        diagnostics.minimumCrossing = minimumCrossing;
        diagnostics.snappedVertices = static_cast<Real>(snappedCount);
        return {std::move(reconstructed), diagnostics};
      }

      template <class Density>
      GradientDiagnostics identifyGradient(HilbertIdentification& identification,
        const Density& density, HilbertIdentification::Field& gradient, const char* name,
        const Math::Vector<Real>* additionalLoad = nullptr)
      {
        const auto diagnostics =
          identification.identify(density, gradient, additionalLoad);
        Alert::Info() << substageHeading(std::string(name) + " gradient")
                      << Alert::NewLine << diagnosticLabel("Linear residual:")
                      << Alert::Notation::Number(diagnostics.residual) << Alert::NewLine
                      << diagnosticLabel("Rotated jump:")
                      << Alert::Notation::Number(diagnostics.jump) << Alert::Raise;
        return diagnostics;
      }

      int m_argc;
      char** m_argv;
  };
}

KelvinBall::KelvinBallOptimization::KelvinBallOptimization(int argc, char** argv)
  : m_implementation(std::make_unique<Implementation>(argc, argv))
{}

KelvinBall::KelvinBallOptimization::~KelvinBallOptimization() = default;

int KelvinBall::KelvinBallOptimization::Implementation::run()
{
  const int argc = m_argc;
  char** argv = m_argv;
  if (!m_options.parse(argc, argv))
    return 0;
  const Real outerRadius = m_options.configuration.outerRadius;
  const Real nitschePenalty = m_options.configuration.nitschePenalty;
  const Real stabilizationFactor = m_options.configuration.stabilizationFactor;
  const Real initialGridSpacing = m_options.configuration.getGridSpacing();
  const Real requestedWelschScale =
    getFittingParameters(m_options.configuration.getGridSpacing(), argc, argv).model.robustScale;
  reportConfiguration(requestedWelschScale);
  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const auto stage1Start = Clock::now();
  announce(m_options.reconstructionMethod == "swift"
      ? "Stage 1: Preparing the background mesh and fitting the initial sphere with "
        "SWIFT."
      : "Stage 1: Discretizing the initial sphere with MMG.");
  Sphere sphere(m_options.configuration);
  auto initial = initialize(sphere, requestedWelschScale);
  MMG::Mesh mesh = std::move(initial.mesh);
  Optional<MMG::Mesh> swiftBackground = std::move(initial.background);
  Real backgroundH = initial.backgroundH;
  ReconstructionDiagnostics reconstruction = initial.diagnostics;
  Optional<IO::XDMF> reconstructionXdmf;
  Optional<IO::XDMF::Grid> reconstructionOutput;
  if (m_options.reconstructionMethod == "mmg")
  {
    reconstructionXdmf.emplace("KelvinBallMMG");
    reconstructionOutput.emplace(reconstructionXdmf->grid("Reconstructed"));
    reconstructionOutput->setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
    reconstructionXdmf->write(Real(0)).flush();
  }
  checkFixedGeometry(mesh, outerRadius);
  checkMaterials(mesh);
  if (m_options.saveMeshDiagnostic)
    mesh.save("KelvinBallInitial.mesh", IO::FileFormat::MEDIT);
  reportSphereGeometry(mesh);
  const Real stage1Seconds = elapsedSeconds(stage1Start);
  reportStageTiming(1, stage1Seconds);
  if (m_options.geometryOnly)
    return 0;
  const Real targetVolume = mesh.getVolume(Obstacle);

  std::ofstream history("kelvin-ball.csv");
  history.precision(17);
  history << "iteration,outer_radius,initial_grid_spacing,hmin_factor,hmax_factor,"
             "requested_hmin,requested_hmax,"
             "h,dt,advection_quadrature_order,"
             "regularization_length,nitsche_penalty,"
             "stabilization_factor,assembly_backend,direct_solver,"
             "vertices,cells,obstacle_cells,fluid_cells,"
             "interface_triangles,mesh_quality_min,mesh_quality_mean,mesh_quality_max,"
             "background_h,"
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
             "stage_9_seconds,"
             "translational_mobility,coupling_mobility,rotational_mobility,"
             "revolution_period,pitch,resistance_radius,"
             "level_set_penalty,thickness_min,normal_smoothing_length,"
             "thickness_penalty,thickness_nominal_penalty,thickness_guard_distance,"
             "thickness_active_multiplier,thickness_active_feasible,"
             "thickness_violating_rays,"
             "thickness_maximum_pair_load,thickness_top_one_percent_load_share,"
             "thickness_descent_normal_rms,thickness_descent_tangent_rms,"
             "thickness_descent_tangent_fraction,theta_tangent_fraction,"
             "thickness_minimum_exit,"
             "thickness_maximum_deficit,thickness_minimum_transversality,"
             "thickness_normal_alignment_min,thickness_normal_alignment_mean,"
             "thickness_curvature_min,thickness_curvature_max,"
             "thickness_normal_magnitude_min,thickness_normal_magnitude_max,"
             "thickness_geometric_seam_jump,thickness_normal_seam_jump,"
             "thickness_ray_seam_jump,"
             "thickness_seam_samples,thickness_seam_unmatched,"
             "actual_delta_thickness,predicted_delta_thickness,d_thickness_theta,"
             "eikonal_rotated_jump,projected_rotated_jump,distance_correction,"
             "interface_shift_max,advection_increment,advected_rotated_jump,"
             "min_crossing_fraction,snapped_vertices,reconstruction_scale,"
             "z_x,z_y,z_z,omega_x,omega_y,omega_z,"
             "transport_attempted_points,transport_omitted_points,"
             "transport_omitted_weight_fraction,transport_fallback_cells,"
             "motion_force_x,motion_force_y,motion_force_z\n";
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
  Optional<Real> previousThicknessPenalty;
  Optional<Real> predictedThicknessChange;
  Optional<Real> predictedVolumeChange;
  Optional<MMG::Mesh> nextMesh;

  for (size_t iteration = 0; iteration < m_options.maxIterations; ++iteration)
  {
    std::array<Real, 9> stageSeconds;
    stageSeconds.fill(nan);
    if (iteration == 0)
      stageSeconds[0] = stage1Seconds;
    if (nextMesh)
    {
      mesh = std::move(*nextMesh);
      nextMesh.reset();
      if (swiftBackground && m_options.configuration.adapt)
      {
        // Refresh only between iterations, after the previous spaces and
        // fields have expired. The background is fixed during each fit.
        swiftBackground.emplace(mesh);
        backgroundH = meanElementSize(*swiftBackground);
      }
    }
    const Real h = meanElementSize(mesh);
    const auto meshDiagnostics = getMeshDiagnostics(mesh);
    const Real hilbertLength = m_options.regularizationFactor * initialGridSpacing;
    const Real normalLength = m_options.normalRegularizationFactor * initialGridSpacing;
    const Real dt = m_options.stepFactor * initialGridSpacing;
    Alert::Info() << " ------------------------------------------------------------"
                  << Alert::NewLine << " Iteration "
                  << Alert::Notation::Number(iteration + 1) << " of "
                  << Alert::Notation::Number(m_options.maxIterations) << Alert::NewLine
                  << " ------------------------------------------------------------"
                  << Alert::Raise;
    Alert::Info() << substageHeading("Current mesh scale") << Alert::NewLine
                  << diagnosticLabel("Mean tetra edge length h:")
                  << Alert::Notation::Number(h) << Alert::NewLine
                  << diagnosticLabel("Regularization length:")
                  << Alert::Notation::Number(hilbertLength) << Alert::NewLine
                  << diagnosticLabel("Normal smoothing length:")
                  << Alert::Notation::Number(normalLength) << Alert::NewLine
                  << diagnosticLabel("Advection step:") << Alert::Notation::Number(dt)
                  << Alert::Raise;
    const auto stage2Start = Clock::now();
    announce("Stage 2: Preparing the chamber fluid mesh.");
    Alert::Info() << substageHeading("Chamber cell counts") << Alert::NewLine
                  << diagnosticLabel("Solid cells:")
                  << Alert::Notation::Number(meshDiagnostics.obstacleCells)
                  << Alert::NewLine << diagnosticLabel("Fluid cells:")
                  << Alert::Notation::Number(meshDiagnostics.fluidCells) << Alert::NewLine
                  << diagnosticLabel("Total cells:")
                  << Alert::Notation::Number(meshDiagnostics.cells) << Alert::Raise;
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

    VelocitySpace Vh(fluid, 3);
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
    const Real resistanceDeterminant = k * q - c * c;
    const Math::SpatialVector<Real> motionTranslation =
      (q / resistanceDeterminant) * m_options.motionForce;
    const Math::SpatialVector<Real> motionAngular =
      (-c / resistanceDeterminant) * m_options.motionForce;
    const Real rho = metrics.rho;
    const Real couplingSymmetry = metrics.couplingSymmetry;
    const Real volume = mesh.getVolume(Obstacle);
    const Real nitscheJump = metrics.nitscheJump;
    const auto chamberBoundingBox = boundingBoxSize(mesh);
    const auto sewnBoundingBox = sewnBoundingBoxSize(mesh);
    const Real actualRhoChange = previousRho ? rho - *previousRho : nan;
    const Real incomingPredictedRhoChange =
      predictedRhoChange ? *predictedRhoChange : nan;
    const Real actualVolumeChange = previousVolume ? volume - *previousVolume : nan;
    const Real incomingPredictedVolumeChange =
      predictedVolumeChange ? *predictedVolumeChange : nan;
    // Diagnostics of the later stages, filled as they run.
    struct
    {
        Real thicknessPenalty = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNominalPenalty = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessGuardDistance = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessActiveMultiplier = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessActiveFeasible = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessViolating = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessMaximumPairLoad = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessTopOnePercentLoadShare = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessDescentNormalRMS = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessDescentTangentRMS = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessDescentTangentFraction = std::numeric_limits<Real>::quiet_NaN();
        Real thetaTangentFraction = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessMinimumExit = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessMaximumDeficit = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessMinimumTransversality = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNormalAlignmentMin = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNormalAlignmentMean = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessCurvatureMin = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessCurvatureMax = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNormalMagnitudeMin = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNormalMagnitudeMax = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessGeometricSeamJump = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessNormalSeamJump = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessRaySeamJump = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessSeamSamples = std::numeric_limits<Real>::quiet_NaN();
        Real thicknessSeamUnmatched = std::numeric_limits<Real>::quiet_NaN();
        Real actualThicknessChange = std::numeric_limits<Real>::quiet_NaN();
        Real predictedThicknessChange = std::numeric_limits<Real>::quiet_NaN();
        Real dThicknessTheta = std::numeric_limits<Real>::quiet_NaN();
        Real eikonalJump = std::numeric_limits<Real>::quiet_NaN();
        Real projectedJump = std::numeric_limits<Real>::quiet_NaN();
        Real distanceCorrection = std::numeric_limits<Real>::quiet_NaN();
        Real interfaceShift = std::numeric_limits<Real>::quiet_NaN();
        Real advectionIncrement = std::numeric_limits<Real>::quiet_NaN();
        Real advectedJump = std::numeric_limits<Real>::quiet_NaN();
        Real transportAttempted = std::numeric_limits<Real>::quiet_NaN();
        Real transportOmitted = std::numeric_limits<Real>::quiet_NaN();
        Real transportOmittedFraction = std::numeric_limits<Real>::quiet_NaN();
        Real transportFallbackCells = std::numeric_limits<Real>::quiet_NaN();
    } stageDiagnostics;
    const auto writeHistory = [&](Real nullSpaceMultiplier, Real xiRhoInfinityNorm,
                                Real thetaInfinityNorm, Real dRhoTheta, Real dVolumeTheta,
                                Real requiredDVolumeTheta,
                                const GradientDiagnostics& rhoGradientDiagnostics,
                                const GradientDiagnostics& volumeGradientDiagnostics) {
      history << iteration << ',' << outerRadius << ',' << initialGridSpacing << ','
              << m_options.configuration.hminFactor << ','
              << m_options.configuration.hmaxFactor << ',' << m_options.configuration.hmin
              << ',' << m_options.configuration.hmax << ',' << h << ',' << dt << ','
              << m_options.advectionQuadratureOrder << ',' << hilbertLength << ','
              << nitschePenalty << ',' << stabilizationFactor << ',' << assemblyBackend
              << ',' << KelvinBall::DirectSolverName << ',' << meshDiagnostics.vertices
              << ',' << meshDiagnostics.cells << ',' << meshDiagnostics.obstacleCells
              << ',' << meshDiagnostics.fluidCells << ','
              << meshDiagnostics.interfaceTriangles << ','
              << meshDiagnostics.minimumQuality << ',' << meshDiagnostics.meanQuality
              << ',' << meshDiagnostics.maximumQuality << ',' << backgroundH << ','
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
      // Free motion under a unit force without torque, R [Z; w] = [F; 0].
      const Real determinant = k * q - c * c;
      history << ',' << q / determinant << ',' << -c / determinant << ','
              << k / determinant << ','
              << (c != 0 ? 2 * M_PI / motionAngular.norm() : nan) << ','
              << (c != 0 ? 2 * M_PI * q / std::abs(c) : nan) << ',' << std::sqrt(q / k)
              << ',' << m_options.levelSetPenalty << ',' << m_options.minimumThickness
              << ',' << normalLength << ',' << stageDiagnostics.thicknessPenalty << ','
              << stageDiagnostics.thicknessNominalPenalty << ','
              << stageDiagnostics.thicknessGuardDistance << ','
              << stageDiagnostics.thicknessActiveMultiplier << ','
              << stageDiagnostics.thicknessActiveFeasible << ','
              << stageDiagnostics.thicknessViolating << ','
              << stageDiagnostics.thicknessMaximumPairLoad << ','
              << stageDiagnostics.thicknessTopOnePercentLoadShare << ','
              << stageDiagnostics.thicknessDescentNormalRMS << ','
              << stageDiagnostics.thicknessDescentTangentRMS << ','
              << stageDiagnostics.thicknessDescentTangentFraction << ','
              << stageDiagnostics.thetaTangentFraction << ','
              << stageDiagnostics.thicknessMinimumExit << ','
              << stageDiagnostics.thicknessMaximumDeficit << ','
              << stageDiagnostics.thicknessMinimumTransversality << ','
              << stageDiagnostics.thicknessNormalAlignmentMin << ','
              << stageDiagnostics.thicknessNormalAlignmentMean << ','
              << stageDiagnostics.thicknessCurvatureMin << ','
              << stageDiagnostics.thicknessCurvatureMax << ','
              << stageDiagnostics.thicknessNormalMagnitudeMin << ','
              << stageDiagnostics.thicknessNormalMagnitudeMax << ','
              << stageDiagnostics.thicknessGeometricSeamJump << ','
              << stageDiagnostics.thicknessNormalSeamJump << ','
              << stageDiagnostics.thicknessRaySeamJump << ','
              << stageDiagnostics.thicknessSeamSamples << ','
              << stageDiagnostics.thicknessSeamUnmatched << ','
              << stageDiagnostics.actualThicknessChange << ','
              << stageDiagnostics.predictedThicknessChange << ','
              << stageDiagnostics.dThicknessTheta << ',' << stageDiagnostics.eikonalJump
              << ',' << stageDiagnostics.projectedJump << ','
              << stageDiagnostics.distanceCorrection << ','
              << stageDiagnostics.interfaceShift << ','
              << stageDiagnostics.advectionIncrement << ','
              << stageDiagnostics.advectedJump << ',' << reconstruction.minimumCrossing
              << ',' << reconstruction.snappedVertices << ',' << reconstruction.scale;
      // Z and omega for the supplied force, including its magnitude.
      for (Eigen::Index component = 0; component < 3; ++component)
        history << ',' << motionTranslation(component);
      for (Eigen::Index component = 0; component < 3; ++component)
        history << ',' << motionAngular(component);
      history << ',' << stageDiagnostics.transportAttempted << ','
              << stageDiagnostics.transportOmitted << ','
              << stageDiagnostics.transportOmittedFraction << ','
              << stageDiagnostics.transportFallbackCells;
      for (Eigen::Index component = 0; component < 3; ++component)
        history << ',' << m_options.motionForce(component);
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
                  << vectorText(chamberBoundingBox) << Alert::NewLine
                  << diagnosticLabel("Sewn mesh bounding box:")
                  << vectorText(sewnBoundingBox) << Alert::NewLine
                  << diagnosticLabel("Nitsche jump:")
                  << Alert::Notation::Number(nitscheJump) << Alert::NewLine
                  << diagnosticLabel("Coupling symmetry residual:")
                  << Alert::Notation::Number(couplingSymmetry) << Alert::Raise;
    {
      Alert::Info() << substageHeading("Free motion under the prescribed force")
                    << Alert::NewLine << diagnosticLabel("Force:")
                    << vectorText(m_options.motionForce) << Alert::NewLine
                    << diagnosticLabel("Translational velocity Z:")
                    << vectorText(motionTranslation) << Alert::NewLine
                    << diagnosticLabel("Angular velocity omega:")
                    << vectorText(motionAngular) << Alert::NewLine
                    << diagnosticLabel("Period of one revolution:")
                    << Alert::Notation::Number(
                         c != 0 ? 2 * M_PI / motionAngular.norm() : nan)
                    << Alert::NewLine << diagnosticLabel("Pitch:")
                    << Alert::Notation::Number(c != 0 ? 2 * M_PI * q / std::abs(c) : nan)
                    << Alert::NewLine << diagnosticLabel("Resistance length sqrt(q/k):")
                    << Alert::Notation::Number(std::sqrt(q / k)) << Alert::Raise;
    }
    stageSeconds[3] = elapsedSeconds(stage4Start);
    reportStageTiming(4, stageSeconds[3]);
    if (m_options.stateOnly)
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
    // The paired first-exit load defines the local thickness correction.
    Math::Vector<Real> thicknessLoad;
    GridFunction thicknessDescent(shapeSpace);
    GridFunction smoothedNormal(shapeSpace);
    smoothedNormal.getData().setZero();
    GridFunction geometricNormal(shapeSpace);
    geometricNormal.getData().setZero();
    GridFunction rayDirection(shapeSpace);
    rayDirection.getData().setZero();
    GridFunction smoothedCurvature(levelSetSpace);
    smoothedCurvature.getData().setZero();
    if (m_options.minimumThickness > 0)
    {
      thicknessLoad = Math::Vector<Real>::Zero(shapeSpace.getSize());
      mesh.getConnectivity().compute(mesh.getDimension() - 1, mesh.getDimension());
      const Real activeGuard = m_options.minimumThickness + Real(2) * dt;
      const KelvinBall::ThicknessPenalty thicknessPenalty(mesh, activeGuard);
      Optional<GridFunction<decltype(shapeSpace), Math::Vector<Real>>> projectedNormal;
      if (m_options.thicknessDiagnostics)
      {
        projectedNormal.emplace(thicknessPenalty.projectNormal(shapeSpace, normalLength));
        smoothedNormal = *projectedNormal;
        for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        {
          const auto dofs = shapeSpace.getDOFs(0, vertex);
          Math::SpatialVector<Real> normal(3);
          for (size_t component = 0; component < 3; ++component)
            normal(component) = smoothedNormal.getData()(dofs(component));
          const Real magnitude = normal.norm();
          if (std::isfinite(magnitude) && magnitude > Real(1e-12))
            for (size_t component = 0; component < 3; ++component)
              smoothedNormal.getData()(dofs(component)) /= magnitude;
          else
            for (size_t component = 0; component < 3; ++component)
              smoothedNormal.getData()(dofs(component)) = 0;
        }
      }
      const auto thickness = thicknessPenalty.evaluate(mesh, shapeSpace, thicknessLoad,
        m_options.minimumThickness, projectedNormal ? &*projectedNormal : nullptr);
      stageDiagnostics.thicknessNominalPenalty = thickness.nominalPenalty;
      stageDiagnostics.thicknessGuardDistance = activeGuard;
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        const auto dofs = shapeSpace.getDOFs(0, vertex);
        for (size_t component = 0; component < 3; ++component)
        {
          rayDirection.getData()(dofs(component)) =
            thickness.rayDirection[vertex](component);
          geometricNormal.getData()(dofs(component)) =
            thickness.geometricNormal[vertex](component);
        }
      }
      const auto geometricSeam =
        interfaceSeamJump(mesh, geometricNormal, shapeCoupling.getLocator());
      const auto normalSeam =
        interfaceSeamJump(mesh, smoothedNormal, shapeCoupling.getLocator());
      const auto raySeam =
        interfaceSeamJump(mesh, rayDirection, shapeCoupling.getLocator());
      stageDiagnostics.thicknessGeometricSeamJump = geometricSeam.maximum;
      stageDiagnostics.thicknessNormalSeamJump =
        m_options.thicknessDiagnostics ? normalSeam.maximum : nan;
      stageDiagnostics.thicknessRaySeamJump = raySeam.maximum;
      stageDiagnostics.thicknessSeamSamples = static_cast<Real>(raySeam.samples);
      stageDiagnostics.thicknessSeamUnmatched = static_cast<Real>(raySeam.unmatched);
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        const auto dofs = levelSetSpace.getDOFs(0, vertex);
        smoothedCurvature.getData()(dofs(0)) = thickness.curvature[vertex];
      }
      stageDiagnostics.thicknessPenalty = thickness.penalty;
      stageDiagnostics.thicknessViolating = static_cast<Real>(thickness.violating);
      stageDiagnostics.thicknessMaximumPairLoad = thickness.maximumPairLoad;
      stageDiagnostics.thicknessTopOnePercentLoadShare = thickness.topOnePercentLoadShare;
      stageDiagnostics.thicknessMinimumExit = thickness.minimumExit;
      stageDiagnostics.thicknessMaximumDeficit = thickness.maximumDeficit;
      stageDiagnostics.thicknessMinimumTransversality = thickness.minimumTransversality;
      stageDiagnostics.thicknessNormalAlignmentMin =
        m_options.thicknessDiagnostics ? thickness.minimumNormalAlignment : nan;
      stageDiagnostics.thicknessNormalAlignmentMean =
        m_options.thicknessDiagnostics ? thickness.meanNormalAlignment : nan;
      stageDiagnostics.thicknessCurvatureMin =
        m_options.thicknessDiagnostics ? thickness.minimumCurvature : nan;
      stageDiagnostics.thicknessCurvatureMax =
        m_options.thicknessDiagnostics ? thickness.maximumCurvature : nan;
      stageDiagnostics.thicknessNormalMagnitudeMin =
        m_options.thicknessDiagnostics ? thickness.minimumNormalMagnitude : nan;
      stageDiagnostics.thicknessNormalMagnitudeMax =
        m_options.thicknessDiagnostics ? thickness.maximumNormalMagnitude : nan;
      stageDiagnostics.actualThicknessChange =
        previousThicknessPenalty ? thickness.penalty - *previousThicknessPenalty : nan;
      stageDiagnostics.predictedThicknessChange =
        predictedThicknessChange ? *predictedThicknessChange : nan;
      Alert::Info() << substageHeading("Thickness penalty") << Alert::NewLine
                    << diagnosticLabel("Minimum thickness:")
                    << Alert::Notation::Number(m_options.minimumThickness)
                    << Alert::NewLine << diagnosticLabel("Guard distance:")
                    << Alert::Notation::Number(activeGuard) << Alert::NewLine
                    << diagnosticLabel("Normal smoothing length:")
                    << Alert::Notation::Number(normalLength) << Alert::NewLine
                    << diagnosticLabel("Separation load:")
                    << "Bounded paired-normal surrogate" << Alert::NewLine
                    << diagnosticLabel("Penalty P:")
                    << Alert::Notation::Number(thickness.penalty) << Alert::NewLine
                    << diagnosticLabel("Nominal penalty:")
                    << Alert::Notation::Number(stageDiagnostics.thicknessNominalPenalty)
                    << Alert::NewLine << diagnosticLabel("Rays exiting before minimum:")
                    << Alert::Notation::Number(thickness.violating) << " of "
                    << Alert::Notation::Number(thickness.samples) << Alert::NewLine
                    << diagnosticLabel("Largest estimated pair load:")
                    << Alert::Notation::Number(thickness.maximumPairLoad)
                    << Alert::NewLine << diagnosticLabel("Top 1% estimated load share:")
                    << Alert::Notation::Number(thickness.topOnePercentLoadShare)
                    << Alert::NewLine << diagnosticLabel("Minimum exit distance:")
                    << Alert::Notation::Number(thickness.minimumExit) << Alert::NewLine
                    << diagnosticLabel("Maximum thickness deficit:")
                    << Alert::Notation::Number(thickness.maximumDeficit) << Alert::NewLine
                    << diagnosticLabel("Minimum exit transversality:")
                    << Alert::Notation::Number(thickness.minimumTransversality)
                    << Alert::NewLine << diagnosticLabel("Geometric normal cut jump:")
                    << Alert::Notation::Number(geometricSeam.maximum) << Alert::NewLine
                    << diagnosticLabel("Ray direction cut jump:")
                    << Alert::Notation::Number(raySeam.maximum) << Alert::NewLine
                    << diagnosticLabel("Interface-cut samples:")
                    << Alert::Notation::Number(raySeam.samples) << Alert::NewLine
                    << diagnosticLabel("Unmatched cut samples:")
                    << Alert::Notation::Number(raySeam.unmatched) << Alert::Raise;
      if (m_options.thicknessDiagnostics)
        Alert::Info() << substageHeading("Thickness normal diagnostics") << Alert::NewLine
                      << diagnosticLabel("Minimum projected alignment:")
                      << Alert::Notation::Number(thickness.minimumNormalAlignment)
                      << Alert::NewLine << diagnosticLabel("Mean projected alignment:")
                      << Alert::Notation::Number(thickness.meanNormalAlignment)
                      << Alert::NewLine << diagnosticLabel("Minimum smoothed curvature:")
                      << Alert::Notation::Number(thickness.minimumCurvature)
                      << Alert::NewLine << diagnosticLabel("Maximum smoothed curvature:")
                      << Alert::Notation::Number(thickness.maximumCurvature)
                      << Alert::NewLine << diagnosticLabel("Raw normal magnitude range:")
                      << Alert::Notation::Number(thickness.minimumNormalMagnitude)
                      << " to "
                      << Alert::Notation::Number(thickness.maximumNormalMagnitude)
                      << Alert::NewLine << diagnosticLabel("Projected normal cut jump:")
                      << Alert::Notation::Number(normalSeam.maximum) << Alert::Raise;
      if (previousThicknessPenalty)
        Alert::Info() << substageHeading("Thickness change (surrogate prediction)")
                      << Alert::NewLine << diagnosticLabel("Actual delta P:")
                      << Alert::Notation::Number(stageDiagnostics.actualThicknessChange)
                      << Alert::NewLine << diagnosticLabel("Surrogate linear delta P:")
                      << Alert::Notation::Number(
                           stageDiagnostics.predictedThicknessChange)
                      << Alert::Raise;
    }
    GradientDiagnostics rhoGradientDiagnostics{}, volumeGradientDiagnostics{};
    {
    // Retain the factorization only while identifying the three differentials.
    // Release it before output, transport, and reconstruction allocate systems.
      HilbertIdentification identification(
        shapeSpace, shapeCoupling, hilbertLength, nitschePenalty);
      if (m_options.minimumThickness > 0)
      {
        identifyGradient(identification, RealFunction{0}, thicknessDescent,
          "Thickness descent", &thicknessLoad);
        const auto components = interfaceComponents(mesh, thicknessDescent);
        stageDiagnostics.thicknessDescentNormalRMS = components.normalRMS;
        stageDiagnostics.thicknessDescentTangentRMS = components.tangentialRMS;
        stageDiagnostics.thicknessDescentTangentFraction = components.tangentialFraction;
        Alert::Info() << substageHeading("Thickness interface direction")
                      << Alert::NewLine << diagnosticLabel("Normal RMS:")
                      << Alert::Notation::Number(components.normalRMS) << Alert::NewLine
                      << diagnosticLabel("Tangential RMS:")
                      << Alert::Notation::Number(components.tangentialRMS)
                      << Alert::NewLine << diagnosticLabel("Tangential fraction:")
                      << Alert::Notation::Number(components.tangentialFraction)
                      << Alert::Raise;
      }
      rhoGradientDiagnostics =
        identifyGradient(identification, -rhoDensity, rhoGradient, "Rho");
      volumeGradientDiagnostics =
        identifyGradient(identification, RealFunction{-1}, volumeGradient, "Volume");
    }
    TestFunction testField(shapeSpace);
    auto interfaceNormal = FaceNormal(mesh);
    interfaceNormal.traceOf(Fluid);
    LinearForm dRho(testField), dVolume(testField);
    dRho = FaceIntegral(-rhoDensity, Dot(interfaceNormal, testField)).over(Gamma);
    dVolume = FaceIntegral(RealFunction{-1}, Dot(interfaceNormal, testField)).over(Gamma);
    dRho.assemble();
    dVolume.assemble();
    const Real volumeMetric = dVolume(volumeGradient);
    if (!(std::isfinite(volumeMetric) && std::abs(volumeMetric) > Real(1e-20)))
      throw std::runtime_error("The volume gradient has zero constraint derivative.");
    Real nullSpaceMultiplier = -dVolume(rhoGradient) / volumeMetric;
    GridFunction xiRho(shapeSpace);
    xiRho = rhoGradient + nullSpaceMultiplier * volumeGradient;
    if (m_options.minimumThickness > 0)
    {
      const Real thicknessVolumeMultiplier = -dVolume(thicknessDescent) / volumeMetric;
      GridFunction thicknessRepair(shapeSpace);
      thicknessRepair = thicknessDescent + thicknessVolumeMultiplier * volumeGradient;
      GridFunction candidate(shapeSpace);
      GridFunction candidateMagnitude(levelSetSpace);
      const auto thicknessRate = [&](Real multiplier) {
        candidate = xiRho + multiplier * thicknessRepair;
        candidateMagnitude = Frobenius(candidate);
        const Real magnitude = std::max(candidateMagnitude.max(), Real(1e-30));
        return -thicknessLoad.dot(candidate.getData()) / magnitude;
      };
      const Real targetRate =
        -stageDiagnostics.thicknessPenalty / m_options.minimumThickness;
      Real multiplier = 0;
      bool feasible = true;
      if (thicknessRate(0) > targetRate)
      {
        candidateMagnitude = Frobenius(thicknessRepair);
        const Real repairMagnitude = std::max(candidateMagnitude.max(), Real(1e-30));
        const Real maximumRepairRate =
          -thicknessLoad.dot(thicknessRepair.getData()) / repairMagnitude;
        if (maximumRepairRate >= targetRate)
        {
          feasible = false;
          multiplier = std::numeric_limits<Real>::infinity();
        }
        else
        {
          Real lower = 0, upper = 1;
          while (thicknessRate(upper) > targetRate && upper < Real(1e6))
            upper *= Real(2);
          if (thicknessRate(upper) > targetRate)
          {
            feasible = false;
            multiplier = std::numeric_limits<Real>::infinity();
          }
          else
          {
            for (size_t search = 0; search < 40; ++search)
            {
              const Real middle = (lower + upper) / Real(2);
              if (thicknessRate(middle) > targetRate)
                lower = middle;
              else
                upper = middle;
            }
            multiplier = upper;
          }
        }
      }
      if (feasible)
      {
        xiRho = rhoGradient + nullSpaceMultiplier * volumeGradient +
          multiplier * thicknessRepair;
        nullSpaceMultiplier += multiplier * thicknessVolumeMultiplier;
      }
      else
      {
        xiRho = thicknessRepair;
        nullSpaceMultiplier = thicknessVolumeMultiplier;
      }
      stageDiagnostics.thicknessActiveMultiplier = multiplier;
      stageDiagnostics.thicknessActiveFeasible = feasible ? Real(1) : Real(0);
      Alert::Info() << substageHeading("Aggregate thickness correction") << Alert::NewLine
                    << diagnosticLabel("Target D guard P [theta]:")
                    << Alert::Notation::Number(targetRate) << Alert::NewLine
                    << diagnosticLabel("Selected multiplier:")
                    << (feasible ? std::to_string(multiplier) : "Repair-only fallback")
                    << Alert::NewLine
                    << diagnosticLabel("Target feasible in tested span:")
                    << (feasible ? "Yes" : "No") << Alert::Raise;
    }
    GridFunction xiRhoNorm(levelSetSpace);
    xiRhoNorm = Frobenius(xiRho);
    const Real xiRhoInfinityNorm = std::max(xiRhoNorm.max(), Real(1e-30));
    GridFunction theta(shapeSpace);
    theta = xiRho;
    theta /= xiRhoInfinityNorm;
    stageDiagnostics.thetaTangentFraction =
      interfaceComponents(mesh, theta).tangentialFraction;
    xiRhoNorm = Frobenius(theta);
    const Real thetaInfinityNorm = xiRhoNorm.max();
    const Real dRhoTheta = dRho(theta);
    const Real dVolumeTheta = dVolume(theta);
    if (m_options.minimumThickness > 0)
    {
      stageDiagnostics.dThicknessTheta = -thicknessLoad.dot(theta.getData());
      previousThicknessPenalty = stageDiagnostics.thicknessPenalty;
      predictedThicknessChange = dt * stageDiagnostics.dThicknessTheta;
    }
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
    if (m_options.minimumThickness > 0)
      Alert::Info() << substageHeading("Thickness directional derivative")
                    << Alert::NewLine << diagnosticLabel("Surrogate D thickness [theta]:")
                    << Alert::Notation::Number(stageDiagnostics.dThicknessTheta)
                    << Alert::Raise;
    previousRho = rho;
    predictedRhoChange = dt * dRhoTheta;
    previousVolume = volume;
    predictedVolumeChange = dt * dVolumeTheta;
    GridFunction distance(levelSetSpace);
    {
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
    // The projection only reconciles the traces on the cuts: Gamma is pinned,
    // so it cannot move the interface the transport starts from.
      distanceProjection = Integral(periodicDistance, distanceTest) +
        shapeCoupling.scalarTrace(
          periodicDistance, distanceTest, m_options.levelSetPenalty) -
        eikonalDistanceLoad + DirichletBC(periodicDistance, RealFunction(0)).on(Gamma);
      Solver::CG(distanceProjection).solve();
      distance = periodicDistance.getSolution();
      const auto& distanceSystem = distanceProjection.getLinearSystem();
      const Real distanceProjectionResidual =
        (distanceSystem.getOperator() * distanceSystem.getSolution() -
          distanceSystem.getVector())
          .norm() /
        std::max(distanceSystem.getVector().norm(), Real(1));
      const Real distanceProjectionCorrection =
        (distance.getData() - eikonalDistance.getData()).lpNorm<Eigen::Infinity>();
    // The Eikonal distance vanishes on Gamma, so the projected value at an
    // interface vertex is how far the projection moves the interface there.
      std::vector<char> onInterface(mesh.getVertexCount(), 0);
      for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
      {
        if (face->getAttribute() == Gamma)
          for (const Index vertex : face->getVertices())
            onInterface[vertex] = 1;
      }
      Real interfaceShiftMaximum = 0;
      Real interfaceShiftSquares = 0;
      size_t interfaceVertices = 0;
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        if (!onInterface[vertex])
          continue;
        const Real shift = std::abs(distance[vertex]);
        interfaceShiftMaximum = std::max(interfaceShiftMaximum, shift);
        interfaceShiftSquares += shift * shift;
        ++interfaceVertices;
      }
      const Real interfaceShiftRMS = interfaceVertices
        ? std::sqrt(interfaceShiftSquares / static_cast<Real>(interfaceVertices))
        : Real(0);
      stageDiagnostics.eikonalJump = shapeCoupling.scalarJump(eikonalDistance);
      stageDiagnostics.projectedJump = shapeCoupling.scalarJump(distance);
      stageDiagnostics.distanceCorrection = distanceProjectionCorrection;
      stageDiagnostics.interfaceShift = interfaceShiftMaximum;
      Alert::Info() << substageHeading("Distance trace projection") << Alert::NewLine
                    << diagnosticLabel("Linear residual:")
                    << Alert::Notation::Number(distanceProjectionResidual)
                    << Alert::NewLine << diagnosticLabel("Eikonal rotated jump:")
                    << Alert::Notation::Number(stageDiagnostics.eikonalJump)
                    << Alert::NewLine << diagnosticLabel("Projected rotated jump:")
                    << Alert::Notation::Number(stageDiagnostics.projectedJump)
                    << Alert::NewLine << diagnosticLabel("Infinity correction:")
                    << Alert::Notation::Number(distanceProjectionCorrection)
                    << Alert::NewLine << diagnosticLabel("Interface shift maximum:")
                    << Alert::Notation::Number(interfaceShiftMaximum) << " = "
                    << Alert::Notation::Number(interfaceShiftMaximum / h) << " h"
                    << Alert::NewLine << diagnosticLabel("Interface shift RMS:")
                    << Alert::Notation::Number(interfaceShiftRMS) << " = "
                    << Alert::Notation::Number(interfaceShiftRMS / h) << " h"
                    << Alert::Raise;
    }

    stageSeconds[5] = elapsedSeconds(stage6Start);
    reportStageTiming(6, stageSeconds[5]);
    MMG::Mesh& advectionMesh = swiftBackground ? *swiftBackground : mesh;
    P1 advectionLevelSetSpace(advectionMesh);
    GridFunction advectedDistance(advectionLevelSetSpace);
    // Only the transported scalar survives this scope. In particular, the
    // full sewn meshes and their fields must expire before reconstruction.
    {
      const auto stage7Start = Clock::now();
      announce("Stage 7: Writing the chamber and sewn fields.");
      chamber.clear();
      chamber.add("Distance", distance, IO::XDMF::Center::Node);
      chamber.add("Theta", theta, IO::XDMF::Center::Node);
      if (m_options.minimumThickness > 0)
        chamber.add("Thickness_Descent", thicknessDescent, IO::XDMF::Center::Node);
      SubMesh<Context::Local>::Builder interfaceBuilder;
      interfaceBuilder.initialize(mesh);
      for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        if (face->getAttribute() == Gamma)
          interfaceBuilder.include(mesh.getDimension() - 1, face->getIndex());
      SubMesh<Context::Local> interfaceMesh = interfaceBuilder.finalize();
      P1 interfaceVectorSpace(interfaceMesh, 3);
      P1 interfaceScalarSpace(interfaceMesh);
      GridFunction interfaceGeometricNormal(interfaceVectorSpace);
      GridFunction interfaceSmoothedNormal(interfaceVectorSpace);
      GridFunction interfaceRayDirection(interfaceVectorSpace);
      GridFunction interfaceThicknessDescent(interfaceVectorSpace);
      GridFunction interfaceCurvature(interfaceScalarSpace);
      const auto& interfaceParentVertices = interfaceMesh.getPolytopeMap(0).left;
      if (m_options.minimumThickness > 0)
        for (Index vertex = 0; vertex < interfaceMesh.getVertexCount(); ++vertex)
        {
          const auto parent = interfaceParentVertices[vertex];
          const auto source = shapeSpace.getDOFs(0, parent);
          const auto target = interfaceVectorSpace.getDOFs(0, vertex);
          for (size_t component = 0; component < 3; ++component)
          {
            interfaceGeometricNormal.getData()(target(component)) =
              geometricNormal.getData()(source(component));
            interfaceSmoothedNormal.getData()(target(component)) =
              smoothedNormal.getData()(source(component));
            interfaceRayDirection.getData()(target(component)) =
              rayDirection.getData()(source(component));
            interfaceThicknessDescent.getData()(target(component)) =
              thicknessDescent.getData()(source(component));
          }
          interfaceCurvature.getData()(interfaceScalarSpace.getDOFs(0, vertex)(0)) =
            smoothedCurvature.getData()(levelSetSpace.getDOFs(0, parent)(0));
        }
      auto interfaceOutput = xdmf.grid("Interface");
      interfaceOutput.clear();
      interfaceOutput.setMesh(interfaceMesh, IO::XDMF::MeshPolicy::Transient);
      if (m_options.minimumThickness > 0)
      {
        interfaceOutput.add(
          "Geometric_Normal", interfaceGeometricNormal, IO::XDMF::Center::Node);
        interfaceOutput.add(
          "Smoothed_Normal", interfaceSmoothedNormal, IO::XDMF::Center::Node);
        interfaceOutput.add(
          "Ray_Direction", interfaceRayDirection, IO::XDMF::Center::Node);
        interfaceOutput.add(
          "Thickness_Descent", interfaceThicknessDescent, IO::XDMF::Center::Node);
        interfaceOutput.add(
          "Smoothed_Curvature", interfaceCurvature, IO::XDMF::Center::Node);
      }
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
      GridFunction chamberMotion(Vh);
      chamberMotion.getData() = motionTranslation(0) * uT0.getData() +
        motionTranslation(1) * uT1.getData() + motionTranslation(2) * uT2.getData() +
        motionAngular(0) * uR0.getData() + motionAngular(1) * uR1.getData() +
        motionAngular(2) * uR2.getData();
      fluidState.add("Motion", chamberMotion, IO::XDMF::Center::Node);

      KelvinBall::SewedOutput sewedDesign(mesh, FlatSet<Attribute>{Gamma, Outer});
      P1 sewedDesignScalar(sewedDesign.getMesh());
      P1 sewedDesignVector(sewedDesign.getMesh(), 3);
      GridFunction sewedDistance(sewedDesignScalar);
      GridFunction sewedVelocity(sewedDesignVector);
      GridFunction sewedGeometricNormal(sewedDesignVector);
      GridFunction sewedNormal(sewedDesignVector);
      GridFunction sewedRayDirection(sewedDesignVector);
      GridFunction sewedThicknessDescent(sewedDesignVector);
      GridFunction sewedCurvature(sewedDesignScalar);
      sewedDesign.setScalar(sewedDistance, distance);
      sewedDesign.setVector(sewedVelocity, theta);
      sewedDesignOutput.clear();
      sewedDesignOutput.setMesh(sewedDesign.getMesh(), IO::XDMF::MeshPolicy::Transient);
      sewedDesignOutput.add("Distance", sewedDistance, IO::XDMF::Center::Node);
      sewedDesignOutput.add("Theta", sewedVelocity, IO::XDMF::Center::Node);
      if (m_options.minimumThickness > 0)
      {
        sewedDesign.setVector(sewedGeometricNormal, geometricNormal);
        sewedDesign.setVector(sewedNormal, smoothedNormal);
        sewedDesign.setVector(sewedRayDirection, rayDirection);
        sewedDesign.setVector(sewedThicknessDescent, thicknessDescent);
        sewedDesign.setScalar(sewedCurvature, smoothedCurvature);
        sewedDesignOutput.add(
          "Thickness_Descent", sewedThicknessDescent, IO::XDMF::Center::Node);
      }

      SubMesh<Context::Local>::Builder sewedInterfaceBuilder;
      sewedInterfaceBuilder.initialize(sewedDesign.getMesh());
      for (auto face =
             sewedDesign.getMesh().getPolytope(sewedDesign.getMesh().getDimension() - 1);
        face; ++face)
        if (face->getAttribute() == Gamma)
          sewedInterfaceBuilder.include(
            sewedDesign.getMesh().getDimension() - 1, face->getIndex());
      SubMesh<Context::Local> sewedInterfaceMesh = sewedInterfaceBuilder.finalize();
      P1 sewedInterfaceVectorSpace(sewedInterfaceMesh, 3);
      P1 sewedInterfaceScalarSpace(sewedInterfaceMesh);
      GridFunction bodyMotion(sewedInterfaceVectorSpace);
      bodyMotion =
        VectorFunction(static_cast<size_t>(3), [&](const Geometry::Point& point) {
          const Math::SpatialVector<Real> velocity =
            motionTranslation + motionAngular.cross(point.getPhysicalCoordinates());
          return velocity;
        });
      GridFunction sewedInterfaceGeometricNormal(sewedInterfaceVectorSpace);
      GridFunction sewedInterfaceNormal(sewedInterfaceVectorSpace);
      GridFunction sewedInterfaceRayDirection(sewedInterfaceVectorSpace);
      GridFunction sewedInterfaceThicknessDescent(sewedInterfaceVectorSpace);
      GridFunction sewedInterfaceCurvature(sewedInterfaceScalarSpace);
      const auto& sewedInterfaceVertices = sewedInterfaceMesh.getPolytopeMap(0).left;
      if (m_options.minimumThickness > 0)
        for (Index vertex = 0; vertex < sewedInterfaceMesh.getVertexCount(); ++vertex)
        {
          const auto parent = sewedInterfaceVertices[vertex];
          const auto source = sewedDesignVector.getDOFs(0, parent);
          const auto target = sewedInterfaceVectorSpace.getDOFs(0, vertex);
          for (size_t component = 0; component < 3; ++component)
          {
            sewedInterfaceGeometricNormal.getData()(target(component)) =
              sewedGeometricNormal.getData()(source(component));
            sewedInterfaceNormal.getData()(target(component)) =
              sewedNormal.getData()(source(component));
            sewedInterfaceRayDirection.getData()(target(component)) =
              sewedRayDirection.getData()(source(component));
            sewedInterfaceThicknessDescent.getData()(target(component)) =
              sewedThicknessDescent.getData()(source(component));
          }
          sewedInterfaceCurvature.getData()(sewedInterfaceScalarSpace.getDOFs(0, vertex)(
            0)) = sewedCurvature.getData()(sewedDesignScalar.getDOFs(0, parent)(0));
        }
      auto sewedInterfaceOutput = sewedXdmf.grid("Interface");
      sewedInterfaceOutput.clear();
      sewedInterfaceOutput.setMesh(sewedInterfaceMesh, IO::XDMF::MeshPolicy::Transient);
      if (m_options.minimumThickness > 0)
      {
        sewedInterfaceOutput.add(
          "Geometric_Normal", sewedInterfaceGeometricNormal, IO::XDMF::Center::Node);
        sewedInterfaceOutput.add(
          "Smoothed_Normal", sewedInterfaceNormal, IO::XDMF::Center::Node);
        sewedInterfaceOutput.add(
          "Ray_Direction", sewedInterfaceRayDirection, IO::XDMF::Center::Node);
        sewedInterfaceOutput.add(
          "Thickness_Descent", sewedInterfaceThicknessDescent, IO::XDMF::Center::Node);
        sewedInterfaceOutput.add(
          "Smoothed_Curvature", sewedInterfaceCurvature, IO::XDMF::Center::Node);
      }

      KelvinBall::SewedOutput sewedFluid(fluid, FlatSet<Attribute>{Gamma, Outer});
      VelocitySpace sewedVelocitySpace(sewedFluid.getMesh(), 3);
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
      GridFunction fluidMotion(sewedVelocitySpace);
      fluidMotion.getData() = motionTranslation(0) * sewedUT0.getData() +
        motionTranslation(1) * sewedUT1.getData() +
        motionTranslation(2) * sewedUT2.getData() +
        motionAngular(0) * sewedUR0.getData() + motionAngular(1) * sewedUR1.getData() +
        motionAngular(2) * sewedUR2.getData();
      sewedFluidOutput.add("Motion", fluidMotion, IO::XDMF::Center::Node);
      sewedInterfaceOutput.add("Motion", bodyMotion, IO::XDMF::Center::Node);

      if (iteration + 1 == m_options.maxIterations)
      {
        xdmf.write(static_cast<Real>(iteration)).flush();
        sewedXdmf.write(static_cast<Real>(iteration)).flush();
        stageSeconds[6] = elapsedSeconds(stage7Start);
        reportStageTiming(7, stageSeconds[6]);
        writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
          dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
          volumeGradientDiagnostics);
        chamber.clear();
        fluidState.clear();
        interfaceOutput.clear();
        sewedDesignOutput.clear();
        sewedFluidOutput.clear();
        sewedInterfaceOutput.clear();
        continue;
      }

      stageSeconds[6] = elapsedSeconds(stage7Start);
      reportStageTiming(7, stageSeconds[6]);
      const auto stage8Start = Clock::now();
      announce("Stage 8: Advecting the level set.");
      auto& advectionConnectivity = advectionMesh.getConnectivity();
      advectionConnectivity.discover(3, 2);
      advectionConnectivity.discover(3, 1);
      advectionConnectivity.restrict(1, 0);
      advectionConnectivity.restrict(2, 0);
      advectionConnectivity.restrict(2, 3);
      advectionConnectivity.discover(0, 0);
      P1 advectionShapeSpace(advectionMesh, 3);
      GridFunction advectionDistance(advectionLevelSetSpace);
      GridFunction advectionDirection(advectionShapeSpace);
      if (swiftBackground)
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
        advectionDirection =
          VectorFunction(static_cast<size_t>(3), [&](const Geometry::Point& point) {
            Math::SpatialVector<Real> value(3);
            value.setZero();
            const auto located = fittedLocator.locate(3, point.getPhysicalCoordinates());
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
        Alert::Info()
          << substageHeading("Fitted-to-background transfer") << Alert::NewLine
          << diagnosticLabel("Nodal-copy error (distance):")
          << Alert::Notation::Number((advectionDistance.getData() - distance.getData())
                 .lpNorm<Eigen::Infinity>())
          << Alert::NewLine << diagnosticLabel("Nodal-copy error (direction):")
          << Alert::Notation::Number(
               (advectionDirection.getData() - theta.getData()).lpNorm<Eigen::Infinity>())
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
      const auto characteristic = Flow(-dt, advectionDistance, advectionDirection,
        Math::RungeKutta::RK4{}, rotationalContinuation);
    // P1 mass products are quadratic. A lower-order sampling rule cannot
    // determine all local coefficients, even with complete coverage.
      const size_t transportQuadratureOrder =
        std::max(size_t(2), m_options.advectionQuadratureOrder);
      const KelvinBall::TransportProjection projection(
        advectionMesh, advectionDistance, characteristic, transportQuadratureOrder);
      const RealFunction retained(
        [&](const Geometry::Point& point) { return projection.mask(point); });
      const RealFunction transported(
        [&](const Geometry::Point& point) { return projection.value(point); });
      auto mass = Integral(retained * advected, test);
      mass.setOrder(transportQuadratureOrder);
      auto transportedDistance = Integral(transported, test);
      transportedDistance.setOrder(transportQuadratureOrder);
      transport = mass +
        advectionCoupling.scalarTrace(advected, test, m_options.levelSetPenalty) -
        transportedDistance;
      // Project the transported field toward a zero periodic jump. Adding the
      // old jump to the load would preserve an existing mismatch instead.
      Solver::CG(transport).solve();
      stageDiagnostics.transportAttempted = projection.getAttemptedCount();
      stageDiagnostics.transportOmitted = projection.getOmittedCount();
      stageDiagnostics.transportOmittedFraction = projection.getOmittedWeightFraction();
      stageDiagnostics.transportFallbackCells = projection.getFallbackCellCount();
      Alert::Info() << substageHeading("Transport coverage") << Alert::NewLine
                    << diagnosticLabel("Projection quadrature order:")
                    << transportQuadratureOrder << Alert::NewLine
                    << diagnosticLabel("Attempted quadrature points:")
                    << projection.getAttemptedCount() << Alert::NewLine
                    << diagnosticLabel("Omitted quadrature points:")
                    << projection.getOmittedCount() << Alert::NewLine
                    << diagnosticLabel("Omitted integration fraction:")
                    << projection.getOmittedWeightFraction() << Alert::NewLine
                    << diagnosticLabel("Previous-distance fallback cells:")
                    << projection.getFallbackCellCount() << Alert::Raise;
      advectedDistance = advected.getSolution();
      const auto& transportSystem = transport.getLinearSystem();
      const Real transportResidual =
        (transportSystem.getOperator() * transportSystem.getSolution() -
          transportSystem.getVector())
          .norm() /
        std::max(transportSystem.getVector().norm(), Real(1));
      const Real advectionIncrement =
        (advectedDistance.getData() - advectionDistance.getData())
          .lpNorm<Eigen::Infinity>();
      const Real distanceJump = advectionCoupling.scalarJump(advectionDistance);
      const Real advectedJump = advectionCoupling.scalarJump(advectedDistance);
      stageDiagnostics.advectionIncrement = advectionIncrement;
      stageDiagnostics.advectedJump = advectedJump;
      Alert::Info() << substageHeading("Advected distance") << Alert::NewLine
                    << diagnosticLabel("Advection step:") << Alert::Notation::Number(dt)
                    << Alert::NewLine << diagnosticLabel("Linear residual:")
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
    // Attribute writers capture fields by reference. Remove those callbacks
    // before their fields expire; the written snapshot metadata is retained.
      chamber.clear();
      fluidState.clear();
      interfaceOutput.clear();
      sewedDesignOutput.clear();
      sewedFluidOutput.clear();
      sewedInterfaceOutput.clear();
    }

    const auto stage9Start = Clock::now();
    announce(m_options.reconstructionMethod == "swift"
        ? "Stage 9: Fitting the advected interface with SWIFT."
        : "Stage 9: Reconstructing the advected interface with MMG.");
    MMGReconstruction result = [&]() {
      if (m_options.reconstructionMethod == "swift")
      {
        MMG::Mesh classified =
          classifyLevelSetForSWIFT(*swiftBackground, advectedDistance);
        P1 classifiedLevelSetSpace(classified);
        GridFunction classifiedLevelSet(classifiedLevelSetSpace);
        classifiedLevelSet.getData() = advectedDistance.getData();
        auto fitted = fitLevelSetSWIFT(classified, classifiedLevelSet, backgroundH,
          initialGridSpacing, outerRadius, argc, argv);
        if (m_options.configuration.adapt)
          adaptSWIFT(fitted, sphere, m_options.configuration, requestedWelschScale);
        return fitted;
      }
      // A failed MMG stage is retried on the same advected level set with
      // every MMG size computed from half the scale; the next iteration
      // starts again from the reference spacing.
      Real scale = initialGridSpacing;
      for (size_t attempt = 0;; ++attempt)
      {
        try
        {
          const Real retryFactor = scale / initialGridSpacing;
          auto reconstructed = discretizeLevelSetMMG(mesh, advectedDistance, scale,
            retryFactor * m_options.configuration.hmin,
            retryFactor * m_options.configuration.hmax, sphere,
            m_options.configuration.adapt, m_options.configuration.mmgSnap,
            requestedWelschScale);
          reconstructed.diagnostics.scale = retryFactor;
          return reconstructed;
        }
        catch (const std::exception& error)
        {
          // Keep what MMG received, so that the failure can be reproduced.
          const std::string failure = "mmg-failure-" + std::to_string(iteration + 1) +
            "-" + std::to_string(attempt);
          mesh.save(failure + ".mesh", IO::FileFormat::MEDIT);
          advectedDistance.save(failure + ".sol", IO::FileFormat::MEDIT);
          Alert::Warning() << "Saved the failed MMG input to " << failure << ".mesh and "
                           << failure << ".sol." << Alert::Raise;
          if (attempt == m_options.configuration.mmgRetries)
            throw;
          Alert::Warning() << "MMG reconstruction failed at scale "
                           << Alert::Notation::Number(scale / initialGridSpacing)
                           << " h0: " << error.what() << Alert::NewLine
                           << "Retrying at scale "
                           << Alert::Notation::Number(scale / (2 * initialGridSpacing))
                           << " h0." << Alert::Raise;
          scale /= 2;
        }
      }
    }();
    reconstruction = result.diagnostics;
    if (reconstructionXdmf)
    {
      reconstructionOutput->clear();
      reconstructionOutput->setMesh(result.mesh, IO::XDMF::MeshPolicy::Transient);
      reconstructionXdmf->write(static_cast<Real>(iteration + 1)).flush();
    }
    checkFixedGeometry(result.mesh, outerRadius);
    checkMaterials(result.mesh);
    nextMesh.emplace(std::move(result.mesh));
    stageSeconds[8] = elapsedSeconds(stage9Start);
    reportStageTiming(9, stageSeconds[8]);
    writeHistory(nullSpaceMultiplier, xiRhoInfinityNorm, thetaInfinityNorm, dRhoTheta,
      dVolumeTheta, requiredDVolumeTheta, rhoGradientDiagnostics,
      volumeGradientDiagnostics);
  }

  Alert::Success() << "Wrote KelvinBall.xdmf, KelvinBallSewed.xdmf, "
                   << (reconstructionXdmf ? "KelvinBallMMG.xdmf, " : "")
                   << "and kelvin-ball.csv" << Alert::Raise;
  return 0;
}

int KelvinBall::KelvinBallOptimization::run()
{
  return m_implementation->run();
}
