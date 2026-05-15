/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include <Rodin/Adaptation.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/IO.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TMOP;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace
{
  static constexpr Real Pi = 3.141592653589793238462643383279502884;
  static constexpr Attribute InterfaceAttribute = 30;
  static constexpr Attribute NegativeAttribute = 1;
  static constexpr Attribute PositiveAttribute = 2;

  struct MovingCircle
  {
    Real radius = 0.18;

    Math::SpatialPoint center(Real t) const
    {
      return {0.50 + 0.14 * std::cos(2 * Pi * t),
              0.50 + 0.11 * std::sin(2 * Pi * t)};
    }

    Real signedDistance(const Math::SpatialPoint& x, Real t) const
    {
      const auto c = center(t);
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      return std::sqrt(dx * dx + dy * dy) - radius;
    }
  };

  struct InterfaceResiduals
  {
    Real p1LInf = 0;
    Real analyticLInf = 0;
    Real analyticL2 = 0;
  };

  struct MeshQualitySummary
  {
    Real minSignedArea = std::numeric_limits<Real>::infinity();
    Real minAbsArea = std::numeric_limits<Real>::infinity();
    Real minQuality = std::numeric_limits<Real>::infinity();
    Index invertedCells = 0;
    Index degenerateCells = 0;
    Index poorQualityCells = 0;
  };

  Real signedTriangleArea(
      const Math::SpatialPoint& x0,
      const Math::SpatialPoint& x1,
      const Math::SpatialPoint& x2)
  {
    return Real(0.5)
      * ((x1[0] - x0[0]) * (x2[1] - x0[1])
       - (x1[1] - x0[1]) * (x2[0] - x0[0]));
  }

  Real triangleQuality(
      const Math::SpatialPoint& x0,
      const Math::SpatialPoint& x1,
      const Math::SpatialPoint& x2)
  {
    const Real area = std::abs(signedTriangleArea(x0, x1, x2));
    const Real l0 = (x1 - x0).norm();
    const Real l1 = (x2 - x1).norm();
    const Real l2 = (x0 - x2).norm();
    const Real denom = l0 * l0 + l1 * l1 + l2 * l2;
    if (denom <= Real(0))
      return Real(0);
    return Real(4) * std::sqrt(Real(3)) * area / denom;
  }

  MeshQualitySummary summarizeMeshQuality(
      const LocalMesh& mesh,
      Real areaTolerance = 1e-14,
      Real qualityTolerance = 1e-8)
  {
    MeshQualitySummary summary;
    const auto& conn = mesh.getConnectivity();

    for (Index cellIndex = 0; cellIndex < static_cast<Index>(mesh.getCellCount()); ++cellIndex)
    {
      if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
        continue;
      const auto& cell = conn.getPolytope(2, cellIndex);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));

      const Real signedArea = signedTriangleArea(x0, x1, x2);
      const Real absArea = std::abs(signedArea);
      const Real quality = triangleQuality(x0, x1, x2);

      summary.minSignedArea = std::min(summary.minSignedArea, signedArea);
      summary.minAbsArea = std::min(summary.minAbsArea, absArea);
      summary.minQuality = std::min(summary.minQuality, quality);
      if (signedArea <= Real(0))
        summary.invertedCells++;
      if (absArea <= areaTolerance)
        summary.degenerateCells++;
      if (quality <= qualityTolerance)
        summary.poorQualityCells++;
    }

    if (!std::isfinite(summary.minSignedArea))
    {
      summary.minSignedArea = 0;
      summary.minAbsArea = 0;
      summary.minQuality = 0;
    }
    return summary;
  }

  Real lastQuadraticRatio(const std::vector<Real>& residuals)
  {
    for (size_t i = residuals.size(); i-- > 1;)
    {
      const Real previous = residuals[i - 1];
      const Real current = residuals[i];
      if (std::isfinite(previous) && std::isfinite(current)
          && previous > Real(0))
        return current / (previous * previous);
    }
    return std::numeric_limits<Real>::quiet_NaN();
  }

  Real maxResidualGrowth(const std::vector<Real>& residuals)
  {
    Real growth = 0;
    for (size_t i = 1; i < residuals.size(); ++i)
    {
      if (std::isfinite(residuals[i]) && std::isfinite(residuals[i - 1])
          && residuals[i - 1] > Real(0))
        growth = std::max(growth, residuals[i] / residuals[i - 1]);
    }
    return growth;
  }

  template <class FES, class Data>
  InterfaceResiduals measureInterfaceResiduals(
      const LocalMesh& cutMesh,
      const GridFunction<FES, Data>& phi,
      const InterfaceGraph& graph,
      const MovingCircle& circle,
      Real time)
  {
    InterfaceResiduals residuals;
    const auto& conn = cutMesh.getConnectivity();

    for (const auto& vertex : graph.vertices)
    {
      Real p1Value = 0;
      if (vertex.kind == InterfaceVertexKind::OriginalVertex)
      {
        assert(vertex.originalVertex);
        p1Value = phi[*vertex.originalVertex];
      }
      else
      {
        assert(vertex.parentEdge);
        assert(vertex.t);
        const auto& ev = conn.getPolytope(1, *vertex.parentEdge);
        p1Value = (1 - *vertex.t) * phi[ev[0]] + *vertex.t * phi[ev[1]];
      }

      const Real analytic = circle.signedDistance(vertex.x, time);
      residuals.p1LInf = std::max(residuals.p1LInf, std::abs(p1Value));
      residuals.analyticLInf = std::max(residuals.analyticLInf, std::abs(analytic));
      residuals.analyticL2 += analytic * analytic;
    }

    if (!graph.vertices.empty())
      residuals.analyticL2 =
        std::sqrt(residuals.analyticL2 / static_cast<Real>(graph.vertices.size()));
    return residuals;
  }

  template <class FES, class Data>
  void applyDisplacement(
      LocalMesh& mesh,
      const GridFunction<FES, Data>& displacement,
      Real scale = 1)
  {
    const auto& data = displacement.getData();
    const Index vertexCount = static_cast<Index>(mesh.getVertexCount());
    const size_t sdim = mesh.getSpaceDimension();
    for (Index vertex = 0; vertex < vertexCount; ++vertex)
    {
      auto x = mesh.getVertexCoordinates(vertex);
      for (size_t c = 0; c < sdim; ++c)
        x[c] += scale * data(vertex + static_cast<Index>(c) * vertexCount);
      mesh.setVertexCoordinates(vertex, x);
    }
  }

  // Analytic fit energy 1/2 * integral over the cut interface of phi^2, with
  // phi the manufactured signed distance. No discretizer report is needed:
  // the interface facets are identified purely by InterfaceAttribute.
  Real interfacePhiEnergy(
      const LocalMesh& mesh, const MovingCircle& circle, Real time)
  {
    const auto& conn = mesh.getConnectivity();
    Real energy = 0;
    for (Index e = 0; e < conn.getCount(1); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != InterfaceAttribute)
        continue;
      const auto& edge = conn.getPolytope(1, e);
      const auto x0 = mesh.getVertexCoordinates(edge(0));
      const auto x1 = mesh.getVertexCoordinates(edge(1));
      const Real length = (x1 - x0).norm();
      const Math::SpatialPoint mid = Real(0.5) * (x0 + x1);
      const Real phi = circle.signedDistance(mid, time);
      energy += length * phi * phi;
    }
    return Real(0.5) * energy;
  }
}

int main(int, char**)
{
  static constexpr size_t resolution = 20;
  static constexpr Index timeSteps = 100;
  static constexpr size_t maxCarryForwardCells = 6000;

  MovingCircle circle;

  IO::XDMF xdmf("LevelSetMovingCircle");

  std::ofstream diagnostics("LevelSetMovingCircle.csv");
  diagnostics << std::setprecision(17);
  diagnostics
    << "step,time,cx,cy,radius,nodal_linf,"
    << "interface_p1_linf,interface_analytic_linf,interface_analytic_l2,"
    << "background_vertices,background_cells,background_min_area,"
    << "background_min_quality,background_inverted_cells,"
    << "graph_vertices,graph_edges,chains,cut_vertices,cut_cells,"
    << "degenerate_cells,pathological_cuts,min_output_area,min_output_quality,"
    << "cut_min_area,cut_min_quality,cut_inverted_cells,cut_poor_quality_cells,"
    << "optimized_min_area,optimized_min_quality,optimized_inverted_cells,"
    << "optimized_degenerate_cells,optimized_poor_quality_cells,"
    << "linear_metric,optimized_linear_metric,fit_initial,fit_final,"
    << "source_fit_initial,source_fit_final,"
    << "geometry_step_scale,"
    << "jacobian_fd_relative_error,newton_max_residual_growth,"
    << "newton_last_quadratic_ratio,"
    << "newton_initial_residual,newton_final_residual,newton_step_norm,"
    << "newton_iterations,newton_converged,optimized_geometry_applied\n";

  SquaredDistanceMetric metric;

  // The uniform lattice is created ONCE at startup as the initial background.
  // Each step projects phi onto the current background, cuts it, TMOP-optimizes
  // the cut mesh, then that optimized cut mesh becomes the next step's
  // background (carry-forward). A clean triangle-only rebuild between steps
  // keeps it re-cuttable. Note: re-cutting an already-cut mesh inserts new
  // interface vertices every step, so the mesh grows over time.
  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, {resolution, resolution});
  background.scale(1.0 / static_cast<Real>(resolution - 1));

  for (Index step = 0; step < timeSteps; ++step)
  {
    const Real time =
      static_cast<Real>(step) / static_cast<Real>(timeSteps - 1);
    const auto c = circle.center(time);

    background.getConnectivity().compute(2, 1);
    background.getConnectivity().compute(1, 0);
    const auto backgroundQuality = summarizeMeshQuality(background);

    P1<Real, LocalMesh> phiSpace(background);
    GridFunction phi(phiSpace);

    // Project phi onto the current background for the current time.
    Real nodalLInf = 0;
    for (Index i = 0; i < background.getVertexCount(); ++i)
    {
      const auto x = background.getVertexCoordinates(i);
      const Real value = circle.signedDistance(x, time);
      phi[i] = value;
      nodalLInf = std::max(
          nodalLInf,
          std::abs(phi[i] - circle.signedDistance(x, time)));
    }

    auto result = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(InterfaceAttribute)
      .setNegativeCellAttribute(NegativeAttribute)
      .setPositiveCellAttribute(PositiveAttribute)
      .discretize();
    const auto cutQuality = summarizeMeshQuality(result.mesh);

    P1 displacementSpace(result.mesh, result.mesh.getSpaceDimension());
    GridFunction displacement(displacementSpace);
    displacement.getData().setZero();

    TrialFunction du(displacementSpace);
    TestFunction v(displacementSpace);

    static constexpr Real qualityWeight = 1.0;
    static constexpr Real deviationWeight = 1000000000000000.0;
    QualityTerm quality(qualityWeight);
    DeviationTerm deviation(deviationWeight);

    // Source-segment fit: this is the Rodin-native level-set fit term used by
    // the nonlinear TMOP problem. It uses cutter provenance to keep the fitted
    // interface near the P1 interface segments that produced it.
    static constexpr Real fitWeight = 100000000000000000.0;
    LevelSetFitTerm fit(result.interfaceGraph, result.report, fitWeight);
    LevelSetFitTerm fitDiagnostic(result.interfaceGraph, result.report, 1.0);

    auto makeTangent = [&]()
    {
      return
          quality.tangent(displacement, du, v)
        + deviation.tangent(du, v)
        + fit.tangent(displacement, du, v);
    };

    auto makeResidual = [&]()
    {
      return
          quality.residual(displacement, v)
        + deviation.residual(displacement, v)
        + fit.residual(displacement, v);
    };

    auto computeJacobianError = [&]()
    {
      BilinearForm tangentForm(du, v);
      tangentForm = makeTangent();
      tangentForm.assemble();

      LinearForm residualForm(v);
      residualForm = makeResidual();
      residualForm.assemble();
      const auto residual0 = residualForm.getVector();

      auto direction = displacement.getData();
      for (Eigen::Index i = 0; i < direction.size(); ++i)
        direction(i) = std::sin(static_cast<Real>(i + 1));
      const Real directionNorm = direction.norm();
      if (directionNorm <= Real(0))
        return std::numeric_limits<Real>::quiet_NaN();
      direction /= directionNorm;

      const auto original = displacement.getData();
      const Real eps = 1e-7;
      displacement.getData() = original + eps * direction;

      LinearForm residualPerturbed(v);
      residualPerturbed = makeResidual();
      residualPerturbed.assemble();
      const auto fd = (residualPerturbed.getVector() - residual0) / eps;
      const auto jd = tangentForm.getOperator() * direction;

      displacement.getData() = original;

      const Real denom = std::max<Real>({Real(1), fd.norm(), jd.norm()});
      return (fd - jd).norm() / denom;
    };

    const Real jacobianError = computeJacobianError();

    Variational::Problem tangentialProblem(du, v);
    // Rodin assembles LFIs with a minus sign in the RHS, so adding residual
    // terms here gives the Newton system J(u_k) du = -R(u_k).
    tangentialProblem =
        quality.tangent(displacement, du, v)
      + deviation.tangent(du, v)
      + fit.tangent(displacement, du, v)
      + quality.residual(displacement, v)
      + deviation.residual(displacement, v)
      + fit.residual(displacement, v);

    SparseLU linearSolver{tangentialProblem};
    NewtonSolver newton(linearSolver);
    newton
      .setMaxIterations(30)
      .setDampingFactor(1)
      .setAbsoluteTolerance(1e-10)
      .setRelativeTolerance(1e-8)
      .setStepTolerance(1e-10);

    std::vector<Real> newtonResiduals;
    newton.setMonitor([&](const auto& report)
    {
      newtonResiduals.push_back(report.final_residual);
    });

    const Real fitInitial = interfacePhiEnergy(result.mesh, circle, time);
    const Real sourceFitInitial =
      fitDiagnostic.sourceSegmentDistanceEnergy(result.mesh);
    newton.solve(displacement);
    const auto newtonReport = newton.getReport();
    const Real newtonGrowth = maxResidualGrowth(newtonResiduals);
    const Real newtonQuadraticRatio = lastQuadraticRatio(newtonResiduals);

    const auto residuals = measureInterfaceResiduals(
        background,
        phi,
        result.interfaceGraph,
        circle,
        time);

    const Real linearMetric = LinearMeshMetricObjective(metric).compute(result.mesh);
    Real geometryStepScale = 1;
    LocalMesh optimizedMesh = result.mesh;
    // Full Newton update only: no line search, damping, or post-solve scaling.
    applyDisplacement(optimizedMesh, displacement, geometryStepScale);
    optimizedMesh.getConnectivity().compute(2, 1);
    optimizedMesh.getConnectivity().compute(1, 0);
    const bool optimizedGeometryApplied = true;
    const Real optimizedLinearMetric =
      LinearMeshMetricObjective(metric).compute(optimizedMesh);

    // Rebuild the finite element space on the optimized cut mesh and reproject
    // phi onto it (exact resample of the manufactured signed distance), so the
    // field is visualized on the fitted geometry instead of the background.
    optimizedMesh.getConnectivity().compute(2, 1);
    optimizedMesh.getConnectivity().compute(1, 0);
    const auto optimizedQuality = summarizeMeshQuality(optimizedMesh);
    const Real fitFinal = interfacePhiEnergy(optimizedMesh, circle, time);
    const Real sourceFitFinal =
      fitDiagnostic.sourceSegmentDistanceEnergy(optimizedMesh);
    P1<Real, LocalMesh> optimizedPhiSpace(optimizedMesh);
    GridFunction optimizedPhi(optimizedPhiSpace);
    for (Index i = 0; i < optimizedMesh.getVertexCount(); ++i)
      optimizedPhi[i] =
        circle.signedDistance(optimizedMesh.getVertexCoordinates(i), time);

    auto fittedGrid = xdmf.grid("fitted");
    fittedGrid.setMesh(optimizedMesh, IO::XDMF::MeshPolicy::Transient);
    fittedGrid.clear();
    fittedGrid.add("phi", optimizedPhi, IO::XDMF::Center::Node);

    xdmf.write(time).flush();

    diagnostics
      << step << ','
      << time << ','
      << c[0] << ','
      << c[1] << ','
      << circle.radius << ','
      << nodalLInf << ','
      << residuals.p1LInf << ','
      << residuals.analyticLInf << ','
      << residuals.analyticL2 << ','
      << background.getVertexCount() << ','
      << background.getCellCount() << ','
      << backgroundQuality.minAbsArea << ','
      << backgroundQuality.minQuality << ','
      << backgroundQuality.invertedCells << ','
      << result.interfaceGraph.vertices.size() << ','
      << result.interfaceGraph.edges.size() << ','
      << result.interfaceGraph.chains.size() << ','
      << result.mesh.getVertexCount() << ','
      << result.mesh.getCellCount() << ','
      << result.report.degenerateCellCount << ','
      << result.report.pathologicalCutCount << ','
      << result.report.minOutputCellArea << ','
      << result.report.minOutputCellQuality << ','
      << cutQuality.minAbsArea << ','
      << cutQuality.minQuality << ','
      << cutQuality.invertedCells << ','
      << cutQuality.poorQualityCells << ','
      << optimizedQuality.minAbsArea << ','
      << optimizedQuality.minQuality << ','
      << optimizedQuality.invertedCells << ','
      << optimizedQuality.degenerateCells << ','
      << optimizedQuality.poorQualityCells << ','
      << linearMetric << ','
      << optimizedLinearMetric << ','
      << fitInitial << ','
      << fitFinal << ','
      << sourceFitInitial << ','
      << sourceFitFinal << ','
      << geometryStepScale << ','
      << jacobianError << ','
      << newtonGrowth << ','
      << newtonQuadraticRatio << ','
      << newtonReport.initial_residual << ','
      << newtonReport.final_residual << ','
      << newtonReport.final_step_norm << ','
      << newtonReport.iterations << ','
      << newtonReport.converged << ','
      << optimizedGeometryApplied << '\n';

    std::cout << "step " << step
              << " t=" << time
              << " graph=(" << result.interfaceGraph.vertices.size()
              << " vertices, " << result.interfaceGraph.edges.size()
              << " segments)"
              << " p1 residual=" << residuals.p1LInf
              << " analytic residual=" << residuals.analyticLInf
              << " tmop residual " << newtonReport.initial_residual
              << " -> " << newtonReport.final_residual
              << " analytic_fit " << fitInitial << " -> " << fitFinal
              << " source_fit " << sourceFitInitial << " -> " << sourceFitFinal
              << " step_scale=" << geometryStepScale
              << " qmin " << cutQuality.minQuality
              << " -> " << optimizedQuality.minQuality
              << " inverted " << cutQuality.invertedCells
              << " -> " << optimizedQuality.invertedCells
              << " jac_fd=" << jacobianError
              << " qratio=" << newtonQuadraticRatio
              << " optimized_geometry_applied=" << optimizedGeometryApplied
              << " bg_in=" << background.getCellCount()
              << " cut_cells=" << result.mesh.getCellCount()
              << std::endl;

    // Carry-forward: the optimized cut mesh becomes the next step's
    // background. Direct move, no rebuild workaround (Connectivity now copies
    // its index correctly, so the cut mesh round-trips through the cutter).
    background = std::move(optimizedMesh);
    if (background.getCellCount() > maxCarryForwardCells)
    {
      std::cout << "Stopping early: carry-forward mesh reached "
                << background.getCellCount() << " cells, above diagnostic cap "
                << maxCarryForwardCells << "." << std::endl;
      break;
    }
  }

  xdmf.close();

  std::cout << "Wrote LevelSetMovingCircle.xdmf and LevelSetMovingCircle.csv."
            << std::endl;
  std::cout << "The nodal error is zero by construction: phi is assigned from"
            << " the manufactured signed-distance function at mesh vertices."
            << std::endl;
  std::cout << "The interface analytic residual measures the expected P1"
            << " polygonal approximation error to the moving circle."
            << std::endl;
  std::cout << "The uniform lattice is created once at startup as the initial"
            << " background. Each step projects phi onto the current"
            << " background, cuts it, solves the active TMOP terms (quality +"
            << " deviation + source-segment level-set fit),"
            << " reprojects phi onto"
            << " the optimized cut mesh for export, and then carries that"
            << " optimized cut mesh forward as the next step's background. The"
            << " mesh grows over time since each re-cut inserts new interface"
            << " vertices; this diagnostic run stops once the carry-forward"
            << " mesh exceeds " << maxCarryForwardCells << " cells."
            << std::endl;
  return 0;
}
