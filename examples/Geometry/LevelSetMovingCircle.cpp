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

  template <class FES, class Data>
  InterfaceResiduals measureInterfaceResiduals(
      const LocalMesh& background,
      const GridFunction<FES, Data>& phi,
      const InterfaceGraph& graph,
      const MovingCircle& circle,
      Real time)
  {
    InterfaceResiduals residuals;
    const auto& conn = background.getConnectivity();

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
      const GridFunction<FES, Data>& displacement)
  {
    const auto& data = displacement.getData();
    const Index vertexCount = static_cast<Index>(mesh.getVertexCount());
    const size_t sdim = mesh.getSpaceDimension();
    for (Index vertex = 0; vertex < vertexCount; ++vertex)
    {
      auto x = mesh.getVertexCoordinates(vertex);
      for (size_t c = 0; c < sdim; ++c)
        x[c] += data(vertex + static_cast<Index>(c) * vertexCount);
      mesh.setVertexCoordinates(vertex, x);
    }
  }
}

int main(int, char**)
{
  static constexpr size_t resolution = 20;
  static constexpr Index timeSteps = 100;

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, {resolution, resolution});
  background.scale(1.0 / static_cast<Real>(resolution - 1));
  background.getConnectivity().compute(2, 1);
  background.getConnectivity().compute(1, 0);

  P1<Real, LocalMesh> phiSpace(background);
  GridFunction phi(phiSpace);
  MovingCircle circle;

  IO::XDMF xdmf("LevelSetMovingCircle");
  auto backgroundGrid = xdmf.grid("background");
  backgroundGrid.setMesh(background, IO::XDMF::MeshPolicy::Static);
  backgroundGrid.add("phi", phi, IO::XDMF::Center::Node);
  auto fittedGrid = xdmf.grid("fitted");

  std::ofstream diagnostics("LevelSetMovingCircle.csv");
  diagnostics << std::setprecision(17);
  diagnostics
    << "step,time,cx,cy,radius,nodal_linf,"
    << "interface_p1_linf,interface_analytic_linf,interface_analytic_l2,"
    << "graph_vertices,graph_edges,chains,cut_vertices,cut_cells,"
    << "degenerate_cells,pathological_cuts,min_output_area,min_output_quality,"
    << "linear_metric,optimized_linear_metric,fit_initial,fit_final,"
    << "newton_initial_residual,newton_final_residual,newton_step_norm,"
    << "newton_iterations,newton_converged,optimized_geometry_applied\n";

  SquaredDistanceMetric metric;

  for (Index step = 0; step < timeSteps; ++step)
  {
    const Real time =
      static_cast<Real>(step) / static_cast<Real>(timeSteps - 1);
    const auto c = circle.center(time);

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

    P1 displacementSpace(result.mesh, result.mesh.getSpaceDimension());
    GridFunction displacement(displacementSpace);
    displacement.getData().setZero();

    TrialFunction du(displacementSpace);
    TestFunction v(displacementSpace);

    QualityTerm quality(1.0);
    DeviationTerm deviation(1e-2);
    LevelSetFitTerm fit(result.interfaceGraph, result.report, 100.0);

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
      .setMaxIterations(3)
      .setDampingFactor(0.01)
      .setAbsoluteTolerance(1e-10)
      .setRelativeTolerance(1e-8)
      .setStepTolerance(1e-10);

    const Real fitInitial = fit.sourceSegmentDistanceEnergy(result.mesh);
    newton.solve(displacement);
    const auto newtonReport = newton.getReport();

    const auto residuals = measureInterfaceResiduals(
        background,
        phi,
        result.interfaceGraph,
        circle,
        time);

    const Real linearMetric = LinearMeshMetricObjective(metric).compute(result.mesh);
    LocalMesh optimizedMesh = result.mesh;
    applyDisplacement(optimizedMesh, displacement);
    const bool optimizedGeometryApplied = true;
    const Real optimizedLinearMetric =
      LinearMeshMetricObjective(metric).compute(optimizedMesh);
    const Real fitFinal = fit.sourceSegmentDistanceEnergy(optimizedMesh);

    fittedGrid.setMesh(optimizedMesh, IO::XDMF::MeshPolicy::Transient);
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
      << result.interfaceGraph.vertices.size() << ','
      << result.interfaceGraph.edges.size() << ','
      << result.interfaceGraph.chains.size() << ','
      << result.mesh.getVertexCount() << ','
      << result.mesh.getCellCount() << ','
      << result.report.degenerateCellCount << ','
      << result.report.pathologicalCutCount << ','
      << result.report.minOutputCellArea << ','
      << result.report.minOutputCellQuality << ','
      << linearMetric << ','
      << optimizedLinearMetric << ','
      << fitInitial << ','
      << fitFinal << ','
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
              << " fit " << fitInitial << " -> " << fitFinal
              << " optimized_geometry_applied=" << optimizedGeometryApplied
              << std::endl;
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
  std::cout << "Each step cuts a fresh fitted mesh, solves a Rodin Newton"
            << " tangential TMOP problem on that fitted topology, and exports"
            << " the displaced fitted mesh."
            << std::endl;
  return 0;
}
