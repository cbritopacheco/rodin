/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [background-source]
#include <Rodin/Adaptation.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MMG.h>
#include <Rodin/QF/QuadratureFormula.h>
#include <Rodin/Variational.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace SWIFT = Rodin::Adaptation::SWIFT;

int main()
{
  try
  {
    constexpr size_t PointsPerAxis = 13;
    constexpr Real ReferenceSpacing = Real(1) / (PointsPerAxis - 1);
    constexpr Real Radius = Real(0.25);
    // Size-field choices in unit-square coordinates; these are tutorial policies.
    constexpr Real NearSize = Real(0.035), FarSize = Real(0.14);
    constexpr Real BandWidth = Real(0.10), Gradation = Real(1.3);
    constexpr Real GeometricTolerance = Real(0.0025);
    constexpr Attribute Inside = 1, Outside = 2, Interface = 10, Boundary = 20;

    //! [background-target]
    RealFunction phi([](const Point& point) {
      return std::hypot(point.x() - Real(0.5), point.y() - Real(0.5)) - Radius;
    });
    VectorFunction gradient(2, [](const Point& point) {
      Math::SpatialVector<Real> radial(2);
      radial(0) = point.x() - Real(0.5);
      radial(1) = point.y() - Real(0.5);
      const Real length = radial.norm();
      if (length == 0)
        return Math::SpatialVector<Real>::Zero(2);
      return Math::SpatialVector<Real>(radial / length);
    });
    //! [background-target]

    MMG::Mesh mesh;
    mesh =
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {PointsPerAxis, PointsPerAxis});
    mesh.scale(ReferenceSpacing);
    mesh.getConnectivity().compute(1, 2);
    for (const Index corner :
      {Index(0), Index(PointsPerAxis - 1), Index(PointsPerAxis * (PointsPerAxis - 1)),
        Index(PointsPerAxis * PointsPerAxis - 1)})
      mesh.setCorner(corner);

    const MMG::Mesh uniform(mesh);
    IO::XDMF output("BackgroundFitting");
    output.grid("uniform").setMesh(uniform);
    const size_t initialCells = mesh.getCellCount();

    //! [background-adapt]
    {
      P1 sizeSpace(mesh);
      MMG::RealGridFunction size(sizeSpace);
      const SWIFT::Loss band(BandWidth);
      size = [&](const Point& point) {
        return FarSize - (FarSize - NearSize) * band.getWeight(phi(point));
      };
      MMG::Adapt().setHMin(NearSize).setHMax(FarSize).setGradation(Gradation).adapt(
        mesh, size);
    }
    // All fields below are constructed on the new mesh.
    //! [background-adapt]
    for (size_t from = 0; from <= 2; ++from)
    {
      for (size_t to = 0; to <= 2; ++to)
        mesh.getConnectivity().compute(from, to);
    }

    //! [background-classify]
    constexpr Real PhaseWidth = Real(1.25), Smoothing = Real(0.01);
    constexpr size_t ClassificationOrder = 2;
    MinSTCut classifier(static_cast<const LocalMesh&>(mesh));
    decltype(classifier)::Parameters classification;
    classification.smoothing = [&](const Polytope& face) {
      const auto& cells = mesh.getConnectivity().getIncidence({1, 2}, face.getIndex());
      Real localSize = std::sqrt(mesh.getPolytope(2, cells[0])->getMeasure());
      for (const Index cell : cells)
      {
        localSize =
          std::min(localSize, std::sqrt(mesh.getPolytope(2, cell)->getMeasure()));
      }
      return Smoothing * localSize;
    };
    classifier.setParameters(classification);
    const auto labels = classifier.classify([&](const Polytope& cell) {
      const Real localSize = std::sqrt(cell.getMeasure());
      const auto& formula =
        QF::QuadratureFormula::get(ClassificationOrder, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(formula);
      Real moment = 0;
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& point = quadrature.getPoint(q);
        moment += formula.getWeight(q) * point.getDistortion() *
          std::tanh(phi(point) / (PhaseWidth * localSize));
      }
      return moment / cell.getMeasure();
    });
    for (const Index cell : labels.inside)
      mesh.setAttribute({2, cell}, Inside);
    for (const Index cell : labels.outside)
      mesh.setAttribute({2, cell}, Outside);
    for (const Index face : labels.cut)
      mesh.setAttribute({1, face}, Interface);
    for (auto face = mesh.getFace(); face; ++face)
    {
      if (mesh.getConnectivity().getIncidence({1, 2}, face->getIndex()).size() == 1)
        mesh.setAttribute({1, face->getIndex()}, Boundary);
    }
    //! [background-classify]
    if (labels.cut.empty() || labels.inside.empty() || labels.outside.empty())
      return EXIT_FAILURE;
    const MMG::Mesh classified(mesh);
    const size_t fittedCells = mesh.getCellCount(),
                 fittedVertices = mesh.getVertexCount();
    output.grid("classified").setMesh(classified);

    //! [background-fit]
    SWIFT::Adapt fitting(mesh);
    SWIFT::Parameters parameters;
    parameters.model.h = NearSize;
    parameters.linear.solver = SWIFT::Parameters::LinearSolver::SparseLU;
    parameters.convergence.tolerance.geometric = GeometricTolerance;
    fitting.setParameters(parameters).setInterfaceAttribute(Interface);
    fitting.getProblem() +=
      DirichletBC(fitting.getTrialFunction(), VectorFunction(Real(0), Real(0)))
        .on(Boundary);
    std::ofstream history("BackgroundFitting.csv");
    history << std::setprecision(17) << "iteration,energy,geometric_error,min_j,max_q\n";
    fitting.getProblem().setMonitor([&](const SWIFT::Report& report) {
      history << report.iterations << ',' << report.energy << ',' << report.geometricSup
              << ',' << report.minJ << ',' << report.maxQRel << '\n';
    });
    const auto report = fitting.execute(phi, gradient);
    //! [background-fit]

    Real boundaryMotion = 0, minimumSignedJacobian = 1, exactInterfaceError = 0;
    for (const Index edge : labels.cut)
    {
      const auto& vertices = mesh.getPolytope(1, edge)->getVertices();
      Math::SpatialVector<Real> a(mesh.getVertexCoordinates(vertices[0]));
      Math::SpatialVector<Real> b(mesh.getVertexCoordinates(vertices[1]));
      a(0) -= Real(0.5);
      a(1) -= Real(0.5);
      b(0) -= Real(0.5);
      b(1) -= Real(0.5);
      const Math::SpatialVector<Real> direction(b - a);
      const Real t =
        std::clamp(-a.dot(direction) / direction.squaredNorm(), Real(0), Real(1));
      const Real nearestRadius = Math::SpatialVector<Real>(a + t * direction).norm();
      exactInterfaceError = std::max({exactInterfaceError, std::abs(a.norm() - Radius),
        std::abs(b.norm() - Radius), std::abs(nearestRadius - Radius)});
    }
    for (auto face = classified.getFace(); face; ++face)
    {
      if (face->getAttribute() == Boundary)
      {
        for (const Index vertex : face->getVertices())
        {
          boundaryMotion = std::max(boundaryMotion,
            (mesh.getVertexCoordinates(vertex) - classified.getVertexCoordinates(vertex))
              .norm());
        }
      }
    }
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& vertices = cell->getVertices();
      const Math::SpatialVector<Real> a(
        mesh.getVertexCoordinates(vertices[1]) - mesh.getVertexCoordinates(vertices[0]));
      const Math::SpatialVector<Real> b(
        mesh.getVertexCoordinates(vertices[2]) - mesh.getVertexCoordinates(vertices[0]));
      minimumSignedJacobian = std::min(minimumSignedJacobian, a(0) * b(1) - a(1) * b(0));
    }
    output.grid("fitted").setMesh(mesh);
    output.write().close();
    history.flush();
    std::cout << std::setprecision(17) << "uniform_cells=" << initialCells
              << " adapted_cells=" << fittedCells
              << " interface_edges=" << labels.cut.size()
              << " exit=" << report.getReasonString() << " error=" << report.geometricSup
              << " min_j=" << report.minJ << " max_q=" << report.maxQRel
              << " boundary_motion=" << boundaryMotion
              << " exact_interface_error=" << exactInterfaceError
              << " min_signed_jacobian=" << minimumSignedJacobian << '\n';
    constexpr Real RoundoffTolerance = Real(1e-12);
    return report.qualityBudgetSatisfied && report.geometricTargetReached &&
        std::isfinite(report.geometricSup) && minimumSignedJacobian > 0 &&
        exactInterfaceError <= GeometricTolerance &&
        std::abs(mesh.getArea() - Real(1)) < RoundoffTolerance &&
        boundaryMotion < RoundoffTolerance && bool(history) &&
        mesh.getCellCount() == fittedCells && mesh.getVertexCount() == fittedVertices
      ? EXIT_SUCCESS
      : EXIT_FAILURE;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
//! [background-source]
