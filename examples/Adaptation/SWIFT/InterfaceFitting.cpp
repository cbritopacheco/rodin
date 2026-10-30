/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Adaptation.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/QF/QuadratureFormula.h>
#include <Rodin/Variational.h>

#include <cmath>
#include <cstdlib>
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
    //! [swift-p1-mesh]
    constexpr size_t PointsPerAxis = 16;
    constexpr Real h = Real(1) / (PointsPerAxis - 1);
    constexpr Attribute Inside = 1, Outside = 2, Interface = 10;
    auto mesh =
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {PointsPerAxis, PointsPerAxis});
    mesh.scale(h);
    for (size_t from = 0; from <= 2; ++from)
    {
      for (size_t to = 0; to <= 2; ++to)
        mesh.getConnectivity().compute(from, to);
    }
    //! [swift-p1-mesh]

    //! [swift-p1-target]
    constexpr Real Radius = Real(0.25);
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
    //! [swift-p1-target]

    //! [swift-p1-classify]
    // Dimensionless phase width and smoothing for this unit-square example.
    constexpr Real PhaseWidth = Real(1.25), Smoothing = Real(0.01);
    constexpr size_t ClassificationOrder = 2;
    MinSTCut classifier(mesh);
    decltype(classifier)::Parameters classification;
    classification.smoothing = [](const Polytope&) { return Smoothing * h; };
    classifier.setParameters(classification);
    const auto labels = classifier.classify([&](const Polytope& cell) {
      const auto& formula =
        QF::QuadratureFormula::get(ClassificationOrder, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(formula);
      Real moment = 0;
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& point = quadrature.getPoint(q);
        moment += formula.getWeight(q) * point.getDistortion() *
          std::tanh(phi(point) / (PhaseWidth * h));
      }
      return moment / cell.getMeasure();
    });
    for (const Index cell : labels.inside)
      mesh.setAttribute({2, cell}, Inside);
    for (const Index cell : labels.outside)
      mesh.setAttribute({2, cell}, Outside);
    for (const Index edge : labels.cut)
      mesh.setAttribute({1, edge}, Interface);
    //! [swift-p1-classify]
    if (labels.cut.empty() || labels.inside.empty() || labels.outside.empty())
      return EXIT_FAILURE;

    //! [swift-p1-problem]
    P1 vh(mesh, 2);
    TrialFunction u(vh);
    TestFunction v(vh);
    SWIFT::Problem fitting(u, v);
    SWIFT::Parameters parameters;
    parameters.model.h = h;
    parameters.linear.solver = SWIFT::Parameters::LinearSolver::SparseLU;
    fitting.setParameters(parameters).setInterfaceAttribute(Interface);
    const auto report = fitting.solve(phi, gradient);
    //! [swift-p1-problem]

    //! [swift-p1-check]
    std::cout << std::setprecision(12) << "exit=" << report.getReasonString()
              << " target_hit=" << report.geometricTargetReached
              << " quality_ok=" << report.qualityBudgetSatisfied
              << " D_inf=" << report.geometricSup
              << " target=" << report.geometricSupTarget << " min_j=" << report.minJ
              << " max_Q=" << report.maxQRel << " outer=" << report.iterations << '\n';
    if (!report.geometricTargetReached || !report.qualityBudgetSatisfied)
      return EXIT_FAILURE;
    //! [swift-p1-check]

    //! [swift-p1-move]
    LocalMesh moved(mesh);
    const auto& displacement = u.getSolution();
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Point point(*vertex, Polytope::Traits(Polytope::Type::Point).getVertex(0));
      moved.setVertexCoordinates(vertex->getIndex(),
        point.getPhysicalCoordinates() + displacement.getValue(point));
    }
    //! [swift-p1-move]

    //! [swift-p1-output]
    IO::XDMF output("InterfaceFitting");
    output.grid("background").setMesh(mesh).add("displacement", displacement);
    output.grid("moved").setMesh(moved);
    output.write().close();
    //! [swift-p1-output]
    return EXIT_SUCCESS;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
