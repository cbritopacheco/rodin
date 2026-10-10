/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <Rodin/Adaptation.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational.h>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <limits>

#include "Options.h"
#include "LobedSphereLevelSet.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  try
  {
    constexpr size_t Degree = RODIN_SWIFT_RECONSTRUCTION_DEGREE;
    const Examples::ReconstructionOptions options(argc, argv);
    if (options.help)
      return 0;
    const size_t n = options.n, dimension = options.dimension;
    constexpr Attribute Inside = 1, Outside = 2, Interface = 10, Boundary = 20;
    const Real h = options.parameters.model.h;
    auto mesh = dimension == 2
      ? Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {n, n})
      : Mesh<Context::Local>::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(h);
    for (size_t from = 0; from <= dimension; ++from)
    {
      for (size_t to = 0; to <= dimension; ++to)
        mesh.getConnectivity().compute(from, to);
    }

    // Targets and analytic gradients are independent of displacement degree.
    Math::SpatialVector<Real> center(dimension);
    center(0) = options.cx;
    center(1) = options.cy;
    if (dimension == 3)
      center(2) = options.cz;
    Examples::LobedSphereLevelSet sphere;
    if (dimension == 3)
      sphere.c = center;
    sphere.R0 = options.radius;
    sphere.amp = options.lobes == 0 ? Real(0) : options.amplitude;
    sphere.lobes = options.lobes;
    sphere.phase = options.phase;
    RealFunction phi([&](const Point& point) -> Real {
      if (dimension == 3)
        return sphere.phi(point.getPhysicalCoordinates());
      const Math::SpatialVector<Real> radial(point.getPhysicalCoordinates() - center);
      const Real angle = std::atan2(radial(1), radial(0));
      const Real perturbation = options.lobes == 0
        ? Real(0)
        : options.amplitude * std::cos(options.lobes * (angle - options.phase));
      return radial.norm() - options.radius - perturbation;
    });
    VectorFunction gradient(dimension, [&](const Point& point) {
      if (dimension == 3)
        return sphere.grad(point.getPhysicalCoordinates());
      const Math::SpatialVector<Real> radial(point.getPhysicalCoordinates() - center);
      // The distance gradient is undefined at the center; it is never on the target.
      const Real r = radial.norm();
      if (r == Real(0))
        return Math::SpatialVector<Real>::Zero(dimension);
      Math::SpatialVector<Real> value(radial / r);
      const Real angle = std::atan2(radial(1), radial(0));
      const Real angular = options.amplitude * options.lobes *
        std::sin(options.lobes * (angle - options.phase));
      value(0) -= angular * radial(1) / (r * r);
      value(1) += angular * radial(0) / (r * r);
      return value;
    });

    // Average a smooth phase indicator, then minimize the binary Potts energy.
    constexpr Real PhaseScaleOverH = Real(1.25);
    constexpr size_t ClassificationOrder = 2;
    MinSTCut classifier(mesh);
    decltype(classifier)::Parameters classificationParameters;
    classificationParameters.fidelity = options.classification.fidelity;
    const Real smoothing = h * options.classification.smoothing;
    classificationParameters.smoothing = [smoothing](const Polytope&) {
      return smoothing;
    };
    classifier.setParameters(classificationParameters);
    const auto classified = classifier.classify([&](const Polytope& cell) -> Real
    {
      const auto& formula =
        QF::PolytopeQuadratureFormula::get(ClassificationOrder, cell.getGeometry());
      const auto& quadrature = cell.getQuadrature(formula);
      Real moment = 0;
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        const auto& point = quadrature.getPoint(q);
        moment += formula.getWeight(q) * point.getDistortion() *
          std::tanh(phi(point) / (PhaseScaleOverH * h));
      }
      return moment / cell.getMeasure();
    });
    for (const Index cell : classified.inside)
      mesh.setAttribute({dimension, cell}, Inside);
    for (const Index cell : classified.outside)
      mesh.setAttribute({dimension, cell}, Outside);
    for (const Index face : classified.cut)
      mesh.setAttribute({dimension - 1, face}, Interface);
    for (auto face = mesh.getFace(); face; ++face)
    {
      const auto& cells =
        mesh.getConnectivity().getIncidence({dimension - 1, dimension}, face->getIndex());
      if (cells.size() == 1)
        mesh.setAttribute({dimension - 1, face->getIndex()}, Boundary);
    }

    Alert::Info() << "MinSTCut classification: " << classified.inside.size()
                  << " solid cells, " << classified.outside.size() << " fluid cells, "
                  << classified.cut.size() << " interface facets." << Alert::Raise;

    // The only degree-dependent choice: all spaces use the same SWIFT problem.
    H1 space(std::integral_constant<size_t, Degree>{}, mesh, dimension);
    TrialFunction u(space);
    TestFunction v(space);
    Adaptation::SWIFT::Problem fitting(u, v);
    fitting.setParameters(options.parameters).setInterfaceAttribute(Interface);
    const auto report = fitting.solve(phi, gradient);

    // A valid best-effort fit may miss the target: report fit and quality separately.
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "exit=" << report.getReasonString() << " energy=" << report.energy
              << " D_inf=" << report.geometricSup << " C=" << report.geometricConstant
              << " target=" << report.geometricSupTarget
              << " target_hit=" << report.geometricTargetReached
              << " quality_ok=" << report.qualityBudgetSatisfied
              << " outer=" << report.iterations << " inner=" << report.innerIterations
              << " inner_max=" << report.maxInnerIterations
              << " inner_last=" << report.lastInnerIterations
              << " inner_residual=" << report.innerResidual << " min_j=" << report.minJ
              << " max_Q=" << report.maxQRel << '\n';
    if (!report.qualityBudgetSatisfied)
      return 1;

    // Apply the solution to a separate mesh, including all higher-order geometry nodes.
    const auto& displacement = u.getSolution();
    Mesh<Context::Local> moved(mesh);
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const Point point(*vertex, Polytope::Traits(Polytope::Type::Point).getVertex(0));
      moved.setVertexCoordinates(vertex->getIndex(),
        point.getPhysicalCoordinates() + displacement.getValue(point));
    }
    for (size_t d = 1; d <= dimension; ++d)
    {
      for (auto entity = mesh.getPolytope(d); entity; ++entity)
      {
        RealH1Element<Degree> geometry(entity->getGeometry());
        PointCloud nodes(dimension, geometry.getCount());
        for (size_t node = 0; node < geometry.getCount(); ++node)
        {
          const Point point(*entity, geometry.getNode(node));
          const Math::SpatialPoint position(
            point.getPhysicalCoordinates() + displacement.getValue(point));
          for (size_t component = 0; component < dimension; ++component)
            nodes(component, node) = position(component);
        }
        moved.setPolytopeTransformation({d, entity->getIndex()},
          new ParametricTransformation<RealH1Element<Degree>>(
            std::move(nodes), geometry));
      }
    }
    const std::filesystem::path stem = options.output.empty()
      ? "swift/ReconstructionP" + std::to_string(Degree)
      : options.output;
    if (stem.has_parent_path())
      std::filesystem::create_directories(stem.parent_path());
    // Keep both configurations as named blocks with fields on their own meshes.
    H1 scalarSpace(std::integral_constant<size_t, Degree>{}, mesh);
    H1 movedScalarSpace(std::integral_constant<size_t, Degree>{}, moved);
    H1 movedVectorSpace(std::integral_constant<size_t, Degree>{}, moved, dimension);
    P0 cellSpace(mesh);
    P0 movedCellSpace(moved);
    GridFunction levelSet(scalarSpace);
    GridFunction movedLevelSet(movedScalarSpace);
    GridFunction movedDisplacement(movedVectorSpace);
    GridFunction classification(cellSpace);
    GridFunction movedClassification(movedCellSpace);
    GridFunction jacobian(cellSpace);
    GridFunction movedJacobian(movedCellSpace);
    GridFunction distortion(cellSpace);
    GridFunction movedDistortion(movedCellSpace);
    levelSet = phi;
    movedLevelSet = phi;
    classification = [&](const Point& point) -> Real {
      return *mesh.getAttribute(dimension, point.getPolytope().getIndex());
    };
    const auto displacementJacobian = Jacobian(displacement);
    jacobian = [&](const Point& point) -> Real {
      Adaptation::CellDeformation state(dimension);
      state.setDisplacementGradient(displacementJacobian.getValue(point));
      return state.getJacobian();
    };
    distortion = [&](const Point& point) -> Real {
      Adaptation::CellDeformation state(dimension);
      state.setDisplacementGradient(displacementJacobian.getValue(point));
      return state.getRelativeDistortion();
    };
    // The copied mesh retains topology and the identical finite-element DOF ordering.
    movedDisplacement.getData() = displacement.getData();
    movedClassification.getData() = classification.getData();
    movedJacobian.getData() = jacobian.getData();
    movedDistortion.getData() = distortion.getData();

    IO::XDMF output(stem.string());
    output.grid("background")
      .setMesh(mesh)
      .add("cell_label", classification, IO::XDMF::Center::Cell)
      .add("displacement", displacement)
      .add("phi", levelSet)
      .add("j", jacobian, IO::XDMF::Center::Cell)
      .add("q_rel", distortion, IO::XDMF::Center::Cell);
    output.grid("moved")
      .setMesh(moved)
      .add("cell_label", movedClassification, IO::XDMF::Center::Cell)
      .add("displacement", movedDisplacement)
      .add("phi", movedLevelSet)
      .add("j", movedJacobian, IO::XDMF::Center::Cell)
      .add("q_rel", movedDistortion, IO::XDMF::Center::Cell);
    output.write().close();
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
