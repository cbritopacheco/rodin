/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <Rodin/Adaptation.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  // Optional arguments: grid points per axis, then spatial dimension (2 or 3).
  const size_t n = argc > 1 ? std::stoul(argv[1]) : 16;
  const size_t dimension = argc > 2 ? std::stoul(argv[2]) : 2;
  if (n < 4 || (dimension != 2 && dimension != 3))
  {
    std::cerr << "Expected at least 4 grid points and dimension 2 or 3.\n";
    return 1;
  }
  constexpr Attribute Inside = 1, Outside = 2, Interface = 10, Boundary = 20;
  constexpr Real Radius = Real(0.25);
  const Real h = Real(1) / Real(n - 1);
  auto mesh = dimension == 2
    ? Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {n, n})
    : Mesh<Context::Local>::UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
  mesh.scale(h);
  for (size_t from = 0; from <= dimension; ++from)
  {
    for (size_t to = 0; to <= dimension; ++to)
      mesh.getConnectivity().compute(from, to);
  }

  // A signed distance and its analytic gradient define a circle or sphere.
  Math::SpatialVector<Real> center(dimension);
  center.setConstant(Real(0.5));
  RealFunction phi([&](const Point& point) -> Real {
    return (point.getPhysicalCoordinates() - center).norm() - Radius;
  });
  VectorFunction gradient(dimension, [&](const Point& point) {
    const Math::SpatialVector<Real> radial(point.getPhysicalCoordinates() - center);
    // The distance gradient is undefined at the center; it is never on the target.
    return radial.norm() > Real(0) ? Math::SpatialVector<Real>(radial / radial.norm())
                                   : Math::SpatialVector<Real>::Zero(dimension);
  });

  // Classify cells by their centroids, then mark the separating mesh facets.
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    const Point centroid(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
    mesh.setAttribute(
      {dimension, cell->getIndex()}, phi(centroid) < Real(0) ? Inside : Outside);
  }
  for (auto face = mesh.getFace(); face; ++face)
  {
    const auto& cells =
      mesh.getConnectivity().getIncidence({dimension - 1, dimension}, face->getIndex());
    if (cells.size() == 1)
      mesh.setAttribute({dimension - 1, face->getIndex()}, Boundary);
    else if (mesh.getAttribute(dimension, cells[0]) !=
      mesh.getAttribute(dimension, cells[1]))
      mesh.setAttribute({dimension - 1, face->getIndex()}, Interface);
  }

  // Adapt owns the P1 space, solves the fit and applies a valid displacement.
  Adaptation::SWIFT::Adapt fitting(mesh);
  Adaptation::SWIFT::Parameters parameters;
  parameters.model.h = h; // The reference spacing stays fixed throughout fitting.
  fitting.setParameters(parameters).setInterfaceAttribute(Interface);
  const Math::Vector<Real> zero = Math::Vector<Real>::Zero(dimension);
  auto& u = fitting.getTrialFunction();
  fitting.getProblem() += DirichletBC(u, VectorFunction(zero)).on(Boundary);
  const auto report = fitting.execute(phi, gradient);

  // A valid best-effort fit may miss the target: report fit and quality separately.
  std::cout << std::scientific << std::setprecision(6)
            << "exit=" << report.getReasonString() << " energy=" << report.energy
            << " D_inf=" << report.geometricSup << " C=" << report.geometricConstant
            << " outer=" << report.iterations << " inner=" << report.innerIterations
            << " min_j=" << report.minJ << " max_Q=" << report.maxQRel << '\n';
  if (!report.qualityBudgetSatisfied)
    return 1;

  std::filesystem::create_directories("swift");
  IO::XDMF output("swift/ReconstructionP1");
  output.setMesh(mesh).write().close();
  return 0;
}
