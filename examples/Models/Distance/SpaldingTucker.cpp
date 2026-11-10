/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [distance-source]
#include <Rodin/Distance/SpaldingTucker.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MMG.h>
#include <Rodin/Solver.h>
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

int main()
{
  try
  {
    constexpr size_t PointsPerAxis = 32;
    constexpr Real MaximumSize = Real(0.02);
    constexpr Attribute BaseMaterial = 1, Inside = 3, Outside = 2, Interface = 5;
    //! [distance-target]
    RealFunction exact([](const Point& point) {
      return std::min(
        {std::hypot(point.x() - Real(0.75), point.y() - Real(0.25)) - Real(0.05),
          std::hypot(point.x() - Real(0.25), point.y() - Real(0.25)) - Real(0.25),
          std::hypot(point.x() - Real(0.75), point.y() - Real(0.75)) - Real(0.10)});
    });
    //! [distance-target]
    MMG::Mesh mesh;
    mesh =
      LocalMesh::UniformGrid(Polytope::Type::Triangle, {PointsPerAxis, PointsPerAxis});
    mesh.scale(Real(1) / (PointsPerAxis - 1));
    for (auto cell = mesh.getCell(); cell; ++cell)
      mesh.setAttribute({2, cell->getIndex()}, BaseMaterial);
    //! [distance-remesh]
    {
      P1 initialSpace(mesh);
      MMG::RealGridFunction levelSet(initialSpace);
      levelSet = exact;
      mesh = MMG::LevelSetDiscretizer()
               .split(BaseMaterial, {Inside, Outside})
               .setBoundaryReference(Interface)
               .setAngleDetection()
               .setHMax(MaximumSize)
               .discretize(levelSet);
    }
    mesh.getConnectivity().compute(1, 2);
    //! [distance-remesh]

    //! [distance-potential]
    P1 space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    Problem potential(u, v);
    potential = Integral(Grad(u), Grad(v)) + Integral(2 * v).over(Inside) - Integral(v) +
      DirichletBC(u, Zero()).on(Interface);
    Solver::CG solver(potential);
    constexpr Real LinearTolerance = Real(1e-10);
    solver.setTolerance(LinearTolerance).solve();
    if (!solver.success())
      return EXIT_FAILURE;
    const auto& system = potential.getLinearSystem();
    const Real relativeResidual =
      (system.getOperator() * system.getSolution() - system.getVector()).norm() /
      std::max(Real(1), system.getVector().norm());
    constexpr Real ResidualTolerance = Real(1e-9);
    if (!std::isfinite(relativeResidual) || relativeResidual > ResidualTolerance)
      return EXIT_FAILURE;
    const auto& poisson = u.getSolution();
    auto normalized = Distance::SpaldingTucker()(poisson);
    GridFunction reference(space);
    reference = exact;
    //! [distance-potential]

    Real maximumError = 0, meanError = 0, interfaceMaximum = 0;
    bool valid = mesh.getArea(Inside) > 0 && mesh.getArea(Outside) > 0;
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const Index dof = space.getDOFs(0, vertex)[0];
      const Real value = normalized.getData()(dof);
      const Real error = std::abs(value - reference.getData()(dof));
      valid = valid && std::isfinite(value) && std::isfinite(poisson.getData()(dof));
      maximumError = std::max(maximumError, error);
      meanError += error;
    }
    meanError /= mesh.getVertexCount();
    size_t interfaceCount = 0;
    for (auto face = mesh.getFace(); face; ++face)
    {
      if (face->getAttribute() == Interface)
      {
        ++interfaceCount;
        for (const Index vertex : face->getVertices())
        {
          const Index dof = space.getDOFs(0, vertex)[0];
          interfaceMaximum =
            std::max(interfaceMaximum, std::abs(normalized.getData()(dof)));
        }
      }
    }
    constexpr Real BoundaryTolerance = Real(1e-12);
    valid = valid && interfaceCount > 0 && interfaceMaximum <= BoundaryTolerance;
    std::ofstream statistics("DistanceNormalization.csv");
    statistics << std::setprecision(17)
               << "vertices,triangles,inside_area,max_nodal_error,mean_nodal_error,"
                  "interface_max,relative_residual\n"
               << mesh.getVertexCount() << ',' << mesh.getCellCount() << ','
               << mesh.getArea(Inside) << ',' << maximumError << ',' << meanError << ','
               << interfaceMaximum << ',' << relativeResidual << '\n';
    statistics.flush();
    IO::XDMF output("DistanceNormalization");
    output.grid()
      .setMesh(mesh)
      .add("potential", poisson)
      .add("normalized", normalized)
      .add("reference", reference);
    output.write().close();
    std::cout << "vertices=" << mesh.getVertexCount()
              << " triangles=" << mesh.getCellCount() << " max_error=" << maximumError
              << " mean_error=" << meanError << " interface_max=" << interfaceMaximum
              << '\n';
    return valid && bool(statistics) ? EXIT_SUCCESS : EXIT_FAILURE;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
//! [distance-source]
