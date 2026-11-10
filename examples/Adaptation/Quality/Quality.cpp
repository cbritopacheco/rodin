/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <Rodin/Adaptation/CellDeformation.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  try
  {
    const std::filesystem::path directory = argc > 1 ? argv[1] : "quality";
    std::filesystem::create_directories(directory);
    constexpr size_t VerticesPerAxis = 3;
    constexpr size_t Dimension = 2;
    constexpr size_t PlotSubdivisions = 32;
    constexpr Real Shear = 0.5;
    constexpr Real Bow = 0.4;
    auto mesh = LocalMesh::UniformGrid(
      Polytope::Type::Triangle, {VerticesPerAxis, VerticesPerAxis});
    mesh.scale(Real(1) / (VerticesPerAxis - 1));
    for (size_t from = 0; from <= 2; ++from)
    {
      for (size_t to = 0; to <= 2; ++to)
        mesh.getConnectivity().compute(from, to);
    }

    // Interpolate one prescribed deformation in each displacement space.
    VectorFunction deformation(Dimension, [](const Point& point) {
      const auto& x = point.getPhysicalCoordinates();
      Math::SpatialVector<Real> value(2);
      value(0) = Shear * x(1);
      value(1) = Bow * x(0) * (1 - x(0));
      return value;
    });
    const auto run = [&]<size_t Degree>() {
      H1 space(std::integral_constant<size_t, Degree>{}, mesh, 2);
      GridFunction displacement(space);
      displacement = deformation;
      const auto gradient = Jacobian(displacement);

      // Construct isoparametric geometry from field values, not coefficients.
      LocalMesh moved(mesh);
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        const Point point(*vertex, Polytope::Traits(Polytope::Type::Point).getVertex(0));
        moved.setVertexCoordinates(vertex->getIndex(),
          point.getPhysicalCoordinates() + displacement.getValue(point));
      }
      for (size_t d = 1; d <= 2; ++d)
      {
        for (auto entity = mesh.getPolytope(d); entity; ++entity)
        {
          RealH1Element<Degree> geometry(entity->getGeometry());
          PointCloud nodes(2, geometry.getCount());
          for (size_t node = 0; node < geometry.getCount(); ++node)
          {
            const Point point(*entity, geometry.getNode(node));
            const Math::SpatialPoint position(
              point.getPhysicalCoordinates() + displacement.getValue(point));
            for (size_t component = 0; component < Dimension; ++component)
              nodes(component, node) = position(component);
          }
          moved.setPolytopeTransformation({d, entity->getIndex()},
            new ParametricTransformation<RealH1Element<Degree>>(
              std::move(nodes), geometry));
        }
      }

      const auto stem = directory / ("P" + std::to_string(Degree));
      H1 movedSpace(std::integral_constant<size_t, Degree>{}, moved, 2);
      GridFunction movedDisplacement(movedSpace);
      movedDisplacement.getData() = displacement.getData();
      IO::XDMF output(stem.string());
      output.grid("background").setMesh(mesh).add("displacement", displacement);
      output.grid("moved").setMesh(moved).add("displacement", movedDisplacement);
      output.write().close();

      // Sample each cell separately: P1 gradients can jump across shared edges.
      std::ofstream samples(stem.string() + ".csv");
      samples.exceptions(std::ios::failbit | std::ios::badbit);
      samples << "cell,i,k,x0,y0,x,y,j,Q\n"
              << std::setprecision(std::numeric_limits<Real>::max_digits10);
      Real minJ = std::numeric_limits<Real>::infinity(), maxQ = 1;
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto movedCell = moved.getCell(cell->getIndex());
        for (size_t i = 0; i <= PlotSubdivisions; ++i)
        {
          for (size_t k = 0; k <= PlotSubdivisions - i; ++k)
          {
            Math::SpatialPoint reference(2);
            reference(0) = Real(i) / PlotSubdivisions;
            reference(1) = Real(k) / PlotSubdivisions;
            const Point point(*cell, reference), image(*movedCell, reference);
            Adaptation::CellDeformation state(2);
            state.setDisplacementGradient(gradient.getValue(point));
            const Real j = state.getJacobian(), Q = state.getRelativeDistortion();
            // Independently verify that the installed geometry has this gradient.
            const Math::SpatialMatrix<Real> F(
              image.getJacobian() * point.getJacobian().inverse());
            constexpr Real VerificationTolerance = 1e-10;
            if ((F - state.getDeformationGradient()).norm() > VerificationTolerance ||
              j <= 0)
              throw std::runtime_error("Deformation geometry or orientation mismatch.");
            const auto& x0 = point.getPhysicalCoordinates();
            const auto& x = image.getPhysicalCoordinates();
            samples << cell->getIndex() << ',' << i << ',' << k << ',' << x0(0) << ','
                    << x0(1) << ',' << x(0) << ',' << x(1) << ',' << j << ',' << Q
                    << '\n';
            minJ = std::min(minJ, j);
            maxQ = std::max(maxQ, Q);
          }
        }
      }
      std::cout << "P" << Degree << ": vertices=" << mesh.getVertexCount()
                << " cells=" << mesh.getCellCount() << " min_j=" << minJ
                << " max_Q=" << maxQ << " output=" << stem << '\n';
    };
    run.template operator()<1>();
    run.template operator()<2>();
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
