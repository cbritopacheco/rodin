/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [fmm-source]
#include <Rodin/Eikonal/FMM.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Variational.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  std::ofstream errors("FMMErrors.csv");
  errors << std::setprecision(17) << "intervals,speed,max_error,mean_error\n";
  bool valid = true;
  for (const Real propagationSpeed : {Real(1), Real(2)})
  {
    Real previousError = std::numeric_limits<Real>::infinity();
    for (const size_t intervals : {16, 32, 64})
    {
      //! [fmm-mesh]
      const size_t n = intervals + 1;
      const Real h = Real(1) / intervals;
      auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {n, n});
      mesh.scale(h);
      mesh.getConnectivity().compute(2, 0);
      mesh.getConnectivity().compute(0, 2);
      mesh.getConnectivity().compute(0, 0);
      //! [fmm-mesh]
      //! [fmm-solve]
      P1 space(mesh);
      GridFunction arrival(space);
      IndexVector seeds;
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        const auto& x = vertex->getCoordinates();
        if (x(0) == Real(0.5) && x(1) == Real(0.5))
          seeds.push_back(vertex->getIndex());
      }
      if (seeds.size() != 1)
        return EXIT_FAILURE;
      const auto speed = [propagationSpeed](const Point&) { return propagationSpeed; };
      Eikonal::FMM marching(arrival, speed);
      marching.seed(seeds).solve();
      //! [fmm-solve]
      //! [fmm-error]
      Real maximum = 0, mean = 0;
      for (auto vertex = mesh.getVertex(); vertex; ++vertex)
      {
        const auto& x = vertex->getCoordinates();
        const auto& dofs = space.getDOFs(0, vertex->getIndex());
        const Real computed = arrival.getData()(dofs[0]);
        const Real exact =
          std::hypot(x(0) - Real(0.5), x(1) - Real(0.5)) / propagationSpeed;
        const Real error = std::abs(computed - exact);
        valid = valid && std::isfinite(computed) && computed >= 0;
        maximum = std::max(maximum, error);
        mean += error;
      }
      mean /= mesh.getVertexCount();
      // A deliberately loose O(h) regression envelope on this mesh sequence.
      // Measured errors, rather than this bound, characterize the observed rate.
      valid =
        valid && maximum < previousError && maximum <= Real(2) * h / propagationSpeed;
      previousError = maximum;
      errors << intervals << ',' << propagationSpeed << ',' << maximum << ',' << mean
             << '\n';
      //! [fmm-error]
      if (intervals == 64 && propagationSpeed == 1)
      {
        GridFunction reference(space);
        reference = [](const Point& point) {
          return std::hypot(point.x() - Real(0.5), point.y() - Real(0.5));
        };
        IO::XDMF output("FMMArrival");
        output.grid().setMesh(mesh).add("arrival", arrival).add("reference", reference);
        output.write().close();
      }
      std::cout << "intervals=" << intervals << " speed=" << propagationSpeed
                << " max_error=" << maximum << " mean_error=" << mean << '\n';
    }
  }
  errors.flush();
  return valid && bool(errors) ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [fmm-source]
