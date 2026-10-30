/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [quadrature-source]
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/QF.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main()
{
  //! [quadrature-reference]
  const auto& formula = QF::QuadratureFormula::get(size_t{4}, Polytope::Type::Triangle);
  Real referenceIntegral = 0;
  for (size_t q = 0; q < formula.getSize(); ++q)
  {
    const auto& point = formula.getPoint(q);
    referenceIntegral += formula.getWeight(q) * point.x() * point.y();
  }
  //! [quadrature-reference]

  //! [quadrature-map]
  Mesh mesh = Mesh<Context::Local>::Builder()
    .initialize(2).nodes(3)
    .vertex({1, 1}).vertex({3, 1}).vertex({1, 4})
    .polytope(Polytope::Type::Triangle, {0, 1, 2}).finalize();
  const auto cell = mesh.getCell();
  PolytopeQuadrature mapped(*cell, formula);
  Real physicalIntegral = 0;
  std::ofstream samples("Quadrature.csv");
  samples.precision(17);
  samples << "r,s,w,x,y,jacobian\n";
  for (size_t q = 0; q < mapped.getSize(); ++q)
  {
    const auto& point = mapped.getPoint(q);
    const auto& reference = formula.getPoint(q);
    const Real weight = formula.getWeight(q);
    const Real jacobian = point.getJacobianDeterminant();
    physicalIntegral += weight * jacobian * point.x() * point.y();
    samples << reference.x() << ',' << reference.y() << ',' << weight << ','
            << point.x() << ',' << point.y() << ',' << jacobian << '\n';
  }
  //! [quadrature-map]

  //! [quadrature-check]
  constexpr Real tolerance = 1e-12;
  std::cout.precision(17);
  std::cout << "Reference integral: " << referenceIntegral << '\n'
            << "Physical integral: " << physicalIntegral << '\n'
            << "Quadrature points: " << formula.getSize() << std::endl;
  const bool valid = std::abs(referenceIntegral - 1.0 / 24.0) < tolerance
    && std::abs(physicalIntegral - 9.5) < tolerance;
  //! [quadrature-check]

  std::ofstream vertices("ReferencePolytopes.csv");
  vertices.precision(17);
  vertices << "type,vertex,dimension,x,y,z\n";
  for (const auto type : {Polytope::Type::Segment, Polytope::Type::Triangle,
         Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
         Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
  {
    const Polytope::Traits traits(type);
    for (size_t i = 0; i < traits.getVertexCount(); ++i)
    {
      const auto& point = traits.getVertex(i);
      vertices << static_cast<int>(type) << ',' << i << ',' << traits.getDimension();
      for (size_t j = 0; j < 3; ++j)
        vertices << ',' << (j < traits.getDimension() ? point[j] : 0);
      vertices << '\n';
    }
  }
  samples.flush();
  vertices.flush();
  return valid && samples && vertices ? EXIT_SUCCESS : EXIT_FAILURE;
}
//! [quadrature-source]
