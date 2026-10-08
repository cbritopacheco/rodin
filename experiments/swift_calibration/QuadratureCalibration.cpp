/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "QuadratureAudit.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Adaptation;

class Calibration
{
  public:
    template <size_t Degree>
    void run(size_t dimension, bool stress) const
    {
      using LocalMesh = Mesh<Context::Local>;
      auto mesh = dimension == 2
        ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3})
        : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
      mesh.scale(Real(0.5));
      for (size_t from = 1; from <= dimension; ++from)
        for (size_t to = 0; to <= dimension; ++to)
          if (from != to)
            mesh.getConnectivity().compute(from, to);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        const auto& qf = QF::PolytopeQuadratureFormula::get(1, cell->getGeometry());
        const auto& point = cell->getQuadrature(qf).getPoint(0);
        mesh.setAttribute({dimension, cell->getIndex()}, point.x() < Real(0.5) ? 1 : 2);
      }
      size_t facets = 0;
      for (auto face = mesh.getFace(); face; ++face)
      {
        const Polytope::Traits traits(face->getGeometry());
        bool interface = true;
        for (size_t i = 0; i < traits.getVertexCount(); ++i)
          interface &= std::abs(Point(*face, traits.getVertex(i)).x() - Real(0.5)) < Real(1e-12);
        if (interface)
        {
          mesh.setAttribute({dimension - 1, face->getIndex()}, 10);
          ++facets;
        }
      }
      H1 space(std::integral_constant<size_t, Degree>{}, mesh, dimension);
      GridFunction current(space);
      current = AnalyticVectorFunction([=](const Point& point) {
        Math::SpatialVector<Real> value = Math::SpatialVector<Real>::Zero(dimension);
        const Real x = point.x(), y = point.y();
        if (stress)
        {
          const Real r = x - Real(0.37);
          value(0) = -Real(0.85) * (r - r * r * r / Real(3));
          value(1) = Real(0.6) * (x - Real(0.5)) * y * (1 - y);
        }
        else
        {
          value(0) = Real(0.08) * x * (1 - x);
          value(1) = Real(0.02) * y * (1 - y);
        }
        return value;
      }, dimension);
      Adaptation::SWIFT::Parameters p;
      p.model.h = Real(0.5);
      p.model.fit = 1;
      p.model.distribution.deviatoric = Real(1e-4);
      p.model.distribution.divergence = Real(1e-2);
      p.model.hinge = 10;
      p.interfaceAttribute = 10;
      std::cout << std::setprecision(17) << "audit case dimension=" << dimension
                << " degree=" << Degree << " stress=" << stress << " facets=" << facets << '\n';
      Experiments::QuadratureAudit audit(current, p);
      audit.volume();
      const RealFunction phi([=](const Point& point) {
        return point.x() - Real(0.55) + Real(0.02) * std::sin(6 * point.y()) +
          (dimension == 3 ? Real(0.01) * std::sin(5 * point.getCoordinates()(2)) : Real(0));
      });
      const AnalyticVectorFunction grad([=](const Point& point) {
        Math::SpatialVector<Real> value = Math::SpatialVector<Real>::Zero(dimension);
        value(0) = 1;
        value(1) = Real(0.12) * std::cos(6 * point.y());
        if (dimension == 3)
          value(2) = Real(0.05) * std::cos(5 * point.getCoordinates()(2));
        return value;
      }, dimension);
      std::cout << "audit target kind=analytic\n";
      audit.surface(phi, grad, Real(0.2), Real(1));
      H1 targetSpace(std::integral_constant<size_t, Degree>{}, mesh);
      GridFunction target(targetSpace);
      target = phi;
      const RealFunction adapter([&](const Point& point) { return target.getValue(point); });
      auto targetGradient = Grad(target);
      targetGradient.traceOf(1);
      std::cout << "audit target kind=fe\n";
      audit.surface(adapter, targetGradient, Real(0.2), Real(1));
    }
};

int main(int argc, char** argv)
{
  if (argc != 4)
    return 2;
  const size_t dimension = std::stoul(argv[1]), degree = std::stoul(argv[2]);
  const bool stress = std::stoi(argv[3]);
  Calibration calibration;
  switch (degree)
  {
    case 1: calibration.run<1>(dimension, stress); break;
    case 2: calibration.run<2>(dimension, stress); break;
    case 3: calibration.run<3>(dimension, stress); break;
    default: return 2;
  }
}
