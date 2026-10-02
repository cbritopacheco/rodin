/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Native-complex Helmholtz P1/P2 convergence on exact P2 geometry. */

#include "Helmholtz.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::Helmholtz
{
  class NativeProblem
  {
    public:
      NativeProblem(Polytope::Type geometry, size_t n)
        : m_mesh(UniformGrid(geometry).makeMesh(n)),
          m_geometry(m_mesh)
      {
        m_geometry.install<2>();
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }
      const auto& getGeometry() const
      {
        return m_geometry;
      }

      template <size_t K>
      ErrorNorms solve(HelmholtzData::Field field, bool omitMass = false,
        size_t order = AssemblyOrder, Real tolerance = 1e-13) const
      {
        const HelmholtzData data(m_mesh.getDimension(), field);
        auto space = [&] {
          if constexpr (K == 1)
            return P1<Complex>(m_mesh);
          else
            return H1<K, Complex>(std::integral_constant<size_t, K>{}, m_mesh);
        }();
        TrialFunction u(space);
        TestFunction v(space);
        auto stiffness = Integral(Grad(u), Grad(v));
        auto mass = Integral((omitMass ? Real(0) : Real(0.25)) * u, v);
        auto load = Integral(data.getSource(), v);
        stiffness.setOrder(order);
        mass.setOrder(order);
        load.setOrder(order);
        Problem problem(u, v);
        problem = stiffness - mass - load + DirichletBC(u, data.getSolution());
        Solver::CG solver(problem);
        solver.setTolerance(tolerance).setMaxIterations(20000).solve();
        EXPECT_TRUE(solver.success());
        const auto& system = problem.getLinearSystem();
        const auto residual =
          system.getOperator() * system.getSolution() - system.getVector();
        const Real relative =
          residual.norm() / std::max(Real(1), system.getVector().norm());
        EXPECT_TRUE(std::isfinite(relative));
        EXPECT_LT(relative, 1e-11);
        return ErrorNorm::compute(
          m_mesh, u.getSolution(), data.getSolution(), data.getGradient(), order + 2);
      }

    private:
      LocalMesh m_mesh;
      CurvedGeometry<LocalMesh> m_geometry;
  };

  using NativeHelmholtzTest = HelmholtzTest<NativeProblem>;
  TEST_P(NativeHelmholtzTest, P1OptimalRates)
  {
    checkRates<1>();
  }
  TEST_P(NativeHelmholtzTest, P2OptimalRates)
  {
    checkRates<2>();
  }
  TEST_P(NativeHelmholtzTest, ConstantP1Patch)
  {
    checkPatch<1>();
  }
  TEST_P(NativeHelmholtzTest, AffinePhysicalP2Patch)
  {
    checkPatch<2>();
  }
  TEST_P(NativeHelmholtzTest, AffineP2PatchRejectsOmittedMass)
  {
    checkNegativeControl();
  }
  TEST_P(NativeHelmholtzTest, P1QuadratureAndSolverSensitivity)
  {
    checkSensitivity<1>();
  }
  TEST_P(NativeHelmholtzTest, P2QuadratureAndSolverSensitivity)
  {
    checkSensitivity<2>();
  }
  TEST_P(NativeHelmholtzTest, GeometryAndPhysicalVolume)
  {
    checkGeometry();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, NativeHelmholtzTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
}
