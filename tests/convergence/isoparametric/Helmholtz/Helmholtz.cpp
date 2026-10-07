/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Native-complex Helmholtz on exact and approximated geometry. */

#include "Helmholtz.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::Helmholtz
{
  template <size_t Q = 2>
  class NativeProblem
  {
    public:
      static constexpr size_t GeometryDegree = Q;

      NativeProblem(Polytope::Type geometry, size_t n,
        CurvedGeometry<LocalMesh>::Map map = CurvedGeometry<LocalMesh>::Map::Quadratic,
        bool lifted = false, Real amplitude = 0.1)
        : m_mesh(UniformGrid(geometry).makeMesh(n)),
          m_reference(),
          m_geometry(m_mesh, map, amplitude)
      {
        if (lifted)
          m_reference.emplace(m_mesh);
        m_geometry.install<Q>();
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }
      const auto& getGeometry() const
      {
        return m_geometry;
      }
      const auto& getReference() const
      {
        return m_reference.value();
      }

      template <size_t K>
      ErrorNorms solve(HelmholtzData::Field field, bool omitMass = false,
        size_t order = AssemblyOrder, Real tolerance = 1e-13, size_t normOrder = 0,
        LiftedErrorNorm::Result* lifted = nullptr) const
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
        const size_t integrationOrder = normOrder == 0 ? order + 2 : normOrder;
        if (lifted)
        {
          assert(m_reference);
          *lifted = LiftedErrorNorm::compute(
            *m_reference, m_mesh, u.getSolution(), data, SineMap(), integrationOrder);
        }
        return ErrorNorm::compute(m_mesh, u.getSolution(), data.getSolution(),
          data.getGradient(), integrationOrder);
      }

    private:
      LocalMesh m_mesh;
      Optional<LocalMesh> m_reference;
      CurvedGeometry<LocalMesh> m_geometry;
  };

  using NativeHelmholtzTest = HelmholtzTest<NativeProblem<>>;
  using NativeHelmholtzQ1Test = HelmholtzTest<NativeProblem<1>>;
  using NativeHelmholtzQ3Test = HelmholtzTest<NativeProblem<3>>;
  TEST_P(NativeHelmholtzQ1Test, LiftedAffineRates)
  {
    checkMatchedGeometryRates();
  }
  TEST_P(NativeHelmholtzQ1Test, LiftedIndependentSensitivity)
  {
    checkMatchedGeometrySensitivity();
  }
  TEST_P(NativeHelmholtzQ3Test, LiftedAffineRates)
  {
    checkMatchedGeometryRates();
  }
  TEST_P(NativeHelmholtzQ3Test, LiftedIndependentSensitivity)
  {
    checkMatchedGeometrySensitivity();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, NativeHelmholtzQ1Test,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
  INSTANTIATE_TEST_SUITE_P(AllGeometries, NativeHelmholtzQ3Test,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
  TEST_P(NativeHelmholtzTest, P1OptimalRates)
  {
    checkRates<1>();
  }
  TEST_P(NativeHelmholtzTest, ComplexLiftedMetricOracle)
  {
    checkComplexLiftedMetric();
  }
  TEST_P(NativeHelmholtzTest, ApproximatedP1Rates)
  {
    checkApproximatedRates<1>();
  }
  TEST_P(NativeHelmholtzTest, ApproximatedP2Rates)
  {
    checkApproximatedRates<2>();
  }
  TEST_P(NativeHelmholtzTest, ApproximatedP1Sensitivity)
  {
    checkApproximatedSensitivity<1>();
  }
  TEST_P(NativeHelmholtzTest, ApproximatedP2Sensitivity)
  {
    checkApproximatedSensitivity<2>();
  }
  TEST_P(NativeHelmholtzTest, ApproximatedAffinePatchRejectsOmittedMass)
  {
    checkApproximatedPatchAndControl();
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
