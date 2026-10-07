/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI complex Helmholtz on exact/approximated maps. */

#include "Helmholtz.h"
#include "Rodin/Assembly.h"
#include "Rodin/PETSc.h"

#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

#ifndef PETSC_USE_COMPLEX
#error "Curved PETSc Helmholtz requires complex-scalar PETSc"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::Helmholtz
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType, size_t Q = 2>
  class PETScProblem
  {
    public:
      static constexpr size_t GeometryDegree = Q;

      PETScProblem(Polytope::Type geometry, size_t n,
        typename CurvedGeometry<Mesh<ContextType>>::Map map =
          CurvedGeometry<Mesh<ContextType>>::Map::Quadratic,
        bool lifted = false, Real amplitude = 0.1)
        : m_mesh(makeMesh(geometry, n)),
          m_reference(),
          m_geometry(m_mesh, map, amplitude)
      {
        if (lifted)
          m_reference.emplace(m_mesh);
        m_geometry.template install<Q>();
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
        const HelmholtzData data(m_mesh.getSpaceDimension(), field);
        H1<K, Complex, Mesh<ContextType>> space(
          std::integral_constant<size_t, K>{}, m_mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        auto stiffness = Integral(Grad(u), Grad(v));
        auto mass = Integral((omitMass ? Real(0) : Real(0.25)) * u, v);
        auto load = Integral(data.getSource(), v);
        stiffness.setOrder(order);
        mass.setOrder(order);
        load.setOrder(order);
        Problem problem(u, v);
        problem = stiffness - mass - load + DirichletBC(u, data.getSolution());
        PETSc::Solver::CG solver(problem);
        solver.setTolerances(tolerance, 1e-14, 1e5, 20000);
        solver.solve();
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
        auto& system = problem.getLinearSystem();
        Vec residual = nullptr;
        EXPECT_EQ(VecDuplicate(system.getVector(), &residual), PETSC_SUCCESS);
        EXPECT_EQ(
          MatMult(system.getOperator(), system.getSolution(), residual), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(residual, -1, system.getVector()), PETSC_SUCCESS);
        PetscReal norm = 0, rhsNorm = 0;
        EXPECT_EQ(VecNorm(residual, NORM_2, &norm), PETSC_SUCCESS);
        EXPECT_EQ(VecNorm(system.getVector(), NORM_2, &rhsNorm), PETSC_SUCCESS);
        const Real relative = norm / std::max(Real(1), rhsNorm);
        EXPECT_TRUE(std::isfinite(relative));
        EXPECT_LT(relative, 1e-11);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
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
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#ifdef RODIN_USE_MPI
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
            .makeMesh(n);
#endif
      }

      Mesh<ContextType> m_mesh;
      Optional<Mesh<ContextType>> m_reference;
      CurvedGeometry<Mesh<ContextType>> m_geometry;
  };

  using PETScHelmholtzLocalQ1Test = HelmholtzTest<PETScProblem<Context::Local, 1>>;
  TEST_P(PETScHelmholtzLocalQ1Test, LiftedAffineRates)
  {
    checkMatchedGeometryRates();
  }
  TEST_P(PETScHelmholtzLocalQ1Test, LiftedIndependentSensitivity)
  {
    checkMatchedGeometrySensitivity();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScHelmholtzLocalQ1Test,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

  using PETScHelmholtzLocalQ3Test = HelmholtzTest<PETScProblem<Context::Local, 3>>;
  TEST_P(PETScHelmholtzLocalQ3Test, LiftedAffineRates)
  {
    checkMatchedGeometryRates();
  }
  TEST_P(PETScHelmholtzLocalQ3Test, LiftedIndependentSensitivity)
  {
    checkMatchedGeometrySensitivity();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScHelmholtzLocalQ3Test,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

  using PETScHelmholtzLocalTest = HelmholtzTest<PETScProblem<Context::Local>>;
  TEST_P(PETScHelmholtzLocalTest, P1OptimalRates)
  {
    checkRates<1>();
  }
  TEST_P(PETScHelmholtzLocalTest, ComplexLiftedMetricOracle)
  {
    checkComplexLiftedMetric();
  }
  TEST_P(PETScHelmholtzLocalTest, ApproximatedP1Rates)
  {
    checkApproximatedRates<1>();
  }
  TEST_P(PETScHelmholtzLocalTest, ApproximatedP2Rates)
  {
    checkApproximatedRates<2>();
  }
  TEST_P(PETScHelmholtzLocalTest, ApproximatedP1Sensitivity)
  {
    checkApproximatedSensitivity<1>();
  }
  TEST_P(PETScHelmholtzLocalTest, ApproximatedP2Sensitivity)
  {
    checkApproximatedSensitivity<2>();
  }
  TEST_P(PETScHelmholtzLocalTest, ApproximatedAffinePatchRejectsOmittedMass)
  {
    checkApproximatedPatchAndControl();
  }
  TEST_P(PETScHelmholtzLocalTest, P2OptimalRates)
  {
    checkRates<2>();
  }
  TEST_P(PETScHelmholtzLocalTest, ConstantP1Patch)
  {
    checkPatch<1>();
  }
  TEST_P(PETScHelmholtzLocalTest, AffinePhysicalP2Patch)
  {
    checkPatch<2>();
  }
  TEST_P(PETScHelmholtzLocalTest, AffineP2PatchRejectsOmittedMass)
  {
    checkNegativeControl();
  }
  TEST_P(PETScHelmholtzLocalTest, P1QuadratureAndSolverSensitivity)
  {
    checkSensitivity<1>();
  }
  TEST_P(PETScHelmholtzLocalTest, P2QuadratureAndSolverSensitivity)
  {
    checkSensitivity<2>();
  }
  TEST_P(PETScHelmholtzLocalTest, GeometryAndPhysicalVolume)
  {
    checkGeometry();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScHelmholtzLocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  using PETScHelmholtzMPIQ1Test = HelmholtzTest<PETScProblem<Context::MPI, 1>>;
  TEST_P(PETScHelmholtzMPIQ1Test, LiftedAffineRates)
  {
    checkMatchedGeometryRates();
  }
  TEST_P(PETScHelmholtzMPIQ1Test, LiftedIndependentSensitivity)
  {
    checkMatchedGeometrySensitivity();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScHelmholtzMPIQ1Test,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

  using PETScHelmholtzMPIQ3Test = HelmholtzTest<PETScProblem<Context::MPI, 3>>;
  TEST_P(PETScHelmholtzMPIQ3Test, LiftedAffineRates)
  {
    checkMatchedGeometryRates();
  }
  TEST_P(PETScHelmholtzMPIQ3Test, LiftedIndependentSensitivity)
  {
    checkMatchedGeometrySensitivity();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScHelmholtzMPIQ3Test,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

  using PETScHelmholtzMPITest = HelmholtzTest<PETScProblem<Context::MPI>>;
  TEST_P(PETScHelmholtzMPITest, P1OptimalRates)
  {
    checkRates<1>();
  }
  TEST_P(PETScHelmholtzMPITest, ComplexLiftedMetricOracle)
  {
    checkComplexLiftedMetric();
  }
  TEST_P(PETScHelmholtzMPITest, ApproximatedP1Rates)
  {
    checkApproximatedRates<1>();
  }
  TEST_P(PETScHelmholtzMPITest, ApproximatedP2Rates)
  {
    checkApproximatedRates<2>();
  }
  TEST_P(PETScHelmholtzMPITest, ApproximatedP1Sensitivity)
  {
    checkApproximatedSensitivity<1>();
  }
  TEST_P(PETScHelmholtzMPITest, ApproximatedP2Sensitivity)
  {
    checkApproximatedSensitivity<2>();
  }
  TEST_P(PETScHelmholtzMPITest, ApproximatedAffinePatchRejectsOmittedMass)
  {
    checkApproximatedPatchAndControl();
  }
  TEST_P(PETScHelmholtzMPITest, P2OptimalRates)
  {
    checkRates<2>();
  }
  TEST_P(PETScHelmholtzMPITest, ConstantP1Patch)
  {
    checkPatch<1>();
  }
  TEST_P(PETScHelmholtzMPITest, AffinePhysicalP2Patch)
  {
    checkPatch<2>();
  }
  TEST_P(PETScHelmholtzMPITest, AffineP2PatchRejectsOmittedMass)
  {
    checkNegativeControl();
  }
  TEST_P(PETScHelmholtzMPITest, P1QuadratureAndSolverSensitivity)
  {
    checkSensitivity<1>();
  }
  TEST_P(PETScHelmholtzMPITest, P2QuadratureAndSolverSensitivity)
  {
    checkSensitivity<2>();
  }
  TEST_P(PETScHelmholtzMPITest, GeometryAndPhysicalVolume)
  {
    checkGeometry();
  }
  TEST_P(PETScHelmholtzMPITest, P2GlobalNormCountsOwnedCurvedCellsOnce)
  {
    PETScProblem<Context::MPI> problem(GetParam(), 2);
    const auto& mesh = problem.getMesh();
    H1<2, Complex, Mesh<Context::MPI>> space(std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::GridFunction u(space);
    u = ComplexFunction(Complex(0));
    const auto zeroGradient = [](const Point& p) {
      return Math::SpatialVector<Complex>::Zero(p.getCoordinates().size());
    };
    const auto error =
      ErrorNorm::compute(mesh, u, ComplexFunction(Complex(1, 1)), zeroGradient, 14);
    const Real volume = UniformGrid(GetParam()).getDimension() == 1 ? 1.1 : 1;
    EXPECT_NEAR(error.getL2(), std::sqrt(2 * volume), 1e-12);
    EXPECT_EQ(error.getH1Seminorm(), 0);
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScHelmholtzMPITest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#endif
}

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::Helmholtz::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::Helmholtz::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
