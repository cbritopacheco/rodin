/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI finite pressure spectra independent of Stokes solves. */
#include "../../PETScStokesProblem.h"
#include "../../CurvedGeometry.h"
#include "../../PETScMixedStability.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::H::PETScStokesStability
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  TEST(PETScMixedStabilityTest, GlobalLogicalConstraintAndIndependentPressureMetric)
  {
    std::array<Mat, 3> operators{nullptr, nullptr, nullptr};
    for (auto& op : operators)
      ASSERT_EQ(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
        2, 2, 2, nullptr, 2, nullptr, &op), PETSC_SUCCESS);
    const std::array<std::array<Real, 4>, 3> values{{
      {4, 0, 0, 9}, {2, 7, -2, -7}, {1, 0, 0, 3}}};
    const PetscInt columns[]{0, 1};
    IndexMap<Real> constrained;
    for (size_t i = 0; i < operators.size(); ++i)
    {
      PetscInt first = 0, last = 0;
      EXPECT_EQ(MatGetOwnershipRange(operators[i], &first, &last), PETSC_SUCCESS);
      for (PetscInt row = first; row < last; ++row)
      {
        EXPECT_EQ(MatSetValues(operators[i], 1, &row, 2, columns,
          values[i].data() + 2 * row, INSERT_VALUES), PETSC_SUCCESS);
        if (i == 0 && row == 1) constrained.emplace(1, 0);
      }
      EXPECT_EQ(MatAssemblyBegin(operators[i], MAT_FINAL_ASSEMBLY), PETSC_SUCCESS);
      EXPECT_EQ(MatAssemblyEnd(operators[i], MAT_FINAL_ASSEMBLY), PETSC_SUCCESS);
    }
    MixedStability::Result result;
    PETScMixedStability::compute(result, operators[0], operators[1], operators[2], constrained);
    for (auto& op : operators) EXPECT_EQ(MatDestroy(&op), PETSC_SUCCESS);
    ASSERT_FALSE(::testing::Test::HasFatalFailure());
    MixedStability::expectConsistent(result);
    ASSERT_EQ(result.freeVelocity, 1);
    ASSERT_EQ(result.zeroMeanPressure, 1);
    // T=(-3,1)^T, A_free=4, B_free=(2,-2)^T: Schur=16, mass=12.
    EXPECT_NEAR(result.eigenvalues(0), Real(4) / 3, MixedStability::ConsistencyTolerance);
  }

  template <class ContextType>
  class StabilityTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      static constexpr size_t AssemblyOrder = 16, RefinedOrder = 20;
      // Relative integration-contamination policy, not a uniform inf-sup bound.
      static constexpr Real QuadratureTolerance = 1e-8;

      void measure(MixedStability::Result& result, size_t n, bool curved,
        size_t order, bool omitDivergence = false) const
      {
        const auto initialize = [curved](auto& mesh) {
          if (curved)
          {
            CurvedGeometry mapping(mesh);
            mapping.template install<2>();
          }
        };
        auto mesh = [&] {
          if constexpr (std::is_same_v<ContextType, Context::Local>)
          {
            auto local = UniformGrid(this->GetParam()).makeMesh(n);
            initialize(local);
            return local;
          }
#ifdef RODIN_USE_MPI
          else
            return DistributedUniformGrid(Context::MPI(*environment, *world), this->GetParam())
              .makeMesh(n, initialize);
#endif
        }();
        H1 velocitySpace(std::integral_constant<size_t, 2>{}, mesh, mesh.getDimension());
        H1 pressureSpace(std::integral_constant<size_t, 1>{}, mesh);
        PETSc::Variational::TrialFunction u(velocitySpace);
        PETSc::Variational::TrialFunction p(pressureSpace);
        PETSc::Variational::TestFunction v(velocitySpace);
        PETSc::Variational::TestFunction q(pressureSpace);
        auto diffusion = Integral(Jacobian(u), Jacobian(v));
        auto divergence = Integral(Real(omitDivergence ? 0 : 1) * Div(u), q);
        auto pressureMass = Integral(p, q);
        diffusion.setOrder(order);
        divergence.setOrder(order);
        pressureMass.setOrder(order);
        BilinearForm a(u, v);
        BilinearForm b(u, q);
        BilinearForm mass(p, q);
        a = diffusion;
        b = divergence;
        mass = pressureMass;
        a.assemble();
        b.assemble();
        mass.assemble();
        auto boundary = DirichletBC(u, Zero());
        boundary.assemble();
        PETScMixedStability::compute(result, a.getOperator(), b.getOperator(),
          mass.getOperator(), std::get<IndexMap<Real>>(boundary.getDOFs()));
      }

      void checkHierarchy() const
      {
        for (bool curved : {false, true})
          for (size_t n : {2u, 3u, 5u})
          {
            SCOPED_TRACE(::testing::Message() << "curved=" << curved << " n=" << n);
            std::array<MixedStability::Result, 2> measurements;
            for (size_t i = 0; i < measurements.size(); ++i)
            {
              measure(measurements[i], n, curved, i == 0 ? AssemblyOrder : RefinedOrder);
              ASSERT_FALSE(::testing::Test::HasFatalFailure());
              MixedStability::expectConsistent(measurements[i]);
              const bool obstructed = n == 2 && this->GetParam() != Polytope::Type::Pyramid;
              EXPECT_EQ(measurements[i].isDimensionObstructed(), obstructed);
              EXPECT_EQ(MixedStability::hasResolvedPositiveSpectrum(measurements[i]), !obstructed);
              SCOPED_TRACE(::testing::Message() << "free velocity=" << measurements[i].freeVelocity
                << "zero-mean pressure=" << measurements[i].zeroMeanPressure
                << "smallest eigenvalue=" << measurements[i].eigenvalues.minCoeff());
            }
            ASSERT_EQ(measurements[0].freeVelocity, measurements[1].freeVelocity);
            ASSERT_EQ(measurements[0].zeroMeanPressure, measurements[1].zeroMeanPressure);
            if (!measurements[0].isDimensionObstructed())
            {
              ASSERT_GT(measurements[0].eigenvalues.minCoeff(), 0);
              EXPECT_LT(std::abs(measurements[1].eigenvalues.minCoeff() /
                measurements[0].eigenvalues.minCoeff() - 1), QuadratureTolerance);
            }
          }
      }

      void checkMissingDivergence() const
      {
        MixedStability::Result correct, wrong;
        measure(correct, 3, true, AssemblyOrder);
        measure(wrong, 3, true, AssemblyOrder, true);
        ASSERT_FALSE(::testing::Test::HasFatalFailure());
        MixedStability::expectConsistent(correct);
        MixedStability::expectConsistent(wrong);
        EXPECT_TRUE(MixedStability::hasResolvedPositiveSpectrum(correct));
        EXPECT_FALSE(MixedStability::hasResolvedPositiveSpectrum(wrong));
        EXPECT_EQ(correct.freeVelocity, wrong.freeVelocity);
        EXPECT_EQ(correct.zeroMeanPressure, wrong.zeroMeanPressure);
        EXPECT_TRUE(wrong.eigenvalues.isZero(0));
      }
  };

  using LocalStabilityTest = StabilityTest<Context::Local>;
  TEST_P(LocalStabilityTest, PressureSpectrumAcrossRefinementLevels) { checkHierarchy(); }
  TEST_P(LocalStabilityTest, MissingDivergenceRejected) { checkMissingDivergence(); }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalStabilityTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) { return std::string(UniformGrid::getGeometryName(info.param)); });
#ifdef RODIN_USE_MPI
  using MPIStabilityTest = StabilityTest<Context::MPI>;
  TEST_P(MPIStabilityTest, PressureSpectrumAcrossRefinementLevels) { checkHierarchy(); }
  TEST_P(MPIStabilityTest, MissingDivergenceRejected) { checkMissingDivergence(); }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPIStabilityTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) { return std::string(UniformGrid::getGeometryName(info.param)); });
#endif
}

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS) return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::H::PETScStokesStability::environment = &env;
    Rodin::Tests::Convergence::H::PETScStokesStability::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS) return 1;
  return result;
}
