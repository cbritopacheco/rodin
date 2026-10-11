/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Mixed traction and nearly incompressible PETSc verification. */

#include "PETScLinearElasticityProblem.h"
#include "FieldConvergence.h"
#ifdef RODIN_ELASTICITY_BOUNDARY_CURVED
#include "CurvedGeometry.h"
#endif
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::PETScLinearElasticityBoundaryTests
{
  using Data = LinearElasticity::ManufacturedSolution;
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n, bool mixed)
  {
    using Attributes = PETScLinearElasticityProblem<1, LocalMesh>;
    const auto initialize = [mixed](LocalMesh& mesh) {
      if (mixed)
        UnitBoxBoundary::labelCoordinatePartition(mesh, 0, Attributes::DirichletAttribute,
          Attributes::TractionAttribute, Attributes::TractionAttribute);
    };
    if constexpr (std::is_same_v<ContextType, Context::Local>)
    {
      auto mesh = UniformGrid(geometry).makeMesh(n);
      initialize(mesh);
#ifdef RODIN_ELASTICITY_BOUNDARY_CURVED
      CurvedGeometry curved(mesh);
      curved.template install<2>();
#endif
      return mesh;
    }
#ifdef RODIN_USE_MPI
    else
    {
      Context::MPI context(*environment, *world);
      auto mesh = DistributedUniformGrid(context, geometry).makeMesh(n, initialize);
#ifdef RODIN_ELASTICITY_BOUNDARY_CURVED
      CurvedGeometry curved(mesh);
      curved.template install<2>();
#endif
      return mesh;
    }
#endif
  }

  template <class ContextType>
  class Fixture : public ::testing::TestWithParam<Polytope::Type>
  {
    public:
#ifdef RODIN_ELASTICITY_BOUNDARY_CURVED
      static constexpr size_t AffinePatchDegree = 2, QuadraticPatchDegree = 4;
#else
      static constexpr size_t AffinePatchDegree = 1, QuadraticPatchDegree = 2;
#endif
      template <size_t K>
      void checkRates(bool nearly = false) const
      {
        const auto geometry = this->GetParam();
        if (nearly && geometry == Polytope::Type::Segment)
          GTEST_SKIP() << "Nonconstant divergence-free displacement needs dimension >= 2";
        std::array<size_t, 3> levels =
          K == 1 ? std::array<size_t, 3>{5, 9, 17} : std::array<size_t, 3>{3, 5, 9};
        if constexpr (K == 1)
          if (geometry == Polytope::Type::Tetrahedron)
            levels = {9, 17, 33};
        FieldConvergence<1> history;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          const auto mesh = makeMesh<ContextType>(geometry, n, !nearly);
          using Problem =
            PETScLinearElasticityProblem<K, std::remove_cvref_t<decltype(mesh)>>;
          const Problem problem(mesh,
            nearly ? Data::Field::DivergenceFree : Data::Field::Exponential,
            AssemblyOrder,
            nearly ? Problem::Boundary::Dirichlet : Problem::Boundary::MixedTraction,
            nearly ? NearlyLambda : Lambda, nearly ? NearlyMu : Mu, PCJACOBI);
          history.append(Real(1) / Real(n - 1), {problem.solve()});
        }
        history.expectAlgebraicFloor(nearly ? Rates{2.25, 1.35}
            : K == 1                        ? Rates{1.65, 0.75}
                                            : Rates{Real(K) + 0.45, Real(K) - 0.45});
      }

      template <size_t K>
      void checkPatch() const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3, true);
        using Problem =
          PETScLinearElasticityProblem<K, std::remove_cvref_t<decltype(mesh)>>;
        const Problem problem(mesh,
          K == AffinePatchDegree ? Data::Field::AsymmetricAffine : Data::Field::Quadratic,
          AssemblyOrder, Problem::Boundary::MixedTraction, Lambda, Mu, PCJACOBI);
        const auto error = problem.solve();
        EXPECT_TRUE(error.isFinite());
        EXPECT_LT(error.getL2(), PatchTolerance);
        EXPECT_LT(error.getH1Seminorm(), PatchTolerance);
      }

      void checkControl(bool omitVolumetric = false) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3, true);
        using Problem =
          PETScLinearElasticityProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const Problem problem(mesh,
#ifdef RODIN_ELASTICITY_BOUNDARY_CURVED
          Data::Field::AsymmetricAffine,
#else
          Data::Field::Quadratic,
#endif
          AssemblyOrder, Problem::Boundary::MixedTraction, Lambda, Mu, PCJACOBI);
        const auto correct = problem.solve();
        const auto wrong =
          problem.solve(omitVolumetric, SolverTolerance, NormOrder, !omitVolumetric);
        EXPECT_LT(correct.getL2(), PatchTolerance);
        EXPECT_LT(correct.getH1Seminorm(), PatchTolerance);
        EXPECT_TRUE(wrong.isFinite());
        EXPECT_GT(wrong.getL2(), WrongValueFloor);
        EXPECT_GT(wrong.getH1Seminorm(), WrongGradientFloor);
      }

      void checkSensitivity(bool nearly = false) const
      {
        const auto geometry = this->GetParam();
        if (nearly && geometry == Polytope::Type::Segment)
          GTEST_SKIP() << "Nonconstant divergence-free displacement needs dimension >= 2";
        const auto mesh = makeMesh<ContextType>(geometry, 3, !nearly);
        using Problem =
          PETScLinearElasticityProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const auto field =
          nearly ? Data::Field::DivergenceFree : Data::Field::Exponential;
        const auto boundary =
          nearly ? Problem::Boundary::Dirichlet : Problem::Boundary::MixedTraction;
        const Real lambda = nearly ? NearlyLambda : Lambda, mu = nearly ? NearlyMu : Mu;
        const Problem base(mesh, field, AssemblyOrder, boundary, lambda, mu, PCJACOBI);
        const Problem higher(
          mesh, field, AssemblyOrder + 2, boundary, lambda, mu, PCJACOBI);
        const auto reference = base.solve();
        ASSERT_GT(reference.getL2(), 0);
        ASSERT_GT(reference.getH1Seminorm(), 0);
        for (const auto& varied :
          {higher.solve(), base.solve(false, SolverTolerance, NormOrder + 2),
            base.solve(false, SolverTolerance / 10)})
        {
          EXPECT_LT(
            std::abs(varied.getL2() / reference.getL2() - 1), SensitivityTolerance);
          EXPECT_LT(std::abs(varied.getH1Seminorm() / reference.getH1Seminorm() - 1),
            SensitivityTolerance);
        }
      }

    private:
      static constexpr Real Lambda = 1.5, Mu = 0.5;
      // Same nearly incompressible parameters as the native h study.
      static constexpr Real NearlyLambda = 1e4, NearlyMu = 1;
      static constexpr size_t AssemblyOrder = 16, NormOrder = 18;
      static constexpr Real SolverTolerance = 1e-13;
      static constexpr Real PatchTolerance = 1e-9;
      static constexpr Real WrongValueFloor = 1e-3, WrongGradientFloor = 1e-2;
      static constexpr Real SensitivityTolerance = 1e-6;
  };

#ifdef RODIN_ELASTICITY_BOUNDARY_CURVED
#define RODIN_ELASTICITY_BOUNDARY_EXTRA(Name)                                            \
  TEST_P(Name, MixedTractionP3Rates)                                                     \
  {                                                                                      \
    checkRates<3>();                                                                     \
  }                                                                                      \
  TEST_P(Name, RejectsMissingVolumetricTerm)                                             \
  {                                                                                      \
    checkControl(true);                                                                  \
  }
#else
#define RODIN_ELASTICITY_BOUNDARY_EXTRA(Name)                                            \
  TEST_P(Name, NearlyIncompressibleP2Rates)                                              \
  {                                                                                      \
    checkRates<2>(true);                                                                 \
  }                                                                                      \
  TEST_P(Name, NearlyIncompressibleIndependentSensitivity)                               \
  {                                                                                      \
    checkSensitivity(true);                                                              \
  }
#endif

#define RODIN_ELASTICITY_BOUNDARY_TESTS(Name)                                            \
  TEST_P(Name, MixedTractionP1Rates)                                                     \
  {                                                                                      \
    checkRates<1>();                                                                     \
  }                                                                                      \
  TEST_P(Name, MixedTractionP2Rates)                                                     \
  {                                                                                      \
    checkRates<2>();                                                                     \
  }                                                                                      \
  TEST_P(Name, AsymmetricAffineTractionPatch)                                            \
  {                                                                                      \
    checkPatch<AffinePatchDegree>();                                                     \
  }                                                                                      \
  TEST_P(Name, QuadraticTractionPatch)                                                   \
  {                                                                                      \
    checkPatch<QuadraticPatchDegree>();                                                  \
  }                                                                                      \
  TEST_P(Name, RejectsMissingTraction)                                                   \
  {                                                                                      \
    checkControl();                                                                      \
  }                                                                                      \
  TEST_P(Name, MixedTractionIndependentSensitivity)                                      \
  {                                                                                      \
    checkSensitivity();                                                                  \
  }                                                                                      \
  RODIN_ELASTICITY_BOUNDARY_EXTRA(Name)                                                  \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,                 \
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,                        \
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge),       \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_ELASTICITY_BOUNDARY_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_ELASTICITY_BOUNDARY_TESTS(MPITest)
#endif
#undef RODIN_ELASTICITY_BOUNDARY_TESTS
#undef RODIN_ELASTICITY_BOUNDARY_EXTRA
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
    Rodin::Tests::Convergence::PETScLinearElasticityBoundaryTests::environment = &env;
    Rodin::Tests::Convergence::PETScLinearElasticityBoundaryTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
