/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/// @file @brief Physical-traction Stokes refinement and pressure-level checks.

#include "../../PETScStokesTractionProblem.h"
#include "../../FieldConvergence.h"
#ifdef RODIN_STOKES_TRACTION_CURVED
#include "../../CurvedGeometry.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::StokesTractionTests
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  constexpr Real PressureOffset = 2;
  constexpr Real PatchBudget = 1e-9;
  constexpr Real SensitivityBudget = 1e-6;

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n)
  {
    using Attributes = PETScStokesTractionProblem<LocalMesh>;
    const auto initialize = [](LocalMesh& mesh) {
      UnitBoxBoundary::labelCoordinatePartition(mesh, 0, Attributes::DirichletAttribute,
        Attributes::TractionAttribute, Attributes::TractionAttribute);
    };
    if constexpr (std::is_same_v<ContextType, Context::Local>)
    {
      auto mesh = UniformGrid(geometry).makeMesh(n);
      initialize(mesh);
#ifdef RODIN_STOKES_TRACTION_CURVED
      CurvedGeometry curved(mesh);
      curved.template install<2>();
#endif
      return mesh;
    }
#ifdef RODIN_USE_MPI
    else
    {
      auto mesh = DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
                    .makeMesh(n, initialize);
#ifdef RODIN_STOKES_TRACTION_CURVED
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
#ifdef RODIN_STOKES_TRACTION_CURVED
      static constexpr size_t QuadraticPatchDegree = 4;
#else
      static constexpr size_t QuadraticPatchDegree = 2;
#endif
      template <size_t K>
      void checkRates() const
      {
        FieldConvergence<1> velocity, pressure;
        for (size_t n : {3, 5, 9})
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          const auto mesh = makeMesh<ContextType>(this->GetParam(), n);
          const StokesData data(mesh.getDimension(),
            K == 2 ? StokesData::Field::Cubic : StokesData::Field::Quartic, 1,
            PressureOffset);
          const auto result = PETScStokesTractionProblem(mesh, data).template solve<K>();
          velocity.append(Real(1) / Real(n - 1), {result.fields.velocity});
          pressure.append(Real(1) / Real(n - 1), {result.fields.pressure});
          EXPECT_LE(std::abs(result.pressureMean - PressureOffset),
            result.fields.pressure.getL2() + PatchBudget);
        }
        velocity.expectAlgebraicFloor({Real(K) + Real(0.45), Real(K) - Real(0.45)});
        pressure.expectAlgebraicFloor({Real(K) - Real(0.55), Real(K) - Real(1.45)});
      }

      template <size_t K>
      void checkPatch(StokesData::Field field) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        const StokesData data(mesh.getDimension(), field, 1, PressureOffset);
        const auto result = PETScStokesTractionProblem(mesh, data).template solve<K>();
        expectPatchFields(result);
        EXPECT_LT(result.fields.pressure.getL2(), PatchBudget);
        EXPECT_NEAR(result.pressureMean, PressureOffset, PatchBudget);
      }

      void checkPressureLevelControl() const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        const StokesData data(
          mesh.getDimension(), StokesData::Field::Quadratic, 1, PressureOffset);
        const auto result = PETScStokesTractionProblem(mesh, data)
                              .template solve<QuadraticPatchDegree>(18, PressureOffset);
        // Adding c*n to traction makes p_h=p_exact-c, not a velocity error.
        expectPatchFields(result);
        EXPECT_NEAR(result.fields.pressure.getL2(), PressureOffset, PatchBudget);
        EXPECT_NEAR(result.pressureMean, 0, PatchBudget);
      }

      void checkSensitivity() const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        const StokesData data(
          mesh.getDimension(), StokesData::Field::Cubic, 1, PressureOffset);
        const PETScStokesTractionProblem base(mesh, data, 16), higher(mesh, data, 18);
        const auto reference = base.template solve<2>();
        for (const auto& varied :
          {higher.template solve<2>(), base.template solve<2>(20)})
        {
          for (const auto& pair :
            {std::pair{reference.fields.velocity, varied.fields.velocity},
              std::pair{reference.fields.pressure, varied.fields.pressure}})
          {
            ASSERT_GT(pair.first.getL2(), 0);
            ASSERT_GT(pair.first.getH1Seminorm(), 0);
            EXPECT_LT(
              std::abs(pair.second.getL2() / pair.first.getL2() - 1), SensitivityBudget);
            EXPECT_LT(
              std::abs(pair.second.getH1Seminorm() / pair.first.getH1Seminorm() - 1),
              SensitivityBudget);
          }
        }
      }

    private:
      void expectPatchFields(const StokesTractionErrors& result) const
      {
        EXPECT_LT(result.fields.velocity.getL2(), PatchBudget);
        EXPECT_LT(result.fields.velocity.getH1Seminorm(), PatchBudget);
        EXPECT_LT(result.fields.pressure.getH1Seminorm(), PatchBudget);
        EXPECT_LT(result.fields.divergence, PatchBudget);
      }
  };

#ifdef RODIN_STOKES_TRACTION_CURVED
#define RODIN_STOKES_TRACTION_PATCHES(Name)                                              \
  TEST_P(Name, QuadraticP4P3Patch)                                                       \
  {                                                                                      \
    checkPatch<4>(StokesData::Field::Quadratic);                                         \
  }                                                                                      \
  TEST_P(Name, CubicP6P5Patch)                                                           \
  {                                                                                      \
    checkPatch<6>(StokesData::Field::Cubic);                                             \
  }
#else
#define RODIN_STOKES_TRACTION_PATCHES(Name)                                              \
  TEST_P(Name, QuadraticP2P1Patch)                                                       \
  {                                                                                      \
    checkPatch<2>(StokesData::Field::Quadratic);                                         \
  }                                                                                      \
  TEST_P(Name, CubicP3P2Patch)                                                           \
  {                                                                                      \
    checkPatch<3>(StokesData::Field::Cubic);                                             \
  }
#endif

#define RODIN_STOKES_TRACTION_TESTS(Name)                                                \
  TEST_P(Name, P2P1Rates)                                                                \
  {                                                                                      \
    checkRates<2>();                                                                     \
  }                                                                                      \
  TEST_P(Name, P3P2Rates)                                                                \
  {                                                                                      \
    checkRates<3>();                                                                     \
  }                                                                                      \
  RODIN_STOKES_TRACTION_PATCHES(Name)                                                    \
  TEST_P(Name, TractionDeterminesPressureLevel)                                          \
  {                                                                                      \
    checkPressureLevelControl();                                                         \
  }                                                                                      \
  TEST_P(Name, IndependentQuadratureSensitivity)                                         \
  {                                                                                      \
    checkSensitivity();                                                                  \
  }                                                                                      \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,           \
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,  \
      Polytope::Type::Wedge),                                                            \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_STOKES_TRACTION_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_STOKES_TRACTION_TESTS(MPITest)
#endif
#undef RODIN_STOKES_TRACTION_TESTS
#undef RODIN_STOKES_TRACTION_PATCHES
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
    Rodin::Tests::Convergence::StokesTractionTests::environment = &env;
    Rodin::Tests::Convergence::StokesTractionTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
