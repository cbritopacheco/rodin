/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI combined mesh/degree Stokes refinement. */

#include "../../PETScStokesProblem.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::HP::PETScStokes
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  class PETScStokesHPTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      auto makeMesh(size_t n) const
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(GetParam()).makeMesh(n);
#ifdef RODIN_USE_MPI
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), GetParam())
            .makeMesh(n);
#endif
      }

      void checkPath() const
      {
        const StokesData data(
          UniformGrid(GetParam()).getDimension(), StokesData::Field::Smooth);
        ErrorHistory velocity, pressure;
        const auto append = [&](Real h, const StokesErrors& errors) {
          velocity.append(h, errors.velocity);
          pressure.append(h, errors.pressure);
        };
        append(0.5, PETScStokesProblem(makeMesh(3), data).template solve<2>());
        append(Real(1) / 3, PETScStokesProblem(makeMesh(4), data).template solve<3>());
        append(0.25, PETScStokesProblem(makeMesh(5), data).template solve<4>());
        for (bool isVelocity : {true, false})
        {
          const auto& history = isVelocity ? velocity : pressure;
          ASSERT_EQ(history.getSize(), 3u);
          for (size_t i = 1; i < history.getSize(); ++i)
          {
            const auto& coarse = history.getSample(i - 1).error;
            const auto& fine = history.getSample(i).error;
            SCOPED_TRACE(::testing::Message()
              << "velocity=" << isVelocity << " interval=" << i << " L2 "
              << coarse.getL2() << " -> " << fine.getL2() << " H1 "
              << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm());
            ASSERT_TRUE(coarse.isFinite());
            ASSERT_TRUE(fine.isFinite());
            ASSERT_GT(fine.getL2(), 0);
            ASSERT_GT(fine.getH1Seminorm(), 0);
            ASSERT_GT(coarse.getL2(), fine.getL2());
            ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            const auto rate = history.getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "effective rates " << rate.getL2() << ", " << rate.getH1Seminorm());
            EXPECT_GT(rate.getL2(), isVelocity ? 2.5 : 1.5);
            EXPECT_GT(rate.getH1Seminorm(), isVelocity ? 1.5 : 0.5);
          }
        }
      }

      void checkPatch() const
      {
        auto mesh = makeMesh(3);
        const auto errors = PETScStokesProblem(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Quartic))
                              .template solve<4>();
        EXPECT_LT(errors.velocity.getL2(), 1e-9);
        EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
        EXPECT_LT(errors.pressure.getL2(), 1e-9);
        EXPECT_LT(errors.pressure.getH1Seminorm(), 1e-9);
        EXPECT_LT(errors.divergence, 1e-9);
      }

      void checkWrongViscosity() const
      {
        auto mesh = makeMesh(3);
        const auto errors = PETScStokesProblem(
          mesh, StokesData(mesh.getDimension(), StokesData::Field::Quadratic))
                              .template solve<4>(2);
        EXPECT_LT(errors.velocity.getL2(), 1e-9);
        EXPECT_LT(errors.velocity.getH1Seminorm(), 1e-9);
        EXPECT_GT(errors.pressure.getL2(), 0.1);
        EXPECT_GT(errors.pressure.getH1Seminorm(), 1);
      }

      void checkQuadrature() const
      {
        auto mesh = makeMesh(3);
        const StokesData data(mesh.getDimension(), StokesData::Field::Smooth);
        const auto a = PETScStokesProblem(mesh, data).template solve<4>();
        const auto b = PETScStokesProblem(mesh, data, 18).template solve<4>();
        for (const auto& pair :
          {std::pair{a.velocity, b.velocity}, std::pair{a.pressure, b.pressure}})
        {
          ASSERT_TRUE(pair.first.isFinite());
          ASSERT_TRUE(pair.second.isFinite());
          ASSERT_GT(pair.first.getL2(), 0);
          ASSERT_GT(pair.first.getH1Seminorm(), 0);
          EXPECT_LT(std::abs(pair.second.getL2() / pair.first.getL2() - 1), 1e-6);
          EXPECT_LT(
            std::abs(pair.second.getH1Seminorm() / pair.first.getH1Seminorm() - 1), 1e-6);
        }
      }
  };

  using PETScStokesHPLocalTest = PETScStokesHPTest<Context::Local>;
  TEST_P(PETScStokesHPLocalTest, BothFieldsImproveAlongCombinedPath)
  {
    checkPath();
  }
  TEST_P(PETScStokesHPLocalTest, P4P3ReproducesQuarticPatch)
  {
    checkPatch();
  }
  TEST_P(PETScStokesHPLocalTest, P4P3RejectsIncorrectViscosity)
  {
    checkWrongViscosity();
  }
  TEST_P(PETScStokesHPLocalTest, P4QuadratureSensitivity)
  {
    checkQuadrature();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScStokesHPLocalTest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#ifdef RODIN_USE_MPI
  using PETScStokesHPMPITest = PETScStokesHPTest<Context::MPI>;
  TEST_P(PETScStokesHPMPITest, BothFieldsImproveAlongCombinedPath)
  {
    checkPath();
  }
  TEST_P(PETScStokesHPMPITest, P4P3ReproducesQuarticPatch)
  {
    checkPatch();
  }
  TEST_P(PETScStokesHPMPITest, P4P3RejectsIncorrectViscosity)
  {
    checkWrongViscosity();
  }
  TEST_P(PETScStokesHPMPITest, P4QuadratureSensitivity)
  {
    checkQuadrature();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, PETScStokesHPMPITest,
    ::testing::Values(Polytope::Type::Triangle, Polytope::Type::Quadrilateral,
      Polytope::Type::Tetrahedron, Polytope::Type::Pyramid, Polytope::Type::Hexahedron,
      Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#endif
}

int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  // Centralized RHS/solution matches the existing PETSc Stokes h workload.
  for (const char* name : {"-mat_mumps_icntl_20", "-mat_mumps_icntl_21"})
  {
    PetscBool set = PETSC_FALSE;
    if (PetscOptionsHasName(nullptr, nullptr, name, &set) != PETSC_SUCCESS)
      return 1;
    if (!set && PetscOptionsSetValue(nullptr, name, "0") != PETSC_SUCCESS)
      return 1;
  }
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::HP::PETScStokes::environment = &env;
    Rodin::Tests::Convergence::HP::PETScStokes::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
