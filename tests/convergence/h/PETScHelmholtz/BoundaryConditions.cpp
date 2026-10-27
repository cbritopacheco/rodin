/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Complex Helmholtz natural/impedance boundary convergence. */

#include "PETScHelmholtzProblem.h"
#include "FieldConvergence.h"
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
#include "CurvedGeometry.h"
#endif
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

#ifndef PETSC_USE_COMPLEX
#error "Complex Helmholtz boundary convergence requires complex PETSc"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::HelmholtzBoundaryTests
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n)
  {
    using Attributes = PETScHelmholtzProblem<1, LocalMesh>;
    const auto initialize = [](LocalMesh& mesh) {
      UnitBoxBoundary::labelCoordinatePartition(mesh, 0, Attributes::DirichletAttribute,
        Attributes::NaturalAttribute, Attributes::NaturalAttribute);
    };
    if constexpr (std::is_same_v<ContextType, Context::Local>)
    {
      auto mesh = UniformGrid(geometry).makeMesh(n);
      initialize(mesh);
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
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
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
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
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
      static constexpr size_t AffinePatchDegree = 2;
#else
      static constexpr size_t AffinePatchDegree = 1;
#endif
      template <size_t K>
      void checkRates(bool impedance) const
      {
        const auto levels =
          K == 1 ? std::array<size_t, 3>{5, 9, 17} : std::array<size_t, 3>{3, 5, 9};
        FieldConvergence<1> history;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          const auto mesh = makeMesh<ContextType>(this->GetParam(), n);
          using Problem = PETScHelmholtzProblem<K, std::remove_cvref_t<decltype(mesh)>>;
          const Problem problem(mesh, HelmholtzData::Field::Smooth, 16,
            impedance ? Problem::Boundary::Impedance : Problem::Boundary::MixedNeumann);
          history.append(Real(1) / Real(n - 1), {problem.solve()});
        }
        history.expectAlgebraicFloor({K == 1 ? Real(1.65) : Real(K) + Real(0.45),
          K == 1 ? Real(0.75) : Real(K) - Real(0.45)});
      }

      template <size_t K>
      void checkPatch(bool impedance, bool omitMass = false) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        using Problem = PETScHelmholtzProblem<K, std::remove_cvref_t<decltype(mesh)>>;
        const auto field =
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
          K == 1   ? HelmholtzData::Field::Constant
          : K == 2 ? HelmholtzData::Field::Affine
                   : HelmholtzData::Field::Quadratic;
#else
          K == 1 ? HelmholtzData::Field::Affine : HelmholtzData::Field::Quadratic;
#endif
        const Problem problem(mesh, field, 16,
          impedance ? Problem::Boundary::Impedance : Problem::Boundary::MixedNeumann);
        const auto error = problem.solve(omitMass);
        EXPECT_TRUE(error.isFinite());
        if (omitMass)
        {
          EXPECT_GT(error.getL2(), Real(1e-3));
          EXPECT_GT(error.getH1Seminorm(), Real(1e-2));
        }
        else
        {
          EXPECT_LT(error.getL2(), Real(1e-9));
          EXPECT_LT(error.getH1Seminorm(), Real(1e-9));
        }
      }

      void checkControl(bool impedance) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        using Problem = PETScHelmholtzProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const Problem problem(mesh,
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
          HelmholtzData::Field::Affine,
#else
          HelmholtzData::Field::Smooth,
#endif
          16, impedance ? Problem::Boundary::Impedance : Problem::Boundary::MixedNeumann);
        const auto correct = problem.solve();
        const auto wrong = problem.solve(false, Real(1e-13), 18, true);
        EXPECT_TRUE(correct.isFinite());
        EXPECT_TRUE(wrong.isFinite());
#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
        // The affine pullback is represented exactly in P2. This separates
        // omitted physical flux from coarse smooth-field discretization error.
        EXPECT_LT(correct.getL2(), Real(1e-9));
        EXPECT_LT(correct.getH1Seminorm(), Real(1e-9));
        EXPECT_GT(wrong.getL2(), Real(1e-3));
        EXPECT_GT(wrong.getH1Seminorm(), Real(1e-3));
#endif
        EXPECT_GT(wrong.getL2(), 5 * correct.getL2());
        EXPECT_GT(wrong.getH1Seminorm(), 5 * correct.getH1Seminorm());
      }

      void checkSensitivity(bool impedance) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        using Problem = PETScHelmholtzProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const auto boundary =
          impedance ? Problem::Boundary::Impedance : Problem::Boundary::MixedNeumann;
        const Problem base(mesh, HelmholtzData::Field::Smooth, 16, boundary);
        const Problem higher(mesh, HelmholtzData::Field::Smooth, 18, boundary);
        const auto reference = base.solve();
        for (const auto& varied : {higher.solve(), base.solve(false, Real(1e-13), 20),
               base.solve(false, Real(1e-14))})
        {
          EXPECT_LT(std::abs(varied.getL2() / reference.getL2() - 1), Real(1e-6));
          EXPECT_LT(
            std::abs(varied.getH1Seminorm() / reference.getH1Seminorm() - 1), Real(1e-6));
        }
      }
  };

#ifdef RODIN_HELMHOLTZ_BOUNDARY_CURVED
#define RODIN_HELMHOLTZ_BOUNDARY_PATCHES(Name, Prefix, Impedance)                        \
  TEST_P(Name, Prefix##ConstantP1Patch)                                                  \
  {                                                                                      \
    checkPatch<1>(Impedance);                                                            \
  }                                                                                      \
  TEST_P(Name, Prefix##AffineP2Patch)                                                    \
  {                                                                                      \
    checkPatch<2>(Impedance);                                                            \
  }                                                                                      \
  TEST_P(Name, Prefix##QuadraticP4Patch)                                                 \
  {                                                                                      \
    checkPatch<4>(Impedance);                                                            \
  }
#else
#define RODIN_HELMHOLTZ_BOUNDARY_PATCHES(Name, Prefix, Impedance)                        \
  TEST_P(Name, Prefix##AffineP1Patch)                                                    \
  {                                                                                      \
    checkPatch<1>(Impedance);                                                            \
  }                                                                                      \
  TEST_P(Name, Prefix##QuadraticP2Patch)                                                 \
  {                                                                                      \
    checkPatch<2>(Impedance);                                                            \
  }
#endif

#define RODIN_HELMHOLTZ_BOUNDARY_CONDITION(Name, Prefix, Impedance)                      \
  TEST_P(Name, Prefix##P1Rates)                                                          \
  {                                                                                      \
    checkRates<1>(Impedance);                                                            \
  }                                                                                      \
  TEST_P(Name, Prefix##P2Rates)                                                          \
  {                                                                                      \
    checkRates<2>(Impedance);                                                            \
  }                                                                                      \
  TEST_P(Name, Prefix##P3Rates)                                                          \
  {                                                                                      \
    checkRates<3>(Impedance);                                                            \
  }                                                                                      \
  RODIN_HELMHOLTZ_BOUNDARY_PATCHES(Name, Prefix, Impedance)                              \
  TEST_P(Name, Prefix##RejectsMissingMass)                                               \
  {                                                                                      \
    checkPatch<AffinePatchDegree>(Impedance, true);                                      \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsMissingFlux)                                               \
  {                                                                                      \
    checkControl(Impedance);                                                             \
  }                                                                                      \
  TEST_P(Name, Prefix##IndependentSensitivity)                                           \
  {                                                                                      \
    checkSensitivity(Impedance);                                                         \
  }

#define RODIN_HELMHOLTZ_BOUNDARY_TESTS(Name)                                             \
  RODIN_HELMHOLTZ_BOUNDARY_CONDITION(Name, MixedNeumann, false)                          \
  RODIN_HELMHOLTZ_BOUNDARY_CONDITION(Name, Impedance, true)                              \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,                 \
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,                        \
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge),       \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_HELMHOLTZ_BOUNDARY_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_HELMHOLTZ_BOUNDARY_TESTS(MPITest)
#endif
#undef RODIN_HELMHOLTZ_BOUNDARY_TESTS
#undef RODIN_HELMHOLTZ_BOUNDARY_CONDITION
#undef RODIN_HELMHOLTZ_BOUNDARY_PATCHES
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
    Rodin::Tests::Convergence::HelmholtzBoundaryTests::environment = &env;
    Rodin::Tests::Convergence::HelmholtzBoundaryTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
