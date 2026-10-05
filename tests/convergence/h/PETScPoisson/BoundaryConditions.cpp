/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI Poisson boundary h-convergence. */

#include "PETScPoissonBoundaryProblem.h"
#include "FieldConvergence.h"
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#include "Rodin/MPI/Variational/P0g/P0g.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::PETScPoissonBoundaryTests
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n, bool pure)
  {
    using Attributes = PETScPoissonBoundaryProblem<1, LocalMesh>;
    const auto initialize = [pure](LocalMesh& mesh) {
      UnitBoxBoundary::labelCoordinatePartition(mesh, 0,
        pure ? Attributes::NaturalAttribute : Attributes::DirichletAttribute,
        Attributes::NaturalAttribute, Attributes::NaturalAttribute);
    };
    if constexpr (std::is_same_v<ContextType, Context::Local>)
    {
      auto mesh = UniformGrid(geometry).makeMesh(n);
      initialize(mesh);
      return mesh;
    }
#ifdef RODIN_USE_MPI
    else
    {
      Context::MPI context(*environment, *world);
      return DistributedUniformGrid(context, geometry).makeMesh(n, initialize);
    }
#endif
  }

  template <class ContextType>
  class Fixture : public ::testing::TestWithParam<Polytope::Type>
  {
    public:
      template <size_t K>
      void checkRates(int condition) const
      {
        const std::array<size_t, 3> levels =
          K == 1 ? std::array<size_t, 3>{5, 9, 17} : std::array<size_t, 3>{3, 5, 9};
        FieldConvergence<1> history;
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          const auto mesh = makeMesh<ContextType>(this->GetParam(), n, condition == 2);
          using Problem =
            PETScPoissonBoundaryProblem<K, std::remove_cvref_t<decltype(mesh)>>;
          const Problem problem(
            mesh, static_cast<typename Problem::Condition>(condition));
          history.append(Real(1) / Real(n - 1), {problem.solve()});
        }
        history.expectAlgebraicFloor({K == 1 ? Real(1.65) : Real(K) + Real(0.45),
          K == 1     ? Real(0.75)
            : K == 2 ? Real(1.55)
                     : Real(2.35)});
      }

      void checkControl(int condition) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3, condition == 2);
        using Problem =
          PETScPoissonBoundaryProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const Problem problem(mesh, static_cast<typename Problem::Condition>(condition));
        const auto correct = problem.solve();
        const auto wrong = problem.solve(true);
        EXPECT_TRUE(correct.isFinite());
        EXPECT_TRUE(wrong.isFinite());
        EXPECT_GT(wrong.getL2(), WrongErrorFactor * correct.getL2());
        EXPECT_GT(wrong.getH1Seminorm(), WrongErrorFactor * correct.getH1Seminorm());
      }

      void checkSensitivity(int condition) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3, condition == 2);
        using Problem =
          PETScPoissonBoundaryProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const auto kind = static_cast<typename Problem::Condition>(condition);
        const Problem base(mesh, kind);
        const Problem higherAssembly(mesh, kind, 18);
        const auto reference = base.solve();
        for (const auto& varied :
          {higherAssembly.solve(), base.solve(false, 20), base.solve(false, 18, 1e-14)})
        {
          EXPECT_LT(std::abs(varied.getL2() - reference.getL2()) / reference.getL2(),
            SensitivityTolerance);
          EXPECT_LT(std::abs(varied.getH1Seminorm() - reference.getH1Seminorm()) /
              reference.getH1Seminorm(),
            SensitivityTolerance);
        }
      }

    private:
      // Missing flux must dominate the resolved discretization error.
      static constexpr Real WrongErrorFactor = 5;
      static constexpr Real SensitivityTolerance = 1e-6;
  };

#define RODIN_POISSON_BOUNDARY_CONDITION(Name, Prefix, Kind)                             \
  TEST_P(Name, Prefix##P1Rates)                                                          \
  {                                                                                      \
    checkRates<1>(Kind);                                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##P2Rates)                                                          \
  {                                                                                      \
    checkRates<2>(Kind);                                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##P3Rates)                                                          \
  {                                                                                      \
    checkRates<3>(Kind);                                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsMissingFlux)                                               \
  {                                                                                      \
    checkControl(Kind);                                                                  \
  }                                                                                      \
  TEST_P(Name, Prefix##IndependentSensitivity)                                           \
  {                                                                                      \
    checkSensitivity(Kind);                                                              \
  }

#ifdef RODIN_POISSON_BOUNDARY_MUMPS
#define RODIN_POISSON_BOUNDARY_PURE(Name)                                                \
  RODIN_POISSON_BOUNDARY_CONDITION(Name, PureNeumann, 2)
#else
#define RODIN_POISSON_BOUNDARY_PURE(Name)
#endif

#define RODIN_POISSON_BOUNDARY_TESTS(Name)                                               \
  RODIN_POISSON_BOUNDARY_CONDITION(Name, MixedNeumann, 0)                                \
  RODIN_POISSON_BOUNDARY_CONDITION(Name, MixedRobin, 1)                                  \
  RODIN_POISSON_BOUNDARY_PURE(Name)                                                      \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,                 \
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,                        \
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge),       \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_POISSON_BOUNDARY_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_POISSON_BOUNDARY_TESTS(MPITest)
#endif
#undef RODIN_POISSON_BOUNDARY_TESTS
#undef RODIN_POISSON_BOUNDARY_PURE
#undef RODIN_POISSON_BOUNDARY_CONDITION
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
    Rodin::Tests::Convergence::PETScPoissonBoundaryTests::environment = &env;
    Rodin::Tests::Convergence::PETScPoissonBoundaryTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
