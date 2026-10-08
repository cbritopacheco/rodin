/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Coupled reaction--diffusion natural-boundary certification. */

#include "../../PETScReactionDiffusionProblem.h"
#include "../../FieldConvergence.h"
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "../../MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::ReactionDiffusionBoundaryTests
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif
  using Attributes = PETScReactionDiffusionProblem<1, LocalMesh>;
  using Boundary = Attributes::Boundary;

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n)
  {
    const auto initialize = [](LocalMesh& mesh) {
      UnitBoxBoundary::labelCoordinatePartition(mesh, 0, Attributes::DirichletAttribute,
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
      return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
        .makeMesh(n, initialize);
#endif
  }

  template <class ContextType>
  class Fixture : public ::testing::TestWithParam<Polytope::Type>
  {
    public:
      template <size_t K>
      void checkRates(Boundary boundary) const
      {
        FieldConvergence<2> history;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                                   : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          const auto mesh = makeMesh<ContextType>(this->GetParam(), n);
          using Problem =
            PETScReactionDiffusionProblem<K, std::remove_cvref_t<decltype(mesh)>>;
          const Problem problem(mesh, ReactionDiffusionData::Field::Smooth, 16,
            static_cast<typename Problem::Boundary>(boundary));
          history.append(Real(1) / Real(n - 1), problem.solve(false, 1e-13, 18));
        }
        history.expectAlgebraicFloor({K == 1 ? Real(1.65) : Real(K) + Real(0.45),
          K == 1 ? Real(0.75) : Real(K) - Real(0.45)});
      }

      template <size_t K>
      void checkPatch(Boundary boundary, bool omitCoupling = false) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        using Problem =
          PETScReactionDiffusionProblem<K, std::remove_cvref_t<decltype(mesh)>>;
        const auto field = K == 1 ? ReactionDiffusionData::Field::Affine
                                  : ReactionDiffusionData::Field::Quadratic;
        const Problem problem(
          mesh, field, 16, static_cast<typename Problem::Boundary>(boundary));
        for (const auto& error : problem.solve(omitCoupling, 1e-13, 18))
        {
          EXPECT_TRUE(error.isFinite());
          if (omitCoupling)
          {
            EXPECT_GT(error.getL2(), Real(1e-3));
            EXPECT_GT(error.getH1Seminorm(), Real(1e-3));
          }
          else
          {
            EXPECT_LT(error.getL2(), Real(1e-9));
            EXPECT_LT(error.getH1Seminorm(), Real(1e-9));
          }
        }
      }

      void checkFluxControl(Boundary boundary) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        using Problem =
          PETScReactionDiffusionProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const Problem problem(mesh, ReactionDiffusionData::Field::Smooth, 16,
          static_cast<typename Problem::Boundary>(boundary));
        const auto correct = problem.solve(false, 1e-13, 18);
        const auto wrong = problem.solve(false, 1e-13, 18, true);
        for (size_t component = 0; component < 2; ++component)
        {
          SCOPED_TRACE(::testing::Message() << "component=" << component);
          EXPECT_TRUE(wrong[component].isFinite());
          EXPECT_GT(wrong[component].getL2(), 5 * correct[component].getL2());
          EXPECT_GT(
            wrong[component].getH1Seminorm(), 5 * correct[component].getH1Seminorm());
        }
      }

      void checkSensitivity(Boundary boundary) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        using Problem =
          PETScReactionDiffusionProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const auto kind = static_cast<typename Problem::Boundary>(boundary);
        const Problem base(mesh, ReactionDiffusionData::Field::Smooth, 16, kind);
        const Problem higher(mesh, ReactionDiffusionData::Field::Smooth, 18, kind);
        const auto reference = base.solve(false, 1e-13, 18);
        for (const auto& varied : {higher.solve(false, 1e-13, 18),
               base.solve(false, 1e-13, 20), base.solve(false, 1e-14, 18)})
        {
          for (size_t component = 0; component < 2; ++component)
          {
            SCOPED_TRACE(::testing::Message() << "component=" << component);
            EXPECT_GT(reference[component].getL2(), 0);
            EXPECT_GT(reference[component].getH1Seminorm(), 0);
            EXPECT_LT(
              std::abs(varied[component].getL2() / reference[component].getL2() - 1),
              Real(1e-6));
            EXPECT_LT(std::abs(varied[component].getH1Seminorm() /
                          reference[component].getH1Seminorm() -
                        1),
              Real(1e-6));
          }
        }
      }
  };

#define RODIN_REACTION_BOUNDARY_CASES(Name, Prefix, Kind)                                \
  TEST_P(Name, Prefix##P1Rates)                                                          \
  {                                                                                      \
    checkRates<1>(Boundary::Kind);                                                       \
  }                                                                                      \
  TEST_P(Name, Prefix##P2Rates)                                                          \
  {                                                                                      \
    checkRates<2>(Boundary::Kind);                                                       \
  }                                                                                      \
  TEST_P(Name, Prefix##P3Rates)                                                          \
  {                                                                                      \
    checkRates<3>(Boundary::Kind);                                                       \
  }                                                                                      \
  TEST_P(Name, Prefix##AffineP1Patch)                                                    \
  {                                                                                      \
    checkPatch<1>(Boundary::Kind);                                                       \
  }                                                                                      \
  TEST_P(Name, Prefix##QuadraticP2Patch)                                                 \
  {                                                                                      \
    checkPatch<2>(Boundary::Kind);                                                       \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsMissingCoupling)                                           \
  {                                                                                      \
    checkPatch<1>(Boundary::Kind, true);                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsMissingFlux)                                               \
  {                                                                                      \
    checkFluxControl(Boundary::Kind);                                                    \
  }                                                                                      \
  TEST_P(Name, Prefix##IndependentSensitivity)                                           \
  {                                                                                      \
    checkSensitivity(Boundary::Kind);                                                    \
  }

#define RODIN_REACTION_BOUNDARY_TESTS(Name)                                              \
  RODIN_REACTION_BOUNDARY_CASES(Name, MixedNeumann, MixedNeumann)                        \
  RODIN_REACTION_BOUNDARY_CASES(Name, Robin, Robin)                                      \
  RODIN_REACTION_BOUNDARY_CASES(Name, PureNeumann, PureNeumann)                          \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,                 \
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,                        \
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge),       \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_REACTION_BOUNDARY_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_REACTION_BOUNDARY_TESTS(MPITest)
#endif
#undef RODIN_REACTION_BOUNDARY_TESTS
#undef RODIN_REACTION_BOUNDARY_CASES
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
    Rodin::Tests::Convergence::ReactionDiffusionBoundaryTests::environment = &env;
    Rodin::Tests::Convergence::ReactionDiffusionBoundaryTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
