/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Semilinear natural-boundary rates and Newton tangent oracles. */

#include "../../PETScNonlinearPoisson.h"
#include "../../FieldConvergence.h"
#ifdef RODIN_USE_MPI
#include <boost/mpi/environment.hpp>
#include "../../MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Convergence::NonlinearBoundaryTests
{
#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif
  using Attributes = PETScNonlinearPoissonProblem<1, LocalMesh>;
  using Boundary = Attributes::Boundary;
  using Field = NonlinearPoissonData::Field;

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

  /** A prvalue factory preserves the workload's internal field/solver references. */
  template <size_t K, class MeshType>
  auto makeProblem(const MeshType& mesh, Boundary boundary, Field field = Field::Sine,
    size_t assemblyOrder = 16, bool omitCubic = false, bool wrongTangent = false,
    bool omitFlux = false, bool wrongBoundaryTangent = false)
  {
    using Problem = PETScNonlinearPoissonProblem<K, MeshType>;
    return Problem(mesh, assemblyOrder, field == Field::Sine ? Real(1) : Real(0.25),
      omitCubic, wrongTangent, true, field,
      static_cast<typename Problem::Boundary>(boundary), omitFlux, wrongBoundaryTangent);
  }

  template <class Problem>
  auto measure(Problem& problem, Real tolerance = 1e-12, size_t normOrder = 18)
  {
    return problem.solve(tolerance, normOrder, [](const auto&, const auto&) {});
  }

  template <class ContextType>
  class Fixture : public ::testing::TestWithParam<Polytope::Type>
  {
    public:
      template <size_t K>
      void checkRates(Boundary boundary) const
      {
        FieldConvergence<1> history;
        const bool resolveNeumannTransient = K == 1 &&
          boundary == Boundary::PureNeumann &&
          this->GetParam() == Polytope::Type::Tetrahedron;
        const auto levels = resolveNeumannTransient
          ? std::initializer_list<size_t>{9, 17, 33}
          : K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                   : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          const auto mesh = makeMesh<ContextType>(this->GetParam(), n);
          auto problem = makeProblem<K>(mesh, boundary);
          history.append(Real(1) / Real(n - 1), {measure(problem)});
        }
        history.expectAlgebraicFloor({K == 1 ? Real(1.65) : Real(K) + Real(0.45),
          K == 1 ? Real(0.75) : Real(K) - Real(0.45)});
      }

      template <size_t K>
      void checkPatch(Boundary boundary, Field field) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        auto problem = makeProblem<K>(mesh, boundary, field);
        const auto error = measure(problem);
        EXPECT_TRUE(error.isFinite());
        EXPECT_LT(error.getL2(), Real(1e-9));
        EXPECT_LT(error.getH1Seminorm(), Real(1e-9));
      }

      void checkCubicControl(Boundary boundary) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        auto problem = makeProblem<2>(mesh, boundary, Field::Affine, 16, true);
        const auto error = measure(problem);
        EXPECT_TRUE(error.isFinite());
        EXPECT_GT(error.getL2(), Real(1e-3));
        EXPECT_GT(error.getH1Seminorm(), Real(1e-3));
      }

      void checkFluxControl(Boundary boundary) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        // An affine patch isolates the missing boundary term from the
        // coarse-mesh approximation error of a non-polynomial solution.
        auto correctProblem = makeProblem<2>(mesh, boundary, Field::Affine);
        auto wrongProblem =
          makeProblem<2>(mesh, boundary, Field::Affine, 16, false, false, true);
        const auto correct = measure(correctProblem), wrong = measure(wrongProblem);
        EXPECT_TRUE(correct.isFinite());
        EXPECT_LT(correct.getL2(), Real(1e-9));
        EXPECT_LT(correct.getH1Seminorm(), Real(1e-9));
        EXPECT_TRUE(wrong.isFinite());
        EXPECT_GT(wrong.getL2(), Real(1e-3));
        EXPECT_GT(wrong.getH1Seminorm(), Real(1e-3));
        EXPECT_GT(wrong.getL2(), 5 * correct.getL2());
        EXPECT_GT(wrong.getH1Seminorm(), 5 * correct.getH1Seminorm());
      }

      template <size_t K>
      void checkTangent(
        Boundary boundary, bool wrongCubic = false, bool wrongRobin = false) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        auto problem = makeProblem<K>(
          mesh, boundary, Field::Sine, 16, false, wrongCubic, false, wrongRobin);
        // x_0 vanishes on Gamma_D but not on the Robin faces: a sine
        // direction would have zero trace and could miss a boundary tangent bug.
        const Variational::RealFunction direction(
          [](const Point& point) { return point.getPhysicalCoordinates()(0); });
        const Real defect = problem.tangentDefect(direction);
        EXPECT_TRUE(std::isfinite(defect));
        if (wrongCubic || wrongRobin)
          EXPECT_GT(defect, Real(1e-3));
        else
          EXPECT_LT(defect, Real(1e-6));
      }

      void checkSensitivity(Boundary boundary) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3);
        auto base = makeProblem<2>(mesh, boundary);
        const auto reference = measure(base);
        EXPECT_GT(reference.getL2(), 0);
        EXPECT_GT(reference.getH1Seminorm(), 0);
        for (size_t variation = 0; variation < 3; ++variation)
        {
          auto variedProblem =
            makeProblem<2>(mesh, boundary, Field::Sine, variation == 0 ? 18 : 16);
          const auto varied = measure(variedProblem,
            variation == 2 ? Real(1e-13) : Real(1e-12), variation == 1 ? 20 : 18);
          EXPECT_LT(std::abs(varied.getL2() / reference.getL2() - 1), Real(1e-6));
          EXPECT_LT(
            std::abs(varied.getH1Seminorm() / reference.getH1Seminorm() - 1), Real(1e-6));
        }
      }
  };

#define RODIN_NONLINEAR_BOUNDARY_CASES(Name, Prefix, Kind)                               \
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
  TEST_P(Name, Prefix##ConstantP1Patch)                                                  \
  {                                                                                      \
    checkPatch<1>(Boundary::Kind, Field::Constant);                                      \
  }                                                                                      \
  TEST_P(Name, Prefix##AffineP1Patch)                                                    \
  {                                                                                      \
    checkPatch<1>(Boundary::Kind, Field::Affine);                                        \
  }                                                                                      \
  TEST_P(Name, Prefix##QuadraticP2Patch)                                                 \
  {                                                                                      \
    checkPatch<2>(Boundary::Kind, Field::Quadratic);                                     \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsMissingCubic)                                              \
  {                                                                                      \
    checkCubicControl(Boundary::Kind);                                                   \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsMissingFlux)                                               \
  {                                                                                      \
    checkFluxControl(Boundary::Kind);                                                    \
  }                                                                                      \
  TEST_P(Name, Prefix##P1ResidualTangent)                                                \
  {                                                                                      \
    checkTangent<1>(Boundary::Kind);                                                     \
  }                                                                                      \
  TEST_P(Name, Prefix##P2ResidualTangent)                                                \
  {                                                                                      \
    checkTangent<2>(Boundary::Kind);                                                     \
  }                                                                                      \
  TEST_P(Name, Prefix##RejectsWrongCubicTangent)                                         \
  {                                                                                      \
    checkTangent<2>(Boundary::Kind, true);                                               \
  }                                                                                      \
  TEST_P(Name, Prefix##IndependentSensitivity)                                           \
  {                                                                                      \
    checkSensitivity(Boundary::Kind);                                                    \
  }

#define RODIN_NONLINEAR_BOUNDARY_TESTS(Name)                                             \
  RODIN_NONLINEAR_BOUNDARY_CASES(Name, MixedNeumann, MixedNeumann)                       \
  RODIN_NONLINEAR_BOUNDARY_CASES(Name, Robin, Robin)                                     \
  RODIN_NONLINEAR_BOUNDARY_CASES(Name, PureNeumann, PureNeumann)                         \
  TEST_P(Name, RejectsMissingRobinTangent)                                               \
  {                                                                                      \
    checkTangent<2>(Boundary::Robin, false, true);                                       \
  }                                                                                      \
  INSTANTIATE_TEST_SUITE_P(AllGeometries, Name,                                          \
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,                 \
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,                        \
      Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge),       \
    [](const auto& info) {                                                               \
      return std::string(UniformGrid::getGeometryName(info.param));                      \
    });

  using LocalTest = Fixture<Context::Local>;
  RODIN_NONLINEAR_BOUNDARY_TESTS(LocalTest)
#ifdef RODIN_USE_MPI
  using MPITest = Fixture<Context::MPI>;
  RODIN_NONLINEAR_BOUNDARY_TESTS(MPITest)
#endif
#undef RODIN_NONLINEAR_BOUNDARY_TESTS
#undef RODIN_NONLINEAR_BOUNDARY_CASES
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
    Rodin::Tests::Convergence::NonlinearBoundaryTests::environment = &env;
    Rodin::Tests::Convergence::NonlinearBoundaryTests::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  return PetscFinalize() == PETSC_SUCCESS ? result : 1;
}
