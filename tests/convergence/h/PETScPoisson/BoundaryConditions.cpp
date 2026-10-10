/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI scalar diffusion boundary h-convergence. */

#include "PETScPoissonBoundaryProblem.h"
#include "FieldConvergence.h"
#ifdef RODIN_DIFFUSION_BOUNDARY_CURVED
#include "CurvedGeometry.h"
#endif
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
  // Both coefficients use the same boundary, norm and solver-budget oracles.
  // Compile-time selection preserves the unweighted Poisson assembly path.
  template <size_t K, class MeshType>
  using BoundaryProblem = PETScDiffusionBoundaryProblem<K, MeshType,
#ifdef RODIN_DIFFUSION_BOUNDARY_VARIABLE
    false
#else
    true
#endif
    >;

#ifdef RODIN_USE_MPI
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  // Reference means are data, not a distributed query. For the exact
  // quadratic map the d>=2 determinant is one; only a one-dimensional
  // exponential integral needs quadrature, independently of the PDE mesh.
  Optional<Real> referenceMean(
    size_t dim, ConductivityData::Field field = ConductivityData::Field::Exponential)
  {
#ifdef RODIN_DIFFUSION_BOUNDARY_CURVED
    constexpr Real amplitude = 0.1;
    // Independent reference-data budgets, not PDE discretization tolerances.
    constexpr size_t ReferenceOrder = 24, CheckOrder = 32;
    constexpr Real ReferenceTolerance = 1e-13;
    // At amplitude 0.1 these truncations bound the positive-series tail
    // below 6e-41, far below the floating-point comparison budget.
    constexpr size_t LinearTerms = 40, QuadraticTerms = 20;
    const Real length = 1 + amplitude;
    switch (field)
    {
      case ConductivityData::Field::Constant:
        return 1;
      case ConductivityData::Field::Affine:
        return dim == 1 ? 1 + length / 2 : 1 + Real(dim) / 2 + amplitude / 3;
      case ConductivityData::Field::Quadratic:
        return dim == 1 ? 1 + length * length / 3
                        : 1 + Real(dim) / 3 + amplitude / 3 + amplitude * amplitude / 5;
      case ConductivityData::Field::Exponential:
      {
        if (dim == 1)
          return std::expm1(length) / length;
        const auto integrate = [](size_t order) {
          const auto& rule =
            QF::PolytopeQuadratureFormula::get(order, Polytope::Type::Segment);
          Real integral = 0;
          for (size_t i = 0; i < rule.getSize(); ++i)
          {
            const Real t = rule.getPoint(i)(0);
            integral += rule.getWeight(i) * std::exp(t + amplitude * t * t);
          }
          return integral;
        };
        const Real value = integrate(ReferenceOrder);
        EXPECT_NEAR(value, integrate(CheckOrder), ReferenceTolerance);
        // Independent positive series for integral_0^1 exp(t+a*t^2) dt:
        // sum_{i,j>=0} a^j / (i! j! (i+2*j+1)). No mapped quadrature,
        // mesh metric or cancellation-prone moment recurrence is reused.
        Real series = 0, quadraticTerm = 1;
        for (size_t j = 0; j <= QuadraticTerms; ++j)
        {
          Real linearTerm = 1;
          for (size_t i = 0; i <= LinearTerms; ++i)
          {
            series += quadraticTerm * linearTerm / Real(i + 2 * j + 1);
            linearTerm /= Real(i + 1);
          }
          quadraticTerm *= amplitude / Real(j + 1);
        }
        EXPECT_NEAR(value, series, ReferenceTolerance);
        return value * std::pow(std::expm1(Real(1)), Real(dim - 1));
      }
      case ConductivityData::Field::Smooth:
        break; // This suite uses exponential rate fields and polynomial patches.
    }
    ADD_FAILURE() << "No curved-domain reference mean for the selected field";
    return 0;
#else
    (void)dim;
    (void)field;
    return std::nullopt;
#endif
  }

  template <class ContextType>
  auto makeMesh(Polytope::Type geometry, size_t n, bool pure)
  {
    using Attributes = BoundaryProblem<1, LocalMesh>;
    const auto initialize = [pure](LocalMesh& mesh) {
      UnitBoxBoundary::labelCoordinatePartition(mesh, 0,
        pure ? Attributes::NaturalAttribute : Attributes::DirichletAttribute,
        Attributes::NaturalAttribute, Attributes::NaturalAttribute);
    };
    if constexpr (std::is_same_v<ContextType, Context::Local>)
    {
      auto mesh = UniformGrid(geometry).makeMesh(n);
      initialize(mesh);
#ifdef RODIN_DIFFUSION_BOUNDARY_CURVED
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
#ifdef RODIN_DIFFUSION_BOUNDARY_CURVED
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
          using Problem = BoundaryProblem<K, std::remove_cvref_t<decltype(mesh)>>;
          const Problem problem(mesh, static_cast<typename Problem::Condition>(condition),
            16, ConductivityData::Field::Exponential, referenceMean(mesh.getDimension()));
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
        using Problem = BoundaryProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const Problem problem(mesh, static_cast<typename Problem::Condition>(condition),
          16, ConductivityData::Field::Exponential, referenceMean(mesh.getDimension()));
        const auto correct = problem.solve();
        const auto wrong = problem.solve(true);
        EXPECT_TRUE(correct.isFinite());
        EXPECT_TRUE(wrong.isFinite());
        EXPECT_GT(wrong.getL2(), WrongErrorFactor * correct.getL2());
        EXPECT_GT(wrong.getH1Seminorm(), WrongErrorFactor * correct.getH1Seminorm());
      }

      template <size_t K>
      void checkPatch(int condition) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3, condition == 2);
        using Problem = BoundaryProblem<K, std::remove_cvref_t<decltype(mesh)>>;
        constexpr auto field =
#ifdef RODIN_DIFFUSION_BOUNDARY_CURVED
          K == 1   ? ConductivityData::Field::Constant
          : K == 2 ? ConductivityData::Field::Affine
                   : ConductivityData::Field::Quadratic;
#else
          K == 1 ? ConductivityData::Field::Affine : ConductivityData::Field::Quadratic;
#endif
        const Problem problem(mesh, static_cast<typename Problem::Condition>(condition),
          16, field, referenceMean(mesh.getDimension(), field));
        const auto error = problem.solve();
        EXPECT_TRUE(error.isFinite());
        EXPECT_LT(error.getL2(), Real(1e-9));
        EXPECT_LT(error.getH1Seminorm(), Real(1e-9));
      }

      void checkSensitivity(int condition) const
      {
        const auto mesh = makeMesh<ContextType>(this->GetParam(), 3, condition == 2);
        using Problem = BoundaryProblem<2, std::remove_cvref_t<decltype(mesh)>>;
        const auto kind = static_cast<typename Problem::Condition>(condition);
        const Problem base(mesh, kind, 16, ConductivityData::Field::Exponential,
          referenceMean(mesh.getDimension()));
        const Problem higherAssembly(mesh, kind, 18, ConductivityData::Field::Exponential,
          referenceMean(mesh.getDimension()));
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

#ifdef RODIN_DIFFUSION_BOUNDARY_CURVED
#define RODIN_DIFFUSION_BOUNDARY_PATCHES(Name, Prefix, Kind)                             \
  TEST_P(Name, Prefix##ConstantP1Patch)                                                  \
  {                                                                                      \
    checkPatch<1>(Kind);                                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##AffineP2Patch)                                                    \
  {                                                                                      \
    checkPatch<2>(Kind);                                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##QuadraticP4Patch)                                                 \
  {                                                                                      \
    checkPatch<4>(Kind);                                                                 \
  }
#elif defined(RODIN_DIFFUSION_BOUNDARY_VARIABLE)
#define RODIN_DIFFUSION_BOUNDARY_PATCHES(Name, Prefix, Kind)                             \
  TEST_P(Name, Prefix##AffineP1Patch)                                                    \
  {                                                                                      \
    checkPatch<1>(Kind);                                                                 \
  }                                                                                      \
  TEST_P(Name, Prefix##QuadraticP2Patch)                                                 \
  {                                                                                      \
    checkPatch<2>(Kind);                                                                 \
  }
#else
#define RODIN_DIFFUSION_BOUNDARY_PATCHES(Name, Prefix, Kind)
#endif

#define RODIN_POISSON_BOUNDARY_CONDITION(Name, Prefix, Kind)                             \
  RODIN_DIFFUSION_BOUNDARY_PATCHES(Name, Prefix, Kind)                                   \
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

#if defined(RODIN_POISSON_BOUNDARY_MUMPS) || defined(RODIN_DIFFUSION_BOUNDARY_MUMPS)
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
#undef RODIN_DIFFUSION_BOUNDARY_PATCHES
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
