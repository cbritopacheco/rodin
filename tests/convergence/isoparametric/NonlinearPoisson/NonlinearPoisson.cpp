/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "../../CurvedGeometry.h"
#include "../../NonlinearPoisson.h"
#include "../../LiftedConvergence.h"
#include "../../SineMap.h"
#ifdef RODIN_CURVED_NONLINEAR_POISSON_PETSC
#include "../../PETScNonlinearPoisson.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif
#endif
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace Rodin::Tests::Convergence::Isoparametric::NonlinearPoisson
{
  using Data = NonlinearPoissonData;
  constexpr size_t AssemblyOrder = 12, RefinedOrder = 16;
  constexpr Real SolveTolerance = 1e-11, RefinedTolerance = 1e-12;
  // Dimensionless patch/derivative tolerances and relative contamination policy.
  constexpr Real PatchTolerance = 1e-9, TangentTolerance = 1e-6;
  constexpr Real WrongTangentMinimum = 1e-3, SensitivityTolerance = 1e-6;
  constexpr Real ControlL2 = 1e-3, ControlH1 = 1e-2;
  constexpr size_t NormOrder = 14, RefinedNormOrder = 18;
#if defined(RODIN_CURVED_NONLINEAR_POISSON_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif
  /** @brief Mapped semilinear workload with physical data and homogeneous corrections.
   * @par Architecture
   * Exact P2 maps are installed before fields are constructed. Shared native
   * Newton/PETSc SNES workloads supply the solve and derivative oracles.
   * A reference bubble supplies homogeneous perturbations on curved traces.
   */
  template <class ContextType>
  class Workload
  {
    public:
      using Map = typename CurvedGeometry<Mesh<ContextType>>::Map;
      Workload(
        Polytope::Type geometry, size_t n, Map map = Map::Quadratic, bool lifted = false)
        : m_mesh(makeMesh(geometry, n)),
          m_geometry(m_mesh, map)
      {
        if (lifted)
          m_reference.emplace(m_mesh);
        m_geometry.template install<2>();
      }
      template <size_t K>
      ErrorNorms solve(Data::Field field = Data::Field::Sine, bool omitCubic = false,
        size_t order = AssemblyOrder, Real tolerance = SolveTolerance, bool lift = true,
        size_t normOrder = 0, LiftedErrorNorm::Result* lifted = nullptr) const
      {
        const auto observe = [&](const auto& state, const auto& data) {
          if (lifted)
          {
            assert(m_reference);
            *lifted = LiftedErrorNorm::compute(*m_reference, m_mesh, state, data,
              SineMap(), normOrder == 0 ? order : normOrder);
          }
        };
#ifdef RODIN_CURVED_NONLINEAR_POISSON_PETSC
        PETScNonlinearPoissonProblem<K, Mesh<ContextType>> problem(
          m_mesh, order, 1, omitCubic, false, lift, field);
        return problem.solve(tolerance, normOrder, observe);
#else
        NonlinearPoissonProblem problem(m_mesh, order, 1, lift, field);
        return problem.template solve<K>(omitCubic, tolerance, normOrder, observe);
#endif
      }
      template <size_t K>
      Real tangent(bool wrong = false) const
      {
        const RealFunction bubble([this](const Point& p) {
          const auto xi =
            m_geometry.referencePosition(p.getPolytope(), p.getReferenceCoordinates());
          Real value = 1;
          for (size_t j = 0; j < xi.size(); ++j)
            value *= std::sin(Math::Constants::pi() * xi(j));
          return value;
        });
#ifdef RODIN_CURVED_NONLINEAR_POISSON_PETSC
        PETScNonlinearPoissonProblem<K, Mesh<ContextType>> problem(
          m_mesh, AssemblyOrder, 1, false, wrong, true);
        return problem.tangentDefect(bubble);
#else
        NonlinearPoissonProblem problem(m_mesh, AssemblyOrder, 1, true);
        return problem.template tangentDefect<K>(bubble, wrong);
#endif
      }

    private:
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_CURVED_NONLINEAR_POISSON_PETSC) && defined(RODIN_USE_MPI)
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
            .makeMesh(n);
#endif
      }
      Mesh<ContextType> m_mesh;
      Optional<Mesh<ContextType>> m_reference;
      CurvedGeometry<Mesh<ContextType>> m_geometry;
  };
  template <class ContextType>
  class CurvedTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      using Map = typename Workload<ContextType>::Map;

      template <size_t K>
      void approximatedRates() const
      {
        LiftedConvergence history;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
          : this->GetParam() == Polytope::Type::Segment
          ? std::initializer_list<size_t>{5, 9, 17, 33}
          : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload<ContextType> problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrorNorm::Result lifted;
          const auto represented = problem.template solve<K>(Data::Field::Sine, false,
            AssemblyOrder, SolveTolerance, true, NormOrder, &lifted);
          history.append(Real(1) / Real(n - 1), represented, lifted);
        }
        history.expectRates(K);
      }

      template <size_t K>
      void approximatedSensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5, Map::Sine, true);
        std::vector<LiftedConvergence::Components> errors;
        for (size_t i = 0; i < 4; ++i)
        {
          LiftedErrorNorm::Result lifted;
          const auto represented = problem.template solve<K>(Data::Field::Sine, false,
            i == 1 ? RefinedOrder : AssemblyOrder,
            i == 2 ? RefinedTolerance : SolveTolerance, true,
            i == 3 ? RefinedNormOrder : NormOrder, &lifted);
          LiftedConvergence::expectDecomposition(lifted);
          errors.push_back(LiftedConvergence::components(represented, lifted));
        }
        for (size_t i = 1; i < errors.size(); ++i)
        {
          SCOPED_TRACE(::testing::Message() << "control=" << i);
          LiftedConvergence::expectSensitivity(errors[0], errors[i]);
        }
      }

      void approximatedPatchAndControl() const
      {
        Workload<ContextType> problem(this->GetParam(), 5, Map::Sine, true);
        LiftedErrorNorm::Result base, wrong;
        const auto patch = problem.template solve<2>(Data::Field::Affine, false,
          AssemblyOrder, SolveTolerance, true, NormOrder, &base);
        const auto incorrect = problem.template solve<2>(Data::Field::Affine, true,
          AssemblyOrder, SolveTolerance, true, NormOrder, &wrong);
        LiftedConvergence::expectDecomposition(base);
        LiftedConvergence::expectDecomposition(wrong);
        for (const auto& e : {patch, base.field})
        {
          EXPECT_LT(e.getL2(), PatchTolerance);
          EXPECT_LT(e.getH1Seminorm(), PatchTolerance);
        }
        EXPECT_EQ(base.geometry.getL2(), wrong.geometry.getL2());
        EXPECT_EQ(base.geometry.getH1Seminorm(), wrong.geometry.getH1Seminorm());
        for (const auto& e : {incorrect, wrong.field})
        {
          EXPECT_GT(e.getL2(), ControlL2);
          EXPECT_GT(e.getH1Seminorm(), ControlH1);
        }
        constexpr Real ControlRatio = 2;
        EXPECT_GT(wrong.total.getL2(), ControlRatio * base.total.getL2());
        EXPECT_GT(wrong.total.getH1Seminorm(), ControlRatio * base.total.getH1Seminorm());
      }

      template <size_t K>
      void rates() const
      {
        ErrorHistory history;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                                   : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload<ContextType> problem(this->GetParam(), n);
          history.append(Real(1) / Real(n - 1), problem.template solve<K>());
        }
        constexpr Real L2Margin = 0.55, H1Margin = 0.45;
        for (size_t i = 1; i < history.getSize(); ++i)
        {
          const auto& coarse = history.getSample(i - 1).error;
          const auto& fine = history.getSample(i).error;
          ASSERT_TRUE(coarse.isFinite());
          ASSERT_TRUE(fine.isFinite());
          ASSERT_GT(fine.getL2(), 0);
          ASSERT_GT(fine.getH1Seminorm(), 0);
          EXPECT_GT(coarse.getL2(), fine.getL2());
          EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
          const auto rate = history.getAlgebraicRates(i);
          SCOPED_TRACE(::testing::Message() << "interval=" << i << " L2=" << rate.getL2()
                                            << " H1=" << rate.getH1Seminorm());
          EXPECT_GT(rate.getL2(), K + 1 - L2Margin);
          EXPECT_LT(rate.getL2(), K + 1 + L2Margin);
          EXPECT_GT(rate.getH1Seminorm(), K - H1Margin);
          EXPECT_LT(rate.getH1Seminorm(), K + H1Margin);
        }
      }
      template <size_t K>
      void patch() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        const auto e =
          problem.template solve<K>(K == 1 ? Data::Field::Constant : Data::Field::Affine);
        EXPECT_LT(e.getL2(), PatchTolerance);
        EXPECT_LT(e.getH1Seminorm(), PatchTolerance);
      }
      template <size_t K>
      void sensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        const auto base = problem.template solve<K>();
        const auto quad =
          problem.template solve<K>(Data::Field::Sine, false, RefinedOrder);
        const auto solver = problem.template solve<K>(
          Data::Field::Sine, false, AssemblyOrder, RefinedTolerance);
        for (const auto& e : {quad, solver})
          for (const auto& pair : {std::pair{base.getL2(), e.getL2()},
                 std::pair{base.getH1Seminorm(), e.getH1Seminorm()}})
          {
            ASSERT_TRUE(std::isfinite(pair.first));
            ASSERT_TRUE(std::isfinite(pair.second));
            ASSERT_GT(pair.first, 0);
            EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
          }
      }
      void tangent(Map map = Map::Quadratic) const
      {
        Workload<ContextType> problem(this->GetParam(), 3, map);
        EXPECT_LT(problem.template tangent<1>(), TangentTolerance);
        EXPECT_LT(problem.template tangent<2>(), TangentTolerance);
        EXPECT_GT(problem.template tangent<2>(true), WrongTangentMinimum);
      }
      void control() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        const auto e = problem.template solve<2>(Data::Field::Affine, true);
        ASSERT_TRUE(e.isFinite());
        EXPECT_GT(e.getL2(), ControlL2);
        EXPECT_GT(e.getH1Seminorm(), ControlH1);
        const auto boundary = problem.template solve<2>(
          Data::Field::Affine, false, AssemblyOrder, SolveTolerance, false);
        ASSERT_TRUE(boundary.isFinite());
        EXPECT_GT(boundary.getL2(), ControlL2);
        EXPECT_GT(boundary.getH1Seminorm(), ControlH1);
      }
  };
  using LocalTest = CurvedTest<Context::Local>;
  TEST_P(LocalTest, ApproximatedResidualTangentConsistency)
  {
    tangent(Map::Sine);
  }
  TEST_P(LocalTest, ApproximatedP1Rates)
  {
    approximatedRates<1>();
  }
  TEST_P(LocalTest, ApproximatedP2Rates)
  {
    approximatedRates<2>();
  }
  TEST_P(LocalTest, ApproximatedP1Sensitivity)
  {
    approximatedSensitivity<1>();
  }
  TEST_P(LocalTest, ApproximatedP2Sensitivity)
  {
    approximatedSensitivity<2>();
  }
  TEST_P(LocalTest, ApproximatedAffinePatchRejectsOmittedCubic)
  {
    approximatedPatchAndControl();
  }
  TEST_P(LocalTest, P1Rates)
  {
    rates<1>();
  }
  TEST_P(LocalTest, P2Rates)
  {
    rates<2>();
  }
  TEST_P(LocalTest, ConstantP1Patch)
  {
    patch<1>();
  }
  TEST_P(LocalTest, AffineP2Patch)
  {
    patch<2>();
  }
  TEST_P(LocalTest, P1Sensitivity)
  {
    sensitivity<1>();
  }
  TEST_P(LocalTest, P2Sensitivity)
  {
    sensitivity<2>();
  }
  TEST_P(LocalTest, TangentConsistency)
  {
    tangent();
  }
  TEST_P(LocalTest, RejectsWrongPhysicsAndBoundary)
  {
    control();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#if defined(RODIN_CURVED_NONLINEAR_POISSON_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = CurvedTest<Context::MPI>;
  TEST_P(MPITest, ApproximatedResidualTangentConsistency)
  {
    tangent(Map::Sine);
  }
  TEST_P(MPITest, ApproximatedP1Rates)
  {
    approximatedRates<1>();
  }
  TEST_P(MPITest, ApproximatedP2Rates)
  {
    approximatedRates<2>();
  }
  TEST_P(MPITest, ApproximatedP1Sensitivity)
  {
    approximatedSensitivity<1>();
  }
  TEST_P(MPITest, ApproximatedP2Sensitivity)
  {
    approximatedSensitivity<2>();
  }
  TEST_P(MPITest, ApproximatedAffinePatchRejectsOmittedCubic)
  {
    approximatedPatchAndControl();
  }
  TEST_P(MPITest, P1Rates)
  {
    rates<1>();
  }
  TEST_P(MPITest, P2Rates)
  {
    rates<2>();
  }
  TEST_P(MPITest, ConstantP1Patch)
  {
    patch<1>();
  }
  TEST_P(MPITest, AffineP2Patch)
  {
    patch<2>();
  }
  TEST_P(MPITest, P1Sensitivity)
  {
    sensitivity<1>();
  }
  TEST_P(MPITest, P2Sensitivity)
  {
    sensitivity<2>();
  }
  TEST_P(MPITest, TangentConsistency)
  {
    tangent();
  }
  TEST_P(MPITest, RejectsWrongPhysicsAndBoundary)
  {
    control();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, MPITest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#endif
}
#ifdef RODIN_CURVED_NONLINEAR_POISSON_PETSC
int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::NonlinearPoisson::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::NonlinearPoisson::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
