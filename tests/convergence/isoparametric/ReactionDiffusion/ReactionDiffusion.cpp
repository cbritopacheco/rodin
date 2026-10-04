/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file
 * @brief Coupled scalar fields on exact quadratic and approximated sine geometry.
 */

#include <gtest/gtest.h>
#include "Rodin/Assembly.h"
#include "Rodin/Solver.h"
#include "../../ReactionDiffusion.h"
#include "../../CurvedGeometry.h"
#include "../../LiftedConvergence.h"
#include "../../SineMap.h"

#ifdef RODIN_CURVED_REACTION_DIFFUSION_PETSC
#include "Rodin/PETSc.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::ReactionDiffusion
{
  using Data = ReactionDiffusionData;
  // Dimensionless symmetric positive-definite reaction matrix.
  constexpr Real Coupling = 0.2;
  // Non-polynomial mapped integrands require sensitivity checks, not a
  // polynomial-exactness claim. Norms use two additional quadrature degrees.
  constexpr size_t AssemblyOrder = 11;
  constexpr size_t RefinedOrder = 16;
  constexpr Real SolverTolerance = 1e-13;
  constexpr Real RefinedTolerance = 1e-14;
  constexpr Real ResidualTolerance = 1e-11;
  constexpr size_t MaxIterations = 50000;
  // PETSc's divergence threshold is deliberately remote from the error budget.
  [[maybe_unused]] constexpr Real DivergenceTolerance = 1e5;
  // Roundoff allowance for dimensionless representable fields and derivatives.
  constexpr Real PatchTolerance = 1e-9;
  // Relative contamination policy, well below the adjacent-rate margins.
  constexpr Real SensitivityTolerance = 1e-6;
  // Absolute dimensionless errors separating omitted coupling from roundoff.
  constexpr Real ControlL2 = 1e-3;
  constexpr Real ControlH1 = 1e-2;
  constexpr size_t NormOrder = 13, RefinedNormOrder = 18;
  constexpr Real ControlRatio = 2;

#if defined(RODIN_CURVED_REACTION_DIFFUSION_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  /** @brief Coupled mapped workload with shared physical data and error oracles.
   * @par Architecture
   * The mesh owns P2 maps on cells and traces: exact for the quadratic map,
   * interpolated for the sine map. Each solve builds fresh
   * scalar spaces and a two-field system. Backend policy is explicit; norms
   * integrate physical errors on owned cells before the MPI reduction.
   */
  template <class ContextType>
  class Workload
  {
    public:
      using Map = typename CurvedGeometry<Mesh<ContextType>>::Map;
      using LiftedErrors = std::array<LiftedErrorNorm::Result, 2>;
      Workload(
        Polytope::Type geometry, size_t n, Map map = Map::Quadratic, bool lifted = false)
        : m_mesh(makeMesh(geometry, n)),
          m_geometry(m_mesh, map)
      {
        if (lifted)
          m_reference.emplace(m_mesh);
        m_geometry.template install<2>();
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }

      template <size_t K>
      std::array<ErrorNorms, 2> solve(Data::Field field, bool omitCoupling = false,
        size_t order = AssemblyOrder, Real tolerance = SolverTolerance,
        size_t normOrder = 0, LiftedErrors* lifted = nullptr) const
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Data data(dim, field);
        const auto exactU = data.getSolution(0), exactW = data.getSolution(1);
        const auto sourceU = data.getSource(0), sourceW = data.getSource(1);
        auto space = [&] {
#ifndef RODIN_CURVED_REACTION_DIFFUSION_PETSC
          if constexpr (K == 1)
            return P1(m_mesh);
          else
#endif
            return H1<K, Real, Mesh<ContextType>>(
              std::integral_constant<size_t, K>{}, m_mesh);
        }();
#ifdef RODIN_CURVED_REACTION_DIFFUSION_PETSC
        PETSc::Variational::TrialFunction u(space), w(space);
        PETSc::Variational::TestFunction v(space), z(space);
#else
        TrialFunction u(space), w(space);
        TestFunction v(space), z(space);
#endif
        const Real alpha = omitCoupling ? Real(0) : Coupling;
        auto diffusionU = Integral(Grad(u), Grad(v));
        auto diffusionW = Integral(Real(2) * Grad(w), Grad(z));
        auto reactionU = Integral(u, v), reactionW = Integral(w, z);
        auto couplingU = Integral(alpha * w, v), couplingW = Integral(alpha * u, z);
        auto bodyU = Integral(sourceU, v), bodyW = Integral(sourceW, z);
        diffusionU.setOrder(order);
        diffusionW.setOrder(order);
        reactionU.setOrder(order);
        reactionW.setOrder(order);
        couplingU.setOrder(order);
        couplingW.setOrder(order);
        bodyU.setOrder(order);
        bodyW.setOrder(order);
        Problem problem(u, w, v, z);
        problem = diffusionU + reactionU + couplingU - bodyU + diffusionW + reactionW +
          couplingW - bodyW + DirichletBC(u, exactU) + DirichletBC(w, exactW);
#ifdef RODIN_CURVED_REACTION_DIFFUSION_PETSC
        PETSc::Solver::CG solver(problem);
        solver.setTolerances(
          tolerance, RefinedTolerance, DivergenceTolerance, MaxIterations);
        solver.solve();
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        EXPECT_GT(reason, 0);
        const auto& system = problem.getLinearSystem();
        Vec residual = nullptr;
        EXPECT_EQ(VecDuplicate(system.getVector(), &residual), PETSC_SUCCESS);
        EXPECT_EQ(
          MatMult(system.getOperator(), system.getSolution(), residual), PETSC_SUCCESS);
        EXPECT_EQ(VecAXPY(residual, -1, system.getVector()), PETSC_SUCCESS);
        PetscReal norm = 0, rhsNorm = 0;
        EXPECT_EQ(VecNorm(residual, NORM_2, &norm), PETSC_SUCCESS);
        EXPECT_EQ(VecNorm(system.getVector(), NORM_2, &rhsNorm), PETSC_SUCCESS);
        const Real relative = norm / std::max(Real(1), rhsNorm);
        EXPECT_EQ(VecDestroy(&residual), PETSC_SUCCESS);
#else
        Solver::CG solver(problem);
        solver.setTolerance(tolerance).setMaxIterations(MaxIterations).solve();
        EXPECT_TRUE(solver.success());
        const auto& system = problem.getLinearSystem();
        const auto residual =
          system.getOperator() * system.getSolution() - system.getVector();
        const Real relative =
          residual.norm() / std::max(Real(1), system.getVector().norm());
#endif
        EXPECT_TRUE(std::isfinite(relative));
        EXPECT_LT(relative, ResidualTolerance);
        const size_t integrationOrder = normOrder == 0 ? order + 2 : normOrder;
        if (lifted)
        {
          assert(m_reference);
          // Adapter selects a scalar component without duplicating analytic data.
          struct Component
          {
              const Data& data;
              size_t index;
              Real getSolution(const Math::SpatialPoint& x) const
              {
                return data.getSolution(x, index);
              }
              auto getGradient(const Math::SpatialPoint& x) const
              {
                return data.getGradient(x, index);
              }
          };
          (*lifted)[0] = LiftedErrorNorm::compute(*m_reference, m_mesh, u.getSolution(),
            Component{data, 0}, SineMap(), integrationOrder);
          (*lifted)[1] = LiftedErrorNorm::compute(*m_reference, m_mesh, w.getSolution(),
            Component{data, 1}, SineMap(), integrationOrder);
        }
        return {ErrorNorm::compute(
                  m_mesh, u.getSolution(), exactU, data.getGradient(0), integrationOrder),
          ErrorNorm::compute(
            m_mesh, w.getSolution(), exactW, data.getGradient(1), integrationOrder)};
      }

    private:
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_CURVED_REACTION_DIFFUSION_PETSC) && defined(RODIN_USE_MPI)
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
  class CurvedReactionDiffusionTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      using Map = typename Workload<ContextType>::Map;
      using LiftedErrors = typename Workload<ContextType>::LiftedErrors;

      template <size_t K>
      void checkApproximatedRates() const
      {
        std::array<LiftedConvergence, 2> histories;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
          : this->GetParam() == Polytope::Type::Segment
          ? std::initializer_list<size_t>{5, 9, 17, 33}
          : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload<ContextType> problem(this->GetParam(), n, Map::Sine, true);
          LiftedErrors lifted;
          const auto represented = problem.template solve<K>(Data::Field::Smooth, false,
            AssemblyOrder, SolverTolerance, NormOrder, &lifted);
          for (size_t component = 0; component < histories.size(); ++component)
          {
            SCOPED_TRACE(::testing::Message() << "field=" << component);
            histories[component].append(
              Real(1) / Real(n - 1), represented[component], lifted[component]);
          }
        }
        for (size_t component = 0; component < histories.size(); ++component)
        {
          SCOPED_TRACE(::testing::Message() << "field=" << component);
          histories[component].expectRates(K);
        }
      }

      template <size_t K>
      void checkApproximatedSensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5, Map::Sine, true);
        std::vector<std::array<LiftedConvergence::Components, 2>> errors;
        for (size_t i = 0; i < 4; ++i)
        {
          LiftedErrors lifted;
          const auto represented = problem.template solve<K>(Data::Field::Smooth, false,
            i == 1 ? RefinedOrder : AssemblyOrder,
            i == 2 ? RefinedTolerance : SolverTolerance,
            i == 3 ? RefinedNormOrder : NormOrder, &lifted);
          for (size_t component = 0; component < 2; ++component)
          {
            LiftedConvergence::expectDecomposition(lifted[component]);
          }
          errors.push_back({LiftedConvergence::components(represented[0], lifted[0]),
            LiftedConvergence::components(represented[1], lifted[1])});
        }
        for (size_t i = 1; i < errors.size(); ++i)
          for (size_t component = 0; component < 2; ++component)
          {
            SCOPED_TRACE(
              ::testing::Message() << "control=" << i << " field=" << component);
            LiftedConvergence::expectSensitivity(
              errors[0][component], errors[i][component]);
          }
      }

      void checkApproximatedPatchAndControl() const
      {
        Workload<ContextType> problem(this->GetParam(), 5, Map::Sine, true);
        LiftedErrors base, wrong;
        const auto patch = problem.template solve<2>(
          Data::Field::Affine, false, AssemblyOrder, SolverTolerance, NormOrder, &base);
        const auto incorrect = problem.template solve<2>(
          Data::Field::Affine, true, AssemblyOrder, SolverTolerance, NormOrder, &wrong);
        for (size_t component = 0; component < 2; ++component)
        {
          SCOPED_TRACE(::testing::Message() << "field=" << component);
          LiftedConvergence::expectDecomposition(base[component]);
          LiftedConvergence::expectDecomposition(wrong[component]);
          for (const auto& e : {patch[component], base[component].field})
          {
            EXPECT_LT(e.getL2(), PatchTolerance);
            EXPECT_LT(e.getH1Seminorm(), PatchTolerance);
          }
          EXPECT_EQ(base[component].geometry.getL2(), wrong[component].geometry.getL2());
          EXPECT_EQ(base[component].geometry.getH1Seminorm(),
            wrong[component].geometry.getH1Seminorm());
          for (const auto& e : {incorrect[component], wrong[component].field})
          {
            EXPECT_GT(e.getL2(), ControlL2);
            EXPECT_GT(e.getH1Seminorm(), ControlH1);
          }
          EXPECT_GT(
            wrong[component].total.getL2(), ControlRatio * base[component].total.getL2());
          EXPECT_GT(wrong[component].total.getH1Seminorm(),
            ControlRatio * base[component].total.getH1Seminorm());
        }
      }

      template <size_t K>
      void checkRates() const
      {
        std::array<ErrorHistory, 2> histories;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                                   : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload<ContextType> problem(this->GetParam(), n);
          const auto errors = problem.template solve<K>(Data::Field::Smooth);
          for (size_t component = 0; component < 2; ++component)
            histories[component].append(Real(1) / Real(n - 1), errors[component]);
        }
        // Finite-resolution policy windows around the theoretical orders.
        constexpr Real L2Margin = 0.55, H1Margin = 0.45;
        for (size_t component = 0; component < 2; ++component)
          for (size_t i = 1; i < histories[component].getSize(); ++i)
          {
            const auto& coarse = histories[component].getSample(i - 1).error;
            const auto& fine = histories[component].getSample(i).error;
            ASSERT_TRUE(coarse.isFinite());
            ASSERT_TRUE(fine.isFinite());
            ASSERT_GT(fine.getL2(), 0);
            ASSERT_GT(fine.getH1Seminorm(), 0);
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
            const auto rate = histories[component].getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << "component=" << component << " interval=" << i << " L2=" << rate.getL2()
              << " H1=" << rate.getH1Seminorm());
            EXPECT_GT(rate.getL2(), K + 1 - L2Margin);
            EXPECT_LT(rate.getL2(), K + 1 + L2Margin);
            EXPECT_GT(rate.getH1Seminorm(), K - H1Margin);
            EXPECT_LT(rate.getH1Seminorm(), K + H1Margin);
          }
      }

      template <size_t K>
      void checkPatch() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        for (const auto& error :
          problem.template solve<K>(K == 1 ? Data::Field::Constant : Data::Field::Affine))
        {
          EXPECT_LT(error.getL2(), PatchTolerance);
          EXPECT_LT(error.getH1Seminorm(), PatchTolerance);
        }
      }

      void checkControl() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        const auto baseline = problem.template solve<2>(Data::Field::Affine);
        const auto wrong = problem.template solve<2>(Data::Field::Affine, true);
        for (size_t component = 0; component < 2; ++component)
        {
          ASSERT_TRUE(baseline[component].isFinite());
          ASSERT_TRUE(wrong[component].isFinite());
          EXPECT_LT(baseline[component].getL2(), PatchTolerance);
          EXPECT_LT(baseline[component].getH1Seminorm(), PatchTolerance);
          EXPECT_GT(wrong[component].getL2(), ControlL2);
          EXPECT_GT(wrong[component].getH1Seminorm(), ControlH1);
        }
      }

      template <size_t K>
      void checkSensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        const auto baseline = problem.template solve<K>(Data::Field::Smooth);
        const auto quadrature =
          problem.template solve<K>(Data::Field::Smooth, false, RefinedOrder);
        const auto solver = problem.template solve<K>(
          Data::Field::Smooth, false, AssemblyOrder, RefinedTolerance);
        for (const auto& refined : {quadrature, solver})
          for (size_t component = 0; component < 2; ++component)
            for (const auto pair :
              {std::pair{baseline[component].getL2(), refined[component].getL2()},
                std::pair{baseline[component].getH1Seminorm(),
                  refined[component].getH1Seminorm()}})
            {
              ASSERT_TRUE(std::isfinite(pair.first));
              ASSERT_TRUE(std::isfinite(pair.second));
              ASSERT_GT(pair.first, 0);
              EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
            }
      }
  };

  using LocalTest = CurvedReactionDiffusionTest<Context::Local>;
  TEST_P(LocalTest, ApproximatedP1Rates)
  {
    checkApproximatedRates<1>();
  }
  TEST_P(LocalTest, ApproximatedP2Rates)
  {
    checkApproximatedRates<2>();
  }
  TEST_P(LocalTest, ApproximatedP1Sensitivity)
  {
    checkApproximatedSensitivity<1>();
  }
  TEST_P(LocalTest, ApproximatedP2Sensitivity)
  {
    checkApproximatedSensitivity<2>();
  }
  TEST_P(LocalTest, ApproximatedAffinePatchRejectsOmittedCoupling)
  {
    checkApproximatedPatchAndControl();
  }
  TEST_P(LocalTest, P1OptimalRates)
  {
    checkRates<1>();
  }
  TEST_P(LocalTest, P2OptimalRates)
  {
    checkRates<2>();
  }
  TEST_P(LocalTest, ConstantP1Patch)
  {
    checkPatch<1>();
  }
  TEST_P(LocalTest, AffinePhysicalP2Patch)
  {
    checkPatch<2>();
  }
  TEST_P(LocalTest, RejectsOmittedCoupling)
  {
    checkControl();
  }
  TEST_P(LocalTest, P1Sensitivity)
  {
    checkSensitivity<1>();
  }
  TEST_P(LocalTest, P2Sensitivity)
  {
    checkSensitivity<2>();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });

#if defined(RODIN_CURVED_REACTION_DIFFUSION_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = CurvedReactionDiffusionTest<Context::MPI>;
  TEST_P(MPITest, ApproximatedP1Rates)
  {
    checkApproximatedRates<1>();
  }
  TEST_P(MPITest, ApproximatedP2Rates)
  {
    checkApproximatedRates<2>();
  }
  TEST_P(MPITest, ApproximatedP1Sensitivity)
  {
    checkApproximatedSensitivity<1>();
  }
  TEST_P(MPITest, ApproximatedP2Sensitivity)
  {
    checkApproximatedSensitivity<2>();
  }
  TEST_P(MPITest, ApproximatedAffinePatchRejectsOmittedCoupling)
  {
    checkApproximatedPatchAndControl();
  }
  TEST_P(MPITest, P1OptimalRates)
  {
    checkRates<1>();
  }
  TEST_P(MPITest, P2OptimalRates)
  {
    checkRates<2>();
  }
  TEST_P(MPITest, ConstantP1Patch)
  {
    checkPatch<1>();
  }
  TEST_P(MPITest, AffinePhysicalP2Patch)
  {
    checkPatch<2>();
  }
  TEST_P(MPITest, RejectsOmittedCoupling)
  {
    checkControl();
  }
  TEST_P(MPITest, P1Sensitivity)
  {
    checkSensitivity<1>();
  }
  TEST_P(MPITest, P2Sensitivity)
  {
    checkSensitivity<2>();
  }
  TEST_P(MPITest, GlobalNormCountsOwnedCurvedCellsOnce)
  {
    Workload<Context::MPI> problem(GetParam(), 5);
    const auto& mesh = problem.getMesh();
    const size_t dim = mesh.getDimension();
    H1<2, Real, Mesh<Context::MPI>> space(std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::GridFunction u(space);
    u = RealFunction([](const Point&) { return Real(0); });
    const auto zeroGradient = VectorFunction(
      dim, [dim](const Point&) { return Math::SpatialVector<Real>::Zero(dim); });
    const Real volume = dim == 1 ? Real(1.1) : Real(1);
    constexpr Real VolumeTolerance = 1e-12; // Analytic volume roundoff allowance.
    for (size_t component = 0; component < 2; ++component)
    {
      const Real c = Real(component + 1);
      const auto error = ErrorNorm::compute(mesh, u,
        RealFunction([c](const Point&) { return c; }), zeroGradient, AssemblyOrder);
      EXPECT_NEAR(error.getL2(), c * std::sqrt(volume), VolumeTolerance);
      EXPECT_EQ(error.getH1Seminorm(), 0);
    }
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

#ifdef RODIN_CURVED_REACTION_DIFFUSION_PETSC
int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::ReactionDiffusion::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::ReactionDiffusion::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
