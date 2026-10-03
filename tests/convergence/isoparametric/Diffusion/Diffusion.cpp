/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Poisson and variable conductivity on exact quadratic geometry. */
#include <gtest/gtest.h>
#include "../../Conductivity.h"
#include "../../CurvedGeometry.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"
#ifdef RODIN_CURVED_DIFFUSION_PETSC
#include "Rodin/PETSc.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::Diffusion
{
  using Data = ConductivityData;
  // Non-polynomial mapped integrands: budgets are checked by sensitivity tests.
  constexpr size_t AssemblyOrder = 11, NormOrder = 13, RefinedOrder = 16;
  constexpr size_t RefinedNormOrder = RefinedOrder + 2;
  constexpr Real SolverTolerance = 1e-13, RefinedTolerance = 1e-14;
  constexpr Real ResidualTolerance = 1e-11, PatchTolerance = 1e-9;
  constexpr Real SensitivityTolerance = 1e-6;
  constexpr size_t MaxIterations = 50000;
  [[maybe_unused]] constexpr Real DivergenceTolerance = 1e5;
  // Finite-resolution acceptance policies, not asymptotic estimates.
  constexpr Real L2Margin = 0.55, H1Margin = 0.45;
  // Absolute dimensionless wrong-operator separation on the stated domains.
  constexpr Real PoissonControlL2 = 0.02, PoissonControlH1 = 0.05;
  constexpr Real ConductivityControlL2 = 1e-3, ConductivityControlH1 = 1e-2;
  constexpr Real ControlRatio = 2;
  constexpr Real MapTolerance = 1e-11, VolumeTolerance = 1e-12;
#if defined(RODIN_CURVED_DIFFUSION_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  /** @brief Shared curved mesh with independent Poisson/conductivity systems.
   * @par Architecture
   * Exact P2 maps are installed before creating spaces. Every solve owns fresh
   * fields, forms and linear algebra. Continuous data select the coefficient
   * before deriving the load. Wrong operators retain that correct load and
   * boundary trace. Physical error integration is independent of assembly;
   * owned-cell squared norms are reduced globally for MPI.
   */
  template <class ContextType>
  class Workload
  {
    public:
      Workload(Polytope::Type geometry, size_t n)
        : m_dimension(UniformGrid(geometry).getDimension()),
          m_mesh(makeMesh(geometry, n)),
          m_geometry(m_mesh)
      {
        m_geometry.template install<2>();
      }
      const auto& getMesh() const
      {
        return m_mesh;
      }
      const auto& getGeometry() const
      {
        return m_geometry;
      }
      size_t getDimension() const
      {
        return m_dimension;
      }

      template <size_t K>
      ErrorNorms solve(bool poisson, Data::Field field, bool wrong = false,
        size_t order = AssemblyOrder, Real tolerance = SolverTolerance,
        size_t normOrder = NormOrder) const
      {
        const Data data(m_dimension, field);
        const auto exact = data.getSolution();
        const auto source = data.getSource(poisson);
        const auto gamma = data.getCoefficient(poisson || wrong);
        auto space = [&] {
#ifndef RODIN_CURVED_DIFFUSION_PETSC
          if constexpr (K == 1)
            return P1(m_mesh);
          else
#endif
            return H1<K, Real, Mesh<ContextType>>(
              std::integral_constant<size_t, K>{}, m_mesh);
        }();
#ifdef RODIN_CURVED_DIFFUSION_PETSC
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
#else
        TrialFunction u(space);
        TestFunction v(space);
#endif
        // Doubling Poisson stiffness or omitting variable conductivity leaves
        // manufactured data unchanged. Both are deliberately incorrect PDEs.
        const Real scale = poisson && wrong ? Real(2) : Real(1);
        auto stiffness = Integral(scale * gamma * Grad(u), Grad(v));
        auto load = Integral(source, v);
        stiffness.setOrder(order);
        load.setOrder(order);
        Problem problem(u, v);
        problem = stiffness - load + DirichletBC(u, exact);
#ifdef RODIN_CURVED_DIFFUSION_PETSC
        PETSc::Solver::CG solver(problem);
        solver
          .setTolerances(tolerance, RefinedTolerance, DivergenceTolerance, MaxIterations)
          .solve();
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
        return ErrorNorm::compute(
          m_mesh, u.getSolution(), exact, data.getGradient(), normOrder);
      }

    private:
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_CURVED_DIFFUSION_PETSC) && defined(RODIN_USE_MPI)
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
            .makeMesh(n);
#endif
      }
      // The continuous dimension is fixed by the grid family. Empty shards
      // may report a local mesh dimension of zero; that does not change data.
      size_t m_dimension;
      Mesh<ContextType> m_mesh;
      CurvedGeometry<Mesh<ContextType>> m_geometry;
  };

  template <class ContextType>
  class DiffusionTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      template <size_t K>
      void rates() const
      {
        std::array<ErrorHistory, 2> histories;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                                   : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload<ContextType> problem(this->GetParam(), n);
          for (bool poisson : {true, false})
          {
            SCOPED_TRACE(poisson ? "Poisson" : "Conductivity");
            histories[poisson ? 0 : 1].append(Real(1) / Real(n - 1),
              problem.template solve<K>(poisson, Data::Field::Smooth));
          }
        }
        for (size_t physics = 0; physics < 2; ++physics)
          for (size_t i = 1; i < histories[physics].getSize(); ++i)
          {
            const auto& coarse = histories[physics].getSample(i - 1).error;
            const auto& fine = histories[physics].getSample(i).error;
            const auto rate = histories[physics].getAlgebraicRates(i);
            SCOPED_TRACE(::testing::Message()
              << (physics == 0 ? "Poisson" : "Conductivity") << " interval=" << i
              << " L2 " << coarse.getL2() << " -> " << fine.getL2() << " H1 "
              << coarse.getH1Seminorm() << " -> " << fine.getH1Seminorm()
              << " rates=" << rate.getL2() << "," << rate.getH1Seminorm());
            ASSERT_TRUE(coarse.isFinite());
            ASSERT_TRUE(fine.isFinite());
            ASSERT_GT(fine.getL2(), 0);
            ASSERT_GT(fine.getH1Seminorm(), 0);
            EXPECT_GT(coarse.getL2(), fine.getL2());
            EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
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
        for (bool poisson : {true, false})
        {
          SCOPED_TRACE(poisson ? "Poisson" : "Conductivity");
          const auto e = problem.template solve<K>(
            poisson, K == 1 ? Data::Field::Constant : Data::Field::Affine);
          EXPECT_LT(e.getL2(), PatchTolerance);
          EXPECT_LT(e.getH1Seminorm(), PatchTolerance);
        }
      }
      void controls() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        for (bool poisson : {true, false})
        {
          SCOPED_TRACE(poisson ? "Poisson" : "Conductivity");
          const auto field = poisson ? Data::Field::Smooth : Data::Field::Affine;
          const auto base = problem.template solve<2>(poisson, field);
          const auto wrong = problem.template solve<2>(poisson, field, true);
          ASSERT_TRUE(base.isFinite());
          ASSERT_TRUE(wrong.isFinite());
          EXPECT_GT(wrong.getL2(), poisson ? PoissonControlL2 : ConductivityControlL2);
          EXPECT_GT(
            wrong.getH1Seminorm(), poisson ? PoissonControlH1 : ConductivityControlH1);
          EXPECT_GT(wrong.getL2(), ControlRatio * base.getL2());
          EXPECT_GT(wrong.getH1Seminorm(), ControlRatio * base.getH1Seminorm());
          if (!poisson)
          {
            EXPECT_LT(base.getL2(), PatchTolerance);
            EXPECT_LT(base.getH1Seminorm(), PatchTolerance);
          }
        }
      }
      template <size_t K>
      void sensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        for (bool poisson : {true, false})
        {
          SCOPED_TRACE(poisson ? "Poisson" : "Conductivity");
          const auto base = problem.template solve<K>(poisson, Data::Field::Smooth);
          const auto quad =
            problem.template solve<K>(poisson, Data::Field::Smooth, false, RefinedOrder);
          const auto solver = problem.template solve<K>(
            poisson, Data::Field::Smooth, false, AssemblyOrder, RefinedTolerance);
          const auto norm = problem.template solve<K>(poisson, Data::Field::Smooth, false,
            AssemblyOrder, SolverTolerance, RefinedNormOrder);
          for (const auto& e : {quad, solver, norm})
            for (const auto pair : {std::pair{base.getL2(), e.getL2()},
                   std::pair{base.getH1Seminorm(), e.getH1Seminorm()}})
            {
              ASSERT_TRUE(std::isfinite(pair.first));
              ASSERT_TRUE(std::isfinite(pair.second));
              ASSERT_GT(pair.first, 0);
              EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
            }
        }
      }
      void geometry() const
      {
        Workload<ContextType> problem(this->GetParam(), 2);
        const auto& mesh = problem.getMesh();
        const auto& mapping = problem.getGeometry();
        const size_t dim = problem.getDimension();
        constexpr size_t GeometryCheckOrder = 4;
        for (size_t d = 1; d <= dim; ++d)
          for (auto polytope = mesh.getPolytope(d); polytope; ++polytope)
          {
            const RealH1Element<2> element(polytope->getGeometry());
            EXPECT_EQ(polytope->getTransformation().getOrder(), element.getOrder());
            const auto& qf = QF::PolytopeQuadratureFormula::get(
              GeometryCheckOrder, polytope->getGeometry());
            const auto& quadrature = polytope->getQuadrature(qf);
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& point = quadrature.getPoint(qp);
              const auto expected = mapping.mapToPhysical(
                mapping.referencePosition(*polytope, point.getReferenceCoordinates()));
              EXPECT_LT((point.getPhysicalCoordinates() - expected).norm(), MapTolerance);
              EXPECT_TRUE(std::isfinite(point.getDistortion()));
              EXPECT_GT(point.getDistortion(), 0);
            }
          }
        EXPECT_NEAR(mesh.getMeasure(dim), dim == 1 ? 1.1 : 1, VolumeTolerance);
      }
  };
  using LocalTest = DiffusionTest<Context::Local>;
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
  TEST_P(LocalTest, AffinePhysicalP2Patch)
  {
    patch<2>();
  }
  TEST_P(LocalTest, RejectsWrongOperators)
  {
    controls();
  }
  TEST_P(LocalTest, P1Sensitivity)
  {
    sensitivity<1>();
  }
  TEST_P(LocalTest, P2Sensitivity)
  {
    sensitivity<2>();
  }
  TEST_P(LocalTest, ExactGeometryAndVolume)
  {
    geometry();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#if defined(RODIN_CURVED_DIFFUSION_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = DiffusionTest<Context::MPI>;
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
  TEST_P(MPITest, AffinePhysicalP2Patch)
  {
    patch<2>();
  }
  TEST_P(MPITest, RejectsWrongOperators)
  {
    controls();
  }
  TEST_P(MPITest, P1Sensitivity)
  {
    sensitivity<1>();
  }
  TEST_P(MPITest, P2Sensitivity)
  {
    sensitivity<2>();
  }
  TEST_P(MPITest, ExactGeometryAndVolume)
  {
    geometry();
  }
  TEST_P(MPITest, GlobalNormCountsOwnedCellsOnce)
  {
    Workload<Context::MPI> problem(GetParam(), 3);
    const auto& mesh = problem.getMesh();
    const size_t dim = problem.getDimension();
    H1<2, Real, Mesh<Context::MPI>> space(std::integral_constant<size_t, 2>{}, mesh);
    PETSc::Variational::GridFunction u(space);
    u = RealFunction(0);
    const auto derivative = VectorFunction(
      dim, [dim](const Point&) { return Math::SpatialVector<Real>::Zero(dim); });
    const auto e = ErrorNorm::compute(mesh, u, RealFunction(1), derivative, NormOrder);
    EXPECT_NEAR(e.getL2(), std::sqrt(dim == 1 ? Real(1.1) : Real(1)), VolumeTolerance);
    EXPECT_EQ(e.getH1Seminorm(), 0);
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
#ifdef RODIN_CURVED_DIFFUSION_PETSC
int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::Diffusion::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::Diffusion::world = &comm;
#endif
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
