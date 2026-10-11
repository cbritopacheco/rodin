/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "../../CurvedGeometry.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"
#ifdef RODIN_CURVED_P0_PETSC
#include "Rodin/PETSc.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/P0/P0.h"
#include "Rodin/MPI/Variational/P0g/P0g.h"
#endif
#endif
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace Rodin::Tests::Convergence::Isoparametric::P0Projection
{
  // Physical affine data have quadratic pullbacks. Higher order independently
  // checks moments and errors; no point values are read from coefficients.
  constexpr size_t AssemblyOrder = 8, OracleOrder = 12, RefinedOrder = 14;
  constexpr Real SolverTolerance = 1e-13, RefinedTolerance = 1e-14;
  constexpr Real ResidualTolerance = 1e-11, MomentTolerance = 1e-10;
  constexpr Real PatchTolerance = 1e-10, SensitivityTolerance = 1e-6;
  // Absolute dimensionless controls on unit-scale domains.
  constexpr Real CentroidDefectMinimum = 1e-7, WrongMassErrorMinimum = 0.1;
  constexpr size_t MaxIterations = 20000;
#if defined(RODIN_CURVED_P0_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif
  struct Errors
  {
      Real l2;
      Real moment;
      Real globalMean;
  };
  /** @brief Curved scalar/vector projection with independently integrated moments.
   * @par Architecture
   * A fresh mapped mesh owns P2 geometry. Each solve assembles a mass problem.
   * Independent physical integration checks cell orthogonality for P0, and the
   * analytic global mean for P0g. MPI moments and norms count owned cells only.
   */
  template <class ContextType>
  class Workload
  {
    public:
      Workload(Polytope::Type geometry, size_t n)
        : m_mesh(makeMesh(geometry, n)),
          m_geometry(m_mesh)
      {
        m_geometry.template install<2>();
      }
      template <class Scalar, bool Vector, bool Global = false>
      Errors project(bool constant = false, size_t order = AssemblyOrder,
        Real tolerance = SolverTolerance, bool centroid = false, Real massScale = 1) const
      {
        const size_t dim = m_mesh.getDimension();
        auto space = [&] {
          using Range = std::conditional_t<Vector, Math::SpatialVector<Scalar>, Scalar>;
          if constexpr (Global)
          {
            if constexpr (Vector)
              return P0g<Range, Mesh<ContextType>>(m_mesh, 2);
            else
              return P0g<Range, Mesh<ContextType>>(m_mesh);
          }
          else
          {
            if constexpr (Vector)
              return P0<Range, Mesh<ContextType>>(m_mesh, 2);
            else
              return P0<Range, Mesh<ContextType>>(m_mesh);
          }
        }();
        const auto exact = field<Scalar, Vector>(dim, constant);
#ifdef RODIN_CURVED_P0_PETSC
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
#else
        TrialFunction u(space);
        TestFunction v(space);
#endif
        auto mass = Integral(massScale * u, v);
        auto load = Integral(exact, v);
        mass.setOrder(order);
        load.setOrder(order);
        Problem problem(u, v);
        problem = mass - load;
#ifdef RODIN_CURVED_P0_PETSC
        PETSc::Solver::CG solver(problem);
        constexpr Real DivergenceTolerance = 1e5;
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
        auto& uh = u.getSolution();
        if (centroid)
          uh = exact;
        Real momentSquared = 0;
        if constexpr (!Global)
          for (auto cell = m_mesh.getCell(); cell; ++cell)
          {
            if constexpr (requires { m_mesh.getShard(); })
            {
              if (!m_mesh.getShard().isOwned(dim, cell->getIndex()))
                continue;
            }
            const auto& qf =
              QF::PolytopeQuadratureFormula::get(OracleOrder, cell->getGeometry());
            const Geometry::PolytopeQuadrature quadrature(*cell, qf);
            auto moment = exact(quadrature.getPoint(0));
            moment *= 0;
            Real volume = 0;
            for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
            {
              const auto& p = quadrature.getPoint(qp);
              const IntegrationPoint ip(p, &qf, qp);
              const Real weight = qf.getWeight(qp) * p.getDistortion();
              volume += weight;
              moment += weight * (uh(ip) - exact(p));
            }
            EXPECT_GT(volume, 0);
            momentSquared += ErrorNorm::squaredMagnitude(moment) / volume;
          }
#ifdef RODIN_USE_MPI
        if constexpr (requires { m_mesh.getShard(); })
          momentSquared = boost::mpi::all_reduce(
            m_mesh.getContext().getCommunicator(), momentSquared, std::plus<Real>());
#endif
        const Real error = ErrorNorm::computeL2(m_mesh, uh, exact, OracleOrder);
        Real meanError = 0;
        if constexpr (Global)
          meanError = ErrorNorm::computeL2(
            m_mesh, uh, field<Scalar, Vector>(dim, constant, true), OracleOrder);
        return {error, std::sqrt(momentSquared), meanError};
      }

    private:
      template <class Scalar>
      static Scalar coefficient(Real real, Real imaginary)
      {
        if constexpr (std::is_same_v<Scalar, Real>)
          return real;
        else
          return Scalar(real, imaginary);
      }
      template <class Scalar, bool Vector>
      static auto field(size_t dim, bool constant, bool mean = false)
      {
        const auto callback = [dim, constant, mean](const Point& p) {
          const Real x = mean ? (dim == 1 ? Real(0.55) : Real(0.5)) : p(0);
          const Real y =
            mean ? (dim == 1 ? Real(0.55) : Real(0.5) + Real(0.1) / 3) : p(dim - 1);
          if constexpr (Vector)
          {
            Math::SpatialVector<Scalar> value(2);
            value(0) = coefficient<Scalar>(1.25, 0.5);
            value(1) = coefficient<Scalar>(-0.75, 0.25);
            if (!constant)
            {
              value(0) += coefficient<Scalar>(1, -0.75) * x;
              value(1) += coefficient<Scalar>(0.8, 0.2) * y;
            }
            return value;
          }
          else
            return coefficient<Scalar>(1.25, 0.5) +
              (constant ? Scalar(0)
                        : coefficient<Scalar>(1, -0.75) * x +
                    coefficient<Scalar>(0.3, 0.2) * y);
        };
        if constexpr (Vector)
          return VectorFunction<decltype(callback)>(2, callback);
        else if constexpr (std::is_same_v<Scalar, Real>)
          return RealFunction(callback);
        else
          return ComplexFunction(callback);
      }
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_CURVED_P0_PETSC) && defined(RODIN_USE_MPI)
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
            .makeMesh(n);
#endif
      }
      Mesh<ContextType> m_mesh;
      CurvedGeometry<Mesh<ContextType>> m_geometry;
  };
  template <class ContextType>
  class ProjectionTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      // Mesh geometry is common to the value cases. Fields, spaces and systems
      // remain fresh per projection; no coefficient or solver state is reused.
      template <class Callback>
      static void values(Callback&& callback)
      {
#if !defined(RODIN_CURVED_P0_PETSC) || !defined(PETSC_USE_COMPLEX)
        {
          SCOPED_TRACE("real scalar");
          callback.template operator()<Real, false>();
        }
        {
          SCOPED_TRACE("real vector");
          callback.template operator()<Real, true>();
        }
#endif
#if !defined(RODIN_CURVED_P0_PETSC) || defined(PETSC_USE_COMPLEX)
        {
          SCOPED_TRACE("complex scalar");
          callback.template operator()<Complex, false>();
        }
        {
          SCOPED_TRACE("complex vector");
          callback.template operator()<Complex, true>();
        }
#endif
      }
      void rates() const
      {
        std::array<NormHistory, 4> histories;
        for (size_t n : {5u, 9u, 17u})
        {
          SCOPED_TRACE(::testing::Message() << "n=" << n);
          Workload<ContextType> problem(this->GetParam(), n);
          size_t sample = 0;
          values([&]<class Scalar, bool Vector>() {
            const auto e = problem.template project<Scalar, Vector>();
            EXPECT_LT(e.moment, MomentTolerance);
            histories[sample++].append(Real(1) / Real(n - 1), e.l2);
          });
        }
        size_t sample = 0;
        constexpr Real RateMargin = 0.2;
        values([&]<class Scalar, bool Vector>() {
          const auto& history = histories[sample++];
          for (size_t i = 1; i < history.getSize(); ++i)
          {
            const Real coarse = history.getSample(i - 1).error,
                       fine = history.getSample(i).error;
            ASSERT_TRUE(std::isfinite(coarse));
            ASSERT_TRUE(std::isfinite(fine));
            ASSERT_GT(fine, 0);
            EXPECT_GT(coarse, fine);
            const Real rate = history.getAlgebraicRate(i);
            SCOPED_TRACE(::testing::Message() << "interval=" << i << " rate=" << rate);
            EXPECT_GT(rate, 1 - RateMargin);
            EXPECT_LT(rate, 1 + RateMargin);
          }
        });
      }
      void constants() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        values([&]<class Scalar, bool Vector>() {
          const auto local = problem.template project<Scalar, Vector>(true);
          const auto global = problem.template project<Scalar, Vector, true>(true);
          EXPECT_LT(local.l2, PatchTolerance);
          EXPECT_LT(local.moment, MomentTolerance);
          EXPECT_LT(global.l2, PatchTolerance);
          EXPECT_LT(global.globalMean, PatchTolerance);
        });
      }
      void globalMean() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        values([&]<class Scalar, bool Vector>() {
          const auto e = problem.template project<Scalar, Vector, true>();
          EXPECT_LT(e.globalMean, PatchTolerance);
          EXPECT_GT(e.l2, 0);
        });
      }
      void controls() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        values([&]<class Scalar, bool Vector>() {
          const auto base = problem.template project<Scalar, Vector>();
          const auto centroid = problem.template project<Scalar, Vector>(
            false, AssemblyOrder, SolverTolerance, true);
          EXPECT_LT(base.moment, MomentTolerance);
          EXPECT_GT(centroid.moment, CentroidDefectMinimum);
          const auto wrong = problem.template project<Scalar, Vector>(
            true, AssemblyOrder, SolverTolerance, false, 2);
          EXPECT_GT(wrong.l2, WrongMassErrorMinimum);
        });
      }
      void sensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        values([&]<class Scalar, bool Vector>() {
          const auto base = problem.template project<Scalar, Vector>();
          const auto quad = problem.template project<Scalar, Vector>(false, RefinedOrder);
          const auto solver = problem.template project<Scalar, Vector>(
            false, AssemblyOrder, RefinedTolerance);
          for (const auto& e : {quad, solver})
          {
            ASSERT_TRUE(std::isfinite(e.l2));
            ASSERT_GT(base.l2, 0);
            EXPECT_LT(std::abs(e.l2 / base.l2 - 1), SensitivityTolerance);
            EXPECT_LT(e.moment, MomentTolerance);
          }
        });
      }
  };
  using LocalTest = ProjectionTest<Context::Local>;
  TEST_P(LocalTest, AllValueTypesRates)
  {
    rates();
  }
  TEST_P(LocalTest, AllValueTypesConstants)
  {
    constants();
  }
  TEST_P(LocalTest, AllValueTypesGlobalMeans)
  {
    globalMean();
  }
  TEST_P(LocalTest, AllValueTypesControls)
  {
    controls();
  }
  TEST_P(LocalTest, AllValueTypesSensitivity)
  {
    sensitivity();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#if defined(RODIN_CURVED_P0_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = ProjectionTest<Context::MPI>;
  TEST_P(MPITest, AllValueTypesRates)
  {
    rates();
  }
  TEST_P(MPITest, AllValueTypesConstants)
  {
    constants();
  }
  TEST_P(MPITest, AllValueTypesGlobalMeans)
  {
    globalMean();
  }
  TEST_P(MPITest, AllValueTypesControls)
  {
    controls();
  }
  TEST_P(MPITest, AllValueTypesSensitivity)
  {
    sensitivity();
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

#ifdef RODIN_CURVED_P0_PETSC
int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::P0Projection::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::P0Projection::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
