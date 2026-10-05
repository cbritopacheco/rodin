/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <array>
#include <gtest/gtest.h>
#include "../../CurvedGeometry.h"
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"
#ifdef RODIN_GEOMETRY_APPROXIMATION_PETSC
#include "Rodin/PETSc.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif
#endif
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace Rodin::Tests::Convergence::Isoparametric::GeometryApproximation
{
  constexpr Real Amplitude = 0.1;
  constexpr size_t MapOrder = 14, RefinedMapOrder = 18;
  constexpr size_t AssemblyOrder = 12, NormOrder = 14, RefinedOrder = 16;
  constexpr Real SolverTolerance = 1e-13, RefinedTolerance = 1e-14;
  constexpr Real PatchTolerance = 1e-9, ResidualTolerance = 1e-11;
  constexpr Real SensitivityTolerance = 1e-6, RateMargin = 0.35;
  constexpr Real WrongMapValueMinimum = 0.05, WrongMapDerivativeMinimum = 0.15;
  constexpr Real WrongTraceMinimum = 0.1;
  // Finite-path improvement policies, not fixed-degree algebraic rate bounds.
  constexpr Real DegreeReduction = 0.5, CoupledReduction = 0.25;
  constexpr size_t MaxIterations = 50000;
  [[maybe_unused]] constexpr Real DivergenceTolerance = 1e5;
#if defined(RODIN_GEOMETRY_APPROXIMATION_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  /** @brief Reference-domain geometry errors and physical-domain Poisson patches.
   * @par Architecture
   * A mesh copy freezes the original charts and exact logical entity indices.
   * Geometry is interpolated only on the second mesh. Map and global-coordinate
   * derivative errors are integrated on the original unit box, with analytic
   * sine-map oracles independent of the installation helper. The matched-degree
   * affine Poisson patch measures field error on the approximated domain.
   * MPI sums owned reference cells only; no coordinate matching is used.
   */
  template <class ContextType, size_t Q>
  class Workload
  {
    public:
      Workload(Polytope::Type geometry, size_t n, bool wrongMap = false)
        : m_dimension(UniformGrid(geometry).getDimension()),
          m_mesh(makeMesh(geometry, n)),
          m_reference(m_mesh),
          m_mapping(m_mesh, CurvedGeometry<Mesh<ContextType>>::Map::Sine,
            wrongMap ? Real(0) : Amplitude)
      {
        m_mapping.template install<Q>();
      }
      ErrorNorms mapErrors(size_t order = MapOrder) const
      {
        Real valueSquared = 0, derivativeSquared = 0;
        const Real pi = Math::Constants::pi();
        for (auto referenceCell = m_reference.getCell(); referenceCell; ++referenceCell)
        {
          if constexpr (requires { m_reference.getShard(); })
            if (!m_reference.getShard().isOwned(m_dimension, referenceCell->getIndex()))
              continue;
          const auto cell = m_mesh.getCell(referenceCell->getIndex());
          EXPECT_EQ(cell->getGeometry(), referenceCell->getGeometry());
          const auto vertices = cell->getVertices();
          const auto originalVertices = referenceCell->getVertices();
          EXPECT_EQ(vertices.size(), originalVertices.size());
          for (size_t local = 0; local < vertices.size(); ++local)
            EXPECT_EQ(vertices[local], originalVertices[local]);
          const auto& qf =
            QF::PolytopeQuadratureFormula::get(order, referenceCell->getGeometry());
          const Geometry::PolytopeQuadrature quadrature(*referenceCell, qf);
          for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
          {
            const auto& referencePoint = quadrature.getPoint(qp);
            const Point mapped(*cell, referencePoint.getReferenceCoordinates());
            auto exact = referencePoint.getPhysicalCoordinates();
            exact(m_dimension - 1) += Amplitude * std::sin(pi * exact(0));
            Math::SpatialMatrix<Real> derivative =
              Math::SpatialMatrix<Real>::Identity(m_dimension, m_dimension);
            derivative(m_dimension - 1, 0) +=
              Amplitude * pi * std::cos(pi * referencePoint(0));
            const auto defect =
              mapped.getJacobian() * referencePoint.getJacobianInverse() - derivative;
            const Real weight = qf.getWeight(qp) * referencePoint.getDistortion();
            valueSquared += weight *
              ErrorNorm::squaredMagnitude(mapped.getPhysicalCoordinates() - exact);
            derivativeSquared += weight * ErrorNorm::squaredMagnitude(defect);
            EXPECT_TRUE(std::isfinite(mapped.getJacobianDeterminant()));
            EXPECT_GT(
              mapped.getJacobianDeterminant() / referencePoint.getJacobianDeterminant(),
              0);
          }
        }
#ifdef RODIN_USE_MPI
        if constexpr (requires { m_reference.getShard(); })
        {
          const auto& comm = m_reference.getContext().getCommunicator();
          valueSquared = boost::mpi::all_reduce(comm, valueSquared, std::plus<Real>());
          derivativeSquared =
            boost::mpi::all_reduce(comm, derivativeSquared, std::plus<Real>());
        }
#endif
        return {std::sqrt(valueSquared), std::sqrt(derivativeSquared)};
      }
      ErrorNorms patch(bool wrongTrace = false, size_t order = AssemblyOrder,
        Real tolerance = SolverTolerance) const
      {
        const auto exact = RealFunction([dim = m_dimension](const Point& p) {
          Real value = 1;
          for (size_t j = 0; j < dim; ++j)
            value += p(j);
          return value;
        });
        const auto gradient =
          VectorFunction(m_dimension, [dim = m_dimension](const Point&) {
            Math::SpatialVector<Real> value(dim);
            value.setConstant(1);
            return value;
          });
        const auto trace = RealFunction([exact, wrongTrace](const Point& p) {
          return wrongTrace ? Real(0) : exact(p);
        });
        H1<Q, Real, Mesh<ContextType>> space(std::integral_constant<size_t, Q>{}, m_mesh);
#ifdef RODIN_GEOMETRY_APPROXIMATION_PETSC
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
#else
        TrialFunction u(space);
        TestFunction v(space);
#endif
        auto stiffness = Integral(Grad(u), Grad(v));
        stiffness.setOrder(order);
        Problem problem(u, v);
        problem = stiffness + DirichletBC(u, trace);
#ifdef RODIN_GEOMETRY_APPROXIMATION_PETSC
        PETSc::Solver::CG solver(problem);
        solver
          .setTolerances(tolerance, RefinedTolerance, DivergenceTolerance, MaxIterations)
          .solve();
        KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
        EXPECT_EQ(KSPGetConvergedReason(solver.getHandle(), &reason), PETSC_SUCCESS);
        // A zero wrong-trace state may terminate without a Krylov iteration.
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
        return ErrorNorm::compute(m_mesh, u.getSolution(), exact, gradient, NormOrder);
      }

    private:
      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_GEOMETRY_APPROXIMATION_PETSC) && defined(RODIN_USE_MPI)
        else
          return DistributedUniformGrid(Context::MPI(*environment, *world), geometry)
            .makeMesh(n);
#endif
      }
      size_t m_dimension;
      Mesh<ContextType> m_mesh, m_reference;
      CurvedGeometry<Mesh<ContextType>> m_mapping;
  };

  template <class ContextType>
  class ApproximationTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      template <bool RefineMesh>
      void degreePath() const
      {
        const std::array<size_t, 3> levels =
          RefineMesh ? std::array<size_t, 3>{2, 3, 5} : std::array<size_t, 3>{3, 3, 3};
        Workload<ContextType, 1> first(this->GetParam(), levels[0]);
        Workload<ContextType, 2> second(this->GetParam(), levels[1]);
        Workload<ContextType, 3> third(this->GetParam(), levels[2]);
        const std::array errors{first.mapErrors(), second.mapErrors(), third.mapErrors()};
        const std::array refined{first.mapErrors(RefinedMapOrder),
          second.mapErrors(RefinedMapOrder), third.mapErrors(RefinedMapOrder)};
        const std::array patches{first.patch(), second.patch(), third.patch()};
        for (size_t i = 0; i < errors.size(); ++i)
        {
          SCOPED_TRACE(
            ::testing::Message() << "geometry degree=" << i + 1 << " n=" << levels[i]);
          ASSERT_TRUE(errors[i].isFinite());
          ASSERT_GT(errors[i].getL2(), 0);
          ASSERT_GT(errors[i].getH1Seminorm(), 0);
          ASSERT_TRUE(refined[i].isFinite());
          EXPECT_LT(
            std::abs(refined[i].getL2() / errors[i].getL2() - 1), SensitivityTolerance);
          EXPECT_LT(std::abs(refined[i].getH1Seminorm() / errors[i].getH1Seminorm() - 1),
            SensitivityTolerance);
          ASSERT_TRUE(patches[i].isFinite());
          EXPECT_LT(patches[i].getL2(), PatchTolerance);
          EXPECT_LT(patches[i].getH1Seminorm(), PatchTolerance);
          if (i == 0)
            continue;
          SCOPED_TRACE(::testing::Message()
            << "map errors " << errors[i - 1].getL2() << " -> " << errors[i].getL2()
            << " derivative " << errors[i - 1].getH1Seminorm() << " -> "
            << errors[i].getH1Seminorm());
          const Real reduction = RefineMesh ? CoupledReduction : DegreeReduction;
          EXPECT_LT(errors[i].getL2(), reduction * errors[i - 1].getL2());
          EXPECT_LT(errors[i].getH1Seminorm(), reduction * errors[i - 1].getH1Seminorm());
        }
      }

      template <size_t Q>
      void rates() const
      {
        ErrorHistory history;
        for (size_t n : {3u, 5u, 9u})
        {
          Workload<ContextType, Q> problem(this->GetParam(), n);
          history.append(Real(1) / Real(n - 1), problem.mapErrors());
        }
        for (size_t i = 1; i < history.getSize(); ++i)
        {
          const auto& coarse = history.getSample(i - 1).error;
          const auto& fine = history.getSample(i).error;
          const auto rate = history.getAlgebraicRates(i);
          SCOPED_TRACE(::testing::Message()
            << "q=" << Q << " interval=" << i << " errors " << coarse.getL2() << " -> "
            << fine.getL2() << " derivative " << coarse.getH1Seminorm() << " -> "
            << fine.getH1Seminorm() << " rates=" << rate.getL2() << ","
            << rate.getH1Seminorm());
          ASSERT_TRUE(coarse.isFinite());
          ASSERT_TRUE(fine.isFinite());
          ASSERT_GT(fine.getL2(), 0);
          ASSERT_GT(fine.getH1Seminorm(), 0);
          EXPECT_GT(coarse.getL2(), fine.getL2());
          EXPECT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
          EXPECT_GT(rate.getL2(), Q + 1 - RateMargin);
          EXPECT_LT(rate.getL2(), Q + 1 + RateMargin);
          EXPECT_GT(rate.getH1Seminorm(), Q - RateMargin);
          EXPECT_LT(rate.getH1Seminorm(), Q + RateMargin);
        }
      }
      template <size_t Q>
      void patch() const
      {
        Workload<ContextType, Q> problem(this->GetParam(), 3);
        const auto map = problem.mapErrors();
        const auto field = problem.patch();
        EXPECT_GT(map.getL2(), 0);
        EXPECT_GT(map.getH1Seminorm(), 0);
        EXPECT_LT(field.getL2(), PatchTolerance);
        EXPECT_LT(field.getH1Seminorm(), PatchTolerance);
      }
      void controls() const
      {
        Workload<ContextType, 3> wrongMap(this->GetParam(), 5, true);
        const auto map = wrongMap.mapErrors();
        EXPECT_GT(map.getL2(), WrongMapValueMinimum);
        EXPECT_GT(map.getH1Seminorm(), WrongMapDerivativeMinimum);
        // An exact physical affine patch does not detect a wrong domain map.
        const auto field = wrongMap.patch();
        EXPECT_LT(field.getL2(), PatchTolerance);
        EXPECT_LT(field.getH1Seminorm(), PatchTolerance);
        Workload<ContextType, 3> correctMap(this->GetParam(), 3);
        const auto wrongField = correctMap.patch(true);
        EXPECT_GT(wrongField.getL2(), WrongTraceMinimum);
        EXPECT_GT(wrongField.getH1Seminorm(), WrongTraceMinimum);
      }
      template <size_t Q>
      void sensitivity() const
      {
        Workload<ContextType, Q> problem(this->GetParam(), 3);
        const auto base = problem.mapErrors(),
                   refined = problem.mapErrors(RefinedMapOrder);
        for (const auto pair : {std::pair{base.getL2(), refined.getL2()},
               std::pair{base.getH1Seminorm(), refined.getH1Seminorm()}})
        {
          ASSERT_GT(pair.first, 0);
          EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
        }
        for (const auto& field : {problem.patch(), problem.patch(false, RefinedOrder),
               problem.patch(false, AssemblyOrder, RefinedTolerance)})
        {
          EXPECT_LT(field.getL2(), PatchTolerance);
          EXPECT_LT(field.getH1Seminorm(), PatchTolerance);
        }
      }
  };
  using LocalTest = ApproximationTest<Context::Local>;
  TEST_P(LocalTest, GeometryDegreeRefinement)
  {
    degreePath<false>();
  }
  TEST_P(LocalTest, CoupledMeshAndGeometryRefinement)
  {
    degreePath<true>();
  }
  TEST_P(LocalTest, Q1Rates)
  {
    rates<1>();
  }
  TEST_P(LocalTest, Q2Rates)
  {
    rates<2>();
  }
  TEST_P(LocalTest, Q3Rates)
  {
    rates<3>();
  }
  TEST_P(LocalTest, Q1Patch)
  {
    patch<1>();
  }
  TEST_P(LocalTest, Q2Patch)
  {
    patch<2>();
  }
  TEST_P(LocalTest, Q3Patch)
  {
    patch<3>();
  }
  TEST_P(LocalTest, IndependentControls)
  {
    controls();
  }
  TEST_P(LocalTest, Q1Sensitivity)
  {
    sensitivity<1>();
  }
  TEST_P(LocalTest, Q2Sensitivity)
  {
    sensitivity<2>();
  }
  TEST_P(LocalTest, Q3Sensitivity)
  {
    sensitivity<3>();
  }
  INSTANTIATE_TEST_SUITE_P(AllGeometries, LocalTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge),
    [](const auto& info) {
      return std::string(UniformGrid::getGeometryName(info.param));
    });
#if defined(RODIN_GEOMETRY_APPROXIMATION_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = ApproximationTest<Context::MPI>;
  TEST_P(MPITest, GeometryDegreeRefinement)
  {
    degreePath<false>();
  }
  TEST_P(MPITest, CoupledMeshAndGeometryRefinement)
  {
    degreePath<true>();
  }
  TEST_P(MPITest, Q1Rates)
  {
    rates<1>();
  }
  TEST_P(MPITest, Q2Rates)
  {
    rates<2>();
  }
  TEST_P(MPITest, Q3Rates)
  {
    rates<3>();
  }
  TEST_P(MPITest, Q1Patch)
  {
    patch<1>();
  }
  TEST_P(MPITest, Q2Patch)
  {
    patch<2>();
  }
  TEST_P(MPITest, Q3Patch)
  {
    patch<3>();
  }
  TEST_P(MPITest, IndependentControls)
  {
    controls();
  }
  TEST_P(MPITest, Q1Sensitivity)
  {
    sensitivity<1>();
  }
  TEST_P(MPITest, Q2Sensitivity)
  {
    sensitivity<2>();
  }
  TEST_P(MPITest, Q3Sensitivity)
  {
    sensitivity<3>();
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
#ifdef RODIN_GEOMETRY_APPROXIMATION_PETSC
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
    Rodin::Tests::Convergence::Isoparametric::GeometryApproximation::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::GeometryApproximation::world = &comm;
#endif
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
