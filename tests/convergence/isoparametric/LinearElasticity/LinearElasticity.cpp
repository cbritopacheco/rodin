/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Curved elasticity convergence and independent vector-lift metric oracles.
 *
 * The physical problem is @f$-\operatorname{div}\sigma(u)=f@f$ with full
 * Dirichlet trace and @f$\sigma(u)=\lambda\operatorname{div}(u)I+
 * 2\mu\varepsilon(u)@f$. Geometry is installed independently of field degree.
 * P1 constants and physical affine P2 fields are representable patches;
 * the exponential field tests approximation rates rather than reproduction.
 */

#include <array>

#include "../../LinearElasticity.h"
#include "../../CurvedGeometry.h"
#include "../../LiftedErrorNorm.h"
#include "../../SineMap.h"

#ifdef RODIN_CURVED_ELASTICITY_PETSC
#include "Rodin/PETSc.h"
#ifdef RODIN_USE_MPI
#include "MPIConvergence.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#endif
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Convergence::Isoparametric::LinearElasticity
{
  using Data = Rodin::Tests::Convergence::LinearElasticity::ManufacturedSolution;
  // Dimensionless Lamé coefficients remain fixed through every study.
  constexpr Real Lambda = 1.5;
  constexpr Real Mu = 0.5;
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
  constexpr Real NegativeControlFactor = 2;

#if defined(RODIN_CURVED_ELASTICITY_PETSC) && defined(RODIN_USE_MPI)
  boost::mpi::environment* environment = nullptr;
  boost::mpi::communicator* world = nullptr;
#endif

  struct Errors
  {
      ErrorNorms displacement;
      Real strain;
      Real stress;
  };

  struct LiftedErrors
  {
      LiftedErrorNorm::Result displacement;
      std::array<Real, 3> strain{}, stress{};
  };

  /**
   * @brief A mapped elasticity workload with backend-independent observables.
   * @par Architecture
   * A fresh mesh owns its P2 transformations, exact for the quadratic map or
   * interpolated for the sine map. Each solve creates fresh
   * vector spaces and a linear system. Shared manufactured physical data are
   * consumed by native/PETSc assembly; the solver policy remains explicit.
   * ErrorNorm integrates physical quantities on owned cells and reduces MPI
   * squared norms once. No coefficient-layout assumptions enter the tests.
   */
  template <class ContextType>
  class Workload
  {
    public:
      using Map = typename CurvedGeometry<Mesh<ContextType>>::Map;
      Workload(Polytope::Type geometry, size_t n, Map map = Map::Quadratic,
        bool lifted = false, Real amplitude = 0.1)
        : m_mesh(makeMesh(geometry, n)),
          m_geometry(m_mesh, map, amplitude)
      {
        if (lifted)
          m_reference.emplace(m_mesh);
        m_geometry.template install<2>();
      }

      const auto& getMesh() const
      {
        return m_mesh;
      }
      const auto& getReference() const
      {
        return m_reference.value();
      }

      template <size_t K>
      Errors solve(Data::Field field, bool omitVolumetric = false,
        size_t order = AssemblyOrder, Real tolerance = SolverTolerance,
        size_t normOrder = 0, LiftedErrors* lifted = nullptr) const
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Data data(dim, Lambda, Mu, field);
        auto space = [&] {
#ifndef RODIN_CURVED_ELASTICITY_PETSC
          if constexpr (K == 1)
            return P1(m_mesh, dim);
          else
#endif
            return H1<K, Math::SpatialVector<Real>, Mesh<ContextType>>(
              std::integral_constant<size_t, K>{}, m_mesh, dim);
        }();
#ifdef RODIN_CURVED_ELASTICITY_PETSC
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
#else
        TrialFunction u(space);
        TestFunction v(space);
#endif
        auto volumetric = Integral((omitVolumetric ? Real(0) : Lambda) * Div(u), Div(v));
        auto shear = Integral(
          Mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
        auto body = Integral(data.forcing, v);
        volumetric.setOrder(order);
        shear.setOrder(order);
        body.setOrder(order);
        Problem problem(u, v);
        problem = volumetric + shear - body + DirichletBC(u, data.exact);
#ifdef RODIN_CURVED_ELASTICITY_PETSC
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
        const auto jacobian = Jacobian(u.getSolution());
        const auto exactJacobian = [&data](const Point& p) { return data.jacobian(p); };
        const auto strain = [&jacobian](const IntegrationPoint& ip) {
          return tensor(jacobian(ip), false);
        };
        const auto stress = [&jacobian](const IntegrationPoint& ip) {
          return tensor(jacobian(ip), true);
        };
        const size_t integrationOrder = normOrder == 0 ? order + 2 : normOrder;
        if (lifted)
        {
          assert(m_reference);
          std::array<Real, 6> squared{};
          // Constitutive norms consume the same owned-cell quadrature traversal.
          // Applying the linear law to derivative defects gives tensor defects.
          const auto observe =
            [&squared](
              const std::array<Math::SpatialMatrix<Real>, 3>& derivatives, Real weight) {
              for (size_t i = 0; i < derivatives.size(); ++i)
              {
                squared[2 * i] +=
                  weight * ErrorNorm::squaredMagnitude(tensor(derivatives[i], false));
                squared[2 * i + 1] +=
                  weight * ErrorNorm::squaredMagnitude(tensor(derivatives[i], true));
              }
            };
          lifted->displacement = LiftedErrorNorm::compute(*m_reference, m_mesh,
            u.getSolution(), data, SineMap(), integrationOrder, observe);
#ifdef RODIN_USE_MPI
          if constexpr (requires { m_mesh.getShard(); })
            for (Real& value : squared)
              value = boost::mpi::all_reduce(
                m_mesh.getContext().getCommunicator(), value, std::plus<Real>());
#endif
          for (size_t i = 0; i < lifted->strain.size(); ++i)
          {
            lifted->strain[i] = std::sqrt(squared[2 * i]);
            lifted->stress[i] = std::sqrt(squared[2 * i + 1]);
          }
        }
        return {ErrorNorm::computeVector(
                  m_mesh, u.getSolution(), data.exact, exactJacobian, integrationOrder),
          ErrorNorm::computeL2(
            m_mesh, strain, [&data](const Point& p) { return data.strain(p); },
            integrationOrder),
          ErrorNorm::computeL2(
            m_mesh, stress, [&data](const Point& p) { return data.stress(p); },
            integrationOrder)};
      }

    private:
      static auto tensor(const Math::SpatialMatrix<Real>& jacobian, bool stress)
      {
        auto value = jacobian;
        Real trace = 0;
        for (size_t i = 0; i < jacobian.rows(); ++i)
          trace += jacobian(i, i);
        for (size_t i = 0; i < jacobian.rows(); ++i)
          for (size_t j = 0; j < jacobian.cols(); ++j)
            value(i, j) = stress
              ? Lambda * trace * (i == j) + Mu * (jacobian(i, j) + jacobian(j, i))
              : (jacobian(i, j) + jacobian(j, i)) / 2;
        return value;
      }

      static Mesh<ContextType> makeMesh(Polytope::Type geometry, size_t n)
      {
        if constexpr (std::is_same_v<ContextType, Context::Local>)
          return UniformGrid(geometry).makeMesh(n);
#if defined(RODIN_CURVED_ELASTICITY_PETSC) && defined(RODIN_USE_MPI)
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
  class CurvedLinearElasticityTest : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      void checkComplexVectorMetric() const
      {
        using Map = typename Workload<ContextType>::Map;
        Workload<ContextType> problem(this->GetParam(), 3, Map::Sine, true, 0);
        const size_t dim = problem.getMesh().getSpaceDimension();
        struct ComplexData
        {
            Data real;
            Math::SpatialVector<Complex> getSolution(const Math::SpatialPoint& x) const
            {
              return Complex(1, 1) * real.getSolution(x);
            }
            Math::SpatialMatrix<Complex> getJacobian(const Math::SpatialPoint& x) const
            {
              return Complex(1, 1) * real.getJacobian(x);
            }
        };
        const ComplexData data{Data(dim, Lambda, Mu, Data::Field::AsymmetricAffine)};
        H1<2, Math::SpatialVector<Complex>, Mesh<ContextType>> space(
          std::integral_constant<size_t, 2>{}, problem.getMesh(), dim);
        GridFunction u(space);
        u = VectorFunction(dim, [data](const Point& p) {
          return data.getSolution(p.getPhysicalCoordinates());
        });
        constexpr size_t MetricOrder = 18;
        constexpr Real Amplitude = 0.1;
        const auto e = LiftedErrorNorm::compute(
          problem.getReference(), problem.getMesh(), u, data, SineMap(), MetricOrder);
        Real columnSquared = 0;
        for (size_t i = 0; i < dim; ++i)
        {
          const Real coefficient = Real(dim * (i + 1)) + (i == 0);
          columnSquared += coefficient * coefficient;
        }
        const Real slope = Amplitude * Math::Constants::pi();
        const Real l2 = Amplitude * std::sqrt(columnSquared);
        const Real h1 = std::sqrt(2 * columnSquared) *
          (dim == 1 ? std::sqrt(Real(1) / std::sqrt(Real(1) - slope * slope) - 1)
                    : slope / std::sqrt(Real(2)));
        EXPECT_TRUE(e.field.isFinite());
        EXPECT_LT(e.field.getL2(), PatchTolerance);
        EXPECT_LT(e.field.getH1Seminorm(), PatchTolerance);
        for (const auto& norm : {e.geometry, e.total})
        {
          EXPECT_TRUE(norm.isFinite());
          EXPECT_NEAR(norm.getL2(), l2, PatchTolerance);
          EXPECT_NEAR(norm.getH1Seminorm(), h1, PatchTolerance);
        }
      }

      void checkLiftedMetric() const
      {
        using Map = typename Workload<ContextType>::Map;
        Workload<ContextType> problem(this->GetParam(), 3, Map::Sine, true, 0);
        LiftedErrors error;
        constexpr size_t MetricOrder = 18;
        const auto represented = problem.template solve<2>(Data::Field::AsymmetricAffine,
          false, AssemblyOrder, SolverTolerance, MetricOrder, &error);
        const size_t dim = problem.getMesh().getSpaceDimension();
        constexpr Real Amplitude = 0.1;
        const Real slope = Amplitude * Math::Constants::pi();
        // Last column of A_ij=(i+1)(j+1)+delta_i0 delta_j,d-1.
        // Its norm differs from the first column and the last row in d>=2.
        Real columnSquared = 0;
        for (size_t i = 0; i < dim; ++i)
        {
          const Real coefficient = Real(dim * (i + 1)) + (i == 0);
          columnSquared += coefficient * coefficient;
        }
        const Real first = Real(dim + 1);
        const Real l2 = Amplitude * std::sqrt(columnSquared / 2);
        const Real h1 = std::sqrt(columnSquared) *
          (dim == 1 ? std::sqrt(Real(1) / std::sqrt(Real(1) - slope * slope) - 1)
                    : slope / std::sqrt(Real(2)));
        const Real strain =
          dim == 1 ? h1 : slope * std::sqrt((columnSquared + first * first) / 4);
        const Real stress = dim == 1 ? (Lambda + 2 * Mu) * h1
                                     : slope / std::sqrt(Real(2)) *
            std::sqrt(2 * Mu * Mu * columnSquared +
              (Real(dim) * Lambda * Lambda + 4 * Lambda * Mu + 2 * Mu * Mu) * first *
                first);
        for (const auto& e : {represented.displacement, error.displacement.field})
        {
          EXPECT_TRUE(e.isFinite());
          EXPECT_LT(e.getL2(), PatchTolerance);
          EXPECT_LT(e.getH1Seminorm(), PatchTolerance);
        }
        EXPECT_LT(represented.strain, PatchTolerance);
        EXPECT_LT(represented.stress, PatchTolerance);
        EXPECT_LT(error.strain[0], PatchTolerance);
        EXPECT_LT(error.stress[0], PatchTolerance);
        for (const auto& e : {error.displacement.geometry, error.displacement.total})
        {
          EXPECT_TRUE(e.isFinite());
          EXPECT_NEAR(e.getL2(), l2, PatchTolerance);
          EXPECT_NEAR(e.getH1Seminorm(), h1, PatchTolerance);
        }
        for (size_t i : {1u, 2u})
        {
          EXPECT_NEAR(error.strain[i], strain, PatchTolerance);
          EXPECT_NEAR(error.stress[i], stress, PatchTolerance);
        }
        const auto& reference = problem.getReference();
        for (auto cell = problem.getMesh().getCell(); cell; ++cell)
        {
          const auto original = reference.getCell(cell->getIndex());
          EXPECT_EQ(cell->getGeometry(), original->getGeometry());
          const auto vertices = cell->getVertices(),
                     originalVertices = original->getVertices();
          ASSERT_EQ(vertices.size(), originalVertices.size());
          for (size_t i = 0; i < vertices.size(); ++i)
            EXPECT_EQ(vertices[i], originalVertices[i]);
          if constexpr (requires { reference.getShard(); })
            EXPECT_EQ(reference.getShard().isOwned(dim, cell->getIndex()),
              problem.getMesh().getShard().isOwned(dim, cell->getIndex()));
        }
      }

      template <size_t K>
      void checkRates() const
      {
        ErrorHistory displacement;
        NormHistory strain, stress;
        const auto levels = K == 1 ? std::initializer_list<size_t>{5, 9, 17}
                                   : std::initializer_list<size_t>{3, 5, 9};
        for (size_t n : levels)
        {
          SCOPED_TRACE(::testing::Message() << "degree=" << K << " n=" << n);
          Workload<ContextType> problem(this->GetParam(), n);
          const auto error = problem.template solve<K>(Data::Field::Exponential);
          const Real h = Real(1) / Real(n - 1);
          displacement.append(h, error.displacement);
          strain.append(h, error.strain);
          stress.append(h, error.stress);
        }
        // Case-specific finite-resolution windows around the theoretical orders.
        constexpr Real L2Margin = 0.55;
        constexpr Real DerivativeMargin = 0.45;
        for (size_t i = 1; i < displacement.getSize(); ++i)
        {
          const auto& coarse = displacement.getSample(i - 1).error;
          const auto& fine = displacement.getSample(i).error;
          ASSERT_TRUE(coarse.isFinite());
          ASSERT_TRUE(fine.isFinite());
          ASSERT_GT(fine.getL2(), 0);
          ASSERT_GT(fine.getH1Seminorm(), 0);
          ASSERT_GT(coarse.getL2(), fine.getL2());
          ASSERT_GT(coarse.getH1Seminorm(), fine.getH1Seminorm());
          const auto rate = displacement.getAlgebraicRates(i);
          SCOPED_TRACE(::testing::Message()
            << "interval=" << i << " L2 rate=" << rate.getL2()
            << " H1 rate=" << rate.getH1Seminorm());
          EXPECT_GT(rate.getL2(), K + 1 - L2Margin);
          EXPECT_LT(rate.getL2(), K + 1 + L2Margin);
          EXPECT_GT(rate.getH1Seminorm(), K - DerivativeMargin);
          EXPECT_LT(rate.getH1Seminorm(), K + DerivativeMargin);
          for (const auto* history : {&strain, &stress})
          {
            const Real coarseError = history->getSample(i - 1).error;
            const Real fineError = history->getSample(i).error;
            ASSERT_TRUE(std::isfinite(coarseError));
            ASSERT_TRUE(std::isfinite(fineError));
            ASSERT_GT(fineError, 0);
            ASSERT_GT(coarseError, fineError);
            const Real tensorRate = history->getAlgebraicRate(i);
            SCOPED_TRACE(::testing::Message() << "tensor errors=" << coarseError << " -> "
                                              << fineError << " rate=" << tensorRate);
            EXPECT_GT(tensorRate, K - DerivativeMargin);
            EXPECT_LT(tensorRate, K + DerivativeMargin);
          }
        }
      }

      template <size_t K>
      void checkPatch() const
      {
        Workload<ContextType> problem(this->GetParam(), 3);
        const auto error =
          problem.template solve<K>(K == 1 ? Data::Field::Constant : Data::Field::Affine);
        EXPECT_LT(error.displacement.getL2(), PatchTolerance);
        EXPECT_LT(error.displacement.getH1Seminorm(), PatchTolerance);
        EXPECT_LT(error.strain, PatchTolerance);
        EXPECT_LT(error.stress, PatchTolerance);
      }

      void checkControl() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        const auto baseline = problem.template solve<2>(Data::Field::Exponential);
        const auto wrong = problem.template solve<2>(Data::Field::Exponential, true);
        ASSERT_TRUE(baseline.displacement.isFinite());
        ASSERT_TRUE(wrong.displacement.isFinite());
        ASSERT_TRUE(std::isfinite(baseline.strain));
        ASSERT_TRUE(std::isfinite(baseline.stress));
        ASSERT_TRUE(std::isfinite(wrong.strain));
        ASSERT_TRUE(std::isfinite(wrong.stress));
        ASSERT_GT(baseline.displacement.getL2(), 0);
        ASSERT_GT(baseline.strain, 0);
        ASSERT_GT(baseline.stress, 0);
        EXPECT_GT(wrong.displacement.getL2(),
          NegativeControlFactor * baseline.displacement.getL2());
        EXPECT_GT(wrong.strain, NegativeControlFactor * baseline.strain);
        EXPECT_GT(wrong.stress, NegativeControlFactor * baseline.stress);
      }

      template <size_t K>
      void checkSensitivity() const
      {
        Workload<ContextType> problem(this->GetParam(), 5);
        const auto baseline = problem.template solve<K>(Data::Field::Exponential);
        const auto quadrature =
          problem.template solve<K>(Data::Field::Exponential, false, RefinedOrder);
        const auto solver = problem.template solve<K>(
          Data::Field::Exponential, false, AssemblyOrder, RefinedTolerance);
        for (const auto& refined : {quadrature, solver})
          for (const auto pair :
            {std::pair{baseline.displacement.getL2(), refined.displacement.getL2()},
              std::pair{baseline.displacement.getH1Seminorm(),
                refined.displacement.getH1Seminorm()},
              std::pair{baseline.strain, refined.strain},
              std::pair{baseline.stress, refined.stress}})
          {
            ASSERT_TRUE(std::isfinite(pair.first));
            ASSERT_TRUE(std::isfinite(pair.second));
            ASSERT_GT(pair.first, 0);
            EXPECT_LT(std::abs(pair.second / pair.first - 1), SensitivityTolerance);
          }
      }
  };

  using LocalTest = CurvedLinearElasticityTest<Context::Local>;
  TEST_P(LocalTest, LiftedAsymmetricAffineMetricOracle)
  {
    checkLiftedMetric();
  }
  TEST_P(LocalTest, LiftedComplexVectorMetricOracle)
  {
    checkComplexVectorMetric();
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
  TEST_P(LocalTest, RejectsOmittedVolumetricTerm)
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

#if defined(RODIN_CURVED_ELASTICITY_PETSC) && defined(RODIN_USE_MPI)
  using MPITest = CurvedLinearElasticityTest<Context::MPI>;
  TEST_P(MPITest, LiftedAsymmetricAffineMetricOracle)
  {
    checkLiftedMetric();
  }
  TEST_P(MPITest, LiftedComplexVectorMetricOracle)
  {
    checkComplexVectorMetric();
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
  TEST_P(MPITest, RejectsOmittedVolumetricTerm)
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
    H1<2, Math::SpatialVector<Real>, Mesh<Context::MPI>> space(
      std::integral_constant<size_t, 2>{}, mesh, dim);
    PETSc::Variational::GridFunction u(space);
    u = VectorFunction(
      dim, [dim](const Point&) { return Math::SpatialVector<Real>::Zero(dim); });
    const VectorFunction exact(dim, [dim](const Point&) {
      Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
      for (size_t i = 0; i < dim; ++i)
        value(i) = 1;
      return value;
    });
    const auto zeroJacobian = [dim](const Point&) {
      Math::SpatialMatrix<Real> value(
        static_cast<std::uint8_t>(dim), static_cast<std::uint8_t>(dim));
      value.setZero();
      return value;
    };
    const auto error =
      ErrorNorm::computeVector(mesh, u, exact, zeroJacobian, AssemblyOrder);
    const Real volume = dim == 1 ? Real(1.1) : Real(1);
    constexpr Real VolumeTolerance = 1e-12; // Analytic volume; quadrature roundoff only.
    EXPECT_NEAR(error.getL2(), std::sqrt(Real(dim) * volume), VolumeTolerance);
    EXPECT_EQ(error.getH1Seminorm(), 0);
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

#ifdef RODIN_CURVED_ELASTICITY_PETSC
int main(int argc, char** argv)
{
  if (PetscInitialize(&argc, &argv, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result;
  {
#ifdef RODIN_USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
    Rodin::Tests::Convergence::Isoparametric::LinearElasticity::environment = &env;
    Rodin::Tests::Convergence::Isoparametric::LinearElasticity::world = &comm;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
#endif
