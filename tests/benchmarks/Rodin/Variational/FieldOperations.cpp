/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Interpolation, point evaluation and L2 projection feature measurements.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"
#include "convergence/Convergence.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool fieldOperationFailed = false;

  /**
   * @brief Measures distinct field operations with a representable affine oracle.
   *
   * For @f$g(x)=1+\sum_{j=0}^{d-1}(j+1)x_j@f$, scalar fields are
   * @f$f=\alpha g@f$ and vector fields have components
   * @f$f_j=(j+1)\alpha g@f$. Here @f$\alpha=1@f$ or @f$1+i/2@f$.
   * On the affine unit box, every selected space represents these fields.
   * Thus interpolation and the orthogonal projection recover the same field,
   * although only projection assembles and solves @f$Mc=b@f$.
   *
   * ## Architecture
   *
   * Meshes, spaces, analytic fields and warmup are untimed. Interpolation
   * alternates an additive field offset to require actual coefficient updates.
   * Projection assembly measures the public problem assembly, solve measures
   * CG setup and iterations, and AssemblyAndSolve includes both. Each solve
   * starts from zero, reset outside timing. Assembly-only cases solve again
   * outside timing: reassembly imports the trial field as the initial guess,
   * not the previous direct linear-system solution. No constraints are used.
   * Numerical guards precede and follow timing: analytic point values on every
   * cell, coefficient error, residual, analytic mass energy and sensitivity to
   * a wrong coefficient vector. Projection checks do not assume a DOF layout.
   * This executable is opt-in and is not a CTest or CI timing gate.
   * PointEvaluation and IntegrationPointEvaluation sweep preconstructed mapped
   * samples through the public GridFunction interface. Physical maps and
   * analytic values are warmed before timing. The integration-point variant
   * carries formula provenance but changes point/cell throughout the sweep;
   * it is not a repeated-point cache-hit measurement. Guards also change the
   * coefficients without changing the space or sample identities.
   *
   * @tparam K Degree of H1, or one for the distinct P1 implementation.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Vector Whether the field has one component per spatial dimension.
   * @tparam Vertex Whether to use P1 rather than H1.
   */
  template <size_t K, class Scalar, bool Vector, bool Vertex = false>
  class FieldOperations
  {
    public:
    /// @brief Dimensionless point-value budget for the representable field.
    static constexpr Real FieldTolerance = 1e-10;
    /// @brief Relative mass-energy budget for quadrature and roundoff.
    static constexpr Real EnergyTolerance = 1e-10;
    /// @brief Algebraic residual budget, looser than the CG stopping tolerance.
    static constexpr Real ResidualTolerance = 1e-11;
    /// @brief Relative coefficient budget allowing mass-matrix conditioning.
    static constexpr Real CoefficientTolerance = 1e-9;
    /// @brief CG stopping tolerance, matching native convergence solves.
    static constexpr Real SolverTolerance = 1e-13;
    /// @brief Iteration cap matching native convergence solves.
    static constexpr size_t MaximumIterations = 20000;
    /// @brief Wrong-field residual floor ensuring a non-vacuous oracle.
    static constexpr Real SensitivityFloor = 1e-4;
    /// @brief Fixed quadrature order defining the measured mass operator.
    static constexpr size_t QuadratureOrder = 8;

    /**
     * @brief Measures one field feature on an affine mesh.
     * @param state Grid points per axis and benchmark iterations.
     * @param geometry Reference-cell family.
     * @param operation Interpolation, projection or point-evaluation stage index.
     */
    static void run(
      benchmark::State& state, Geometry::Polytope::Type geometry, size_t operation)
    {
      using namespace Variational;
      using Range = std::conditional_t<Vector, Math::SpatialVector<Scalar>, Scalar>;
      using Space = std::conditional_t<Vertex, P1<Range>, H1<K, Range>>;
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
      const size_t dim = mesh.getDimension();
      auto space = [&]() {
        if constexpr (Vertex && Vector)
          return Space(mesh, dim);
        else if constexpr (Vertex)
          return Space(mesh);
        else if constexpr (Vector)
          return Space(std::integral_constant<size_t, K>{}, mesh, dim);
        else
          return Space(std::integral_constant<size_t, K>{}, mesh);
      }();
      const Scalar phase = []() {
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Scalar(1, Real(1) / 2);
        else
          return Scalar(1);
      }();
      Real offset = 0;
      const auto scalar = [dim, phase, &offset](const Geometry::Point& p) {
        Real value = 1 + offset;
        for (size_t j = 0; j < dim; ++j)
          value += Real(j + 1) * p(j);
        return phase * value;
      };
      const auto field = [&]() {
        if constexpr (Vector)
          return VectorFunction(dim, [dim, scalar](const Geometry::Point& p) {
            Range value(dim);
            for (size_t j = 0; j < dim; ++j)
              value(j) = Real(j + 1) * scalar(p);
            return value;
          });
        else if constexpr (std::is_same_v<Scalar, Complex>)
          return ComplexFunction(scalar);
        else
          return RealFunction(scalar);
      }();
      const auto& formula = QF::PolytopeQuadratureFormula::get(QuadratureOrder, geometry);
      GridFunction expected(space);
      expected = field;
      const auto fieldError = [&](const auto& measured) {
        Real error = 0;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& points = cell->getQuadrature(formula);
          for (size_t qp = 0; qp < points.getSize(); ++qp)
          {
            const auto& point = points.getPoint(qp);
            const auto actual = measured.getValue(point);
            const auto exact = field(point);
            const Real difference = [&]() {
              if constexpr (Vector)
                return (actual - exact).norm();
              else
                return std::abs(actual - exact);
            }();
            if (!std::isfinite(difference))
              return std::numeric_limits<Real>::infinity();
            error = std::max(error, difference);
          }
        }
        return error;
      };
      if (fieldError(expected) > FieldTolerance)
      {
        fieldOperationFailed = true;
        state.SkipWithError("Interpolated oracle disagrees with analytic field");
        return;
      }
      state.counters["dimension"] = dim;
      state.counters["cells"] = mesh.getCellCount();
      state.counters["dofs"] = space.getSize();
      state.counters["local_dofs_per_cell"] = space.getDOFs(dim, 0).size();
      state.counters["degree"] = K;
      state.counters["components"] = Vector ? dim : 1;
      state.counters["oracle_quadrature_order"] = QuadratureOrder;
      state.counters["oracle_points_per_cell"] = formula.getSize();
      if (operation == 4 || operation == 5)
      {
        std::vector<IntegrationPoint> points;
        std::vector<Range> exact, measured;
        for (auto cell = mesh.getCell(); cell; ++cell)
        {
          const auto& quadrature = cell->getQuadrature(formula);
          for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
          {
            points.emplace_back(quadrature.getPoint(qp), &formula, qp);
            exact.emplace_back(field(quadrature.getPoint(qp)));
          }
        }
        measured.resize(points.size());
        const auto evaluate = [&]() {
          for (size_t i = 0; i < points.size(); ++i)
          {
            measured[i] = operation == 4 ? expected.getValue(points[i].getPoint())
                                         : expected.getValue(points[i]);
          }
        };
        const auto error = [&]() {
          Real maximum = 0;
          for (size_t i = 0; i < points.size(); ++i)
          {
            const Real difference = [&]() {
              if constexpr (Vector)
                return (measured[i] - exact[i]).norm();
              else
                return std::abs(measured[i] - exact[i]);
            }();
            if (!std::isfinite(difference))
              return std::numeric_limits<Real>::infinity();
            maximum = std::max(maximum, difference);
          }
          return maximum;
        };
        evaluate();
        const Real before = error();
        if constexpr (Vector)
          measured.front()(0) += Scalar(1);
        else
          measured.front() += Scalar(1);
        const Real sensitivity = error();
        offset = 1;
        expected = field;
        evaluate();
        for (size_t i = 0; i < points.size(); ++i)
          exact[i] = field(points[i].getPoint());
        const Real changed = error();
        offset = 0;
        expected = field;
        for (size_t i = 0; i < points.size(); ++i)
          exact[i] = field(points[i].getPoint());
        evaluate();
        if (before > FieldTolerance || changed > FieldTolerance ||
          error() > FieldTolerance || !(sensitivity > SensitivityFloor))
        {
          fieldOperationFailed = true;
          state.SkipWithError("Point-evaluation oracle or coefficient refresh failed");
          return;
        }
        for (auto _ : state)
        {
          evaluate();
          benchmark::DoNotOptimize(measured.data());
          benchmark::ClobberMemory();
        }
        const Real after = error();
        state.counters["absolute_field_error"] = after;
        state.counters["point_evaluations"] = points.size();
        state.counters["scalar_basis_evaluations_per_point"] =
          space.getDOFs(dim, 0).size() / (operation == 4 && Vector ? dim : 1);
        state.counters["scalar_output_entries"] = points.size() * (Vector ? dim : 1);
        if (after > FieldTolerance)
        {
          fieldOperationFailed = true;
          state.SkipWithError("Timed point evaluation changed the represented field");
        }
        return;
      }
      if (operation == 0)
      {
        GridFunction measured(space);
        measured = field;
        const Real before = fieldError(measured);
        measured.getData()(0) += Scalar(1);
        const Real wrongFieldError = fieldError(measured);
        measured = field;
        for (auto _ : state)
        {
          offset = 1 - offset;
          measured = field;
          benchmark::DoNotOptimize(measured.getData().data());
          benchmark::ClobberMemory();
        }
        const Real error = fieldError(measured);
        state.counters["absolute_field_error"] = error;
        if (before > FieldTolerance || error > FieldTolerance ||
          !(wrongFieldError > SensitivityFloor))
        {
          fieldOperationFailed = true;
          state.SkipWithError("Changing-field interpolation failed analytic oracle");
        }
        return;
      }
      TrialFunction u(space);
      state.counters["quadrature_order"] = QuadratureOrder;
      state.counters["quadrature_points_per_cell"] = formula.getSize();
      TestFunction v(space);
      auto mass = Integral(u, v);
      auto load = Integral(field, v);
      mass.setOrder(QuadratureOrder);
      load.setOrder(QuadratureOrder);
      Problem projection(u, v);
      projection = mass - load;
      projection.assemble();
      auto& system = projection.getLinearSystem();
      const auto baselineOperator = system.getOperator();
      const auto baselineLoad = system.getVector();
      const auto coefficients = expected.getData().eval();
      const Real rhsScale = std::max(Real(1), system.getVector().norm());
      const Real coefficientScale = std::max(Real(1), coefficients.norm());
      const Real sum = Real(dim * (dim + 1)) / 2;
      const Real squareSum = Real(dim * (dim + 1) * (2 * dim + 1)) / 6;
      const Real energy = std::norm(phase) * (Vector ? squareSum : 1) *
        (1 + sum + sum * sum / 4 + squareSum / 12);
      Solver::CG solver(projection);
      solver.setTolerance(SolverTolerance).setMaxIterations(MaximumIterations);
      system.getSolution().setZero();
      solver.solve(system);
      const size_t iterations = solver.getIterationNumber();
      const auto residual = [&]() {
        return (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          rhsScale;
      };
      const auto coefficientError = [&]() {
        return (system.getSolution() - coefficients).norm() / coefficientScale;
      };
      const auto valid = [&]() {
        const Scalar action = coefficients.dot(system.getOperator() * coefficients);
        return solver.success() && solver.getIterationNumber() == iterations &&
          std::isfinite(residual()) && residual() < ResidualTolerance &&
          std::isfinite(coefficientError()) &&
          coefficientError() < CoefficientTolerance &&
          std::abs(action / energy - Scalar(1)) < EnergyTolerance &&
          (system.getOperator() - baselineOperator).norm() <
          ResidualTolerance * std::max(Real(1), baselineOperator.norm()) &&
          (system.getVector() - baselineLoad).norm() < ResidualTolerance * rhsScale &&
          (system.getOperator() * coefficients - system.getVector()).norm() / rhsScale <
          ResidualTolerance;
      };
      const Math::Vector<Scalar> wrong = coefficients.array() + Scalar(1);
      const Real sensitivity =
        (system.getOperator() * wrong - system.getVector()).norm() / rhsScale;
      if (!valid() || !(sensitivity > SensitivityFloor))
      {
        fieldOperationFailed = true;
        state.SkipWithError("Projection warmup failed mass, energy or solution oracle");
        return;
      }
      bool solvesValid = true;
      for (auto _ : state)
      {
        if (operation != 1)
        {
          state.PauseTiming();
          system.getSolution().setZero();
          state.ResumeTiming();
        }
        if (operation != 2)
          projection.assemble();
        if (operation != 1)
          solver.solve(system);
        benchmark::DoNotOptimize(system.getSolution().data());
        benchmark::ClobberMemory();
        if (operation != 1)
        {
          state.PauseTiming();
          solvesValid =
            solvesValid && solver.success() && solver.getIterationNumber() == iterations;
          state.ResumeTiming();
        }
      }
      if (operation == 1)
      {
        system.getSolution().setZero();
        solver.solve(system);
      }
      u.getSolution().getData() = system.getSolution();
      const Real error = fieldError(u.getSolution());
      if (!solvesValid || !valid() || error > FieldTolerance)
      {
        fieldOperationFailed = true;
        state.SkipWithError("Timed projection failed numerical or repeatability oracle");
      }
      state.counters["nnz"] = system.getOperator().nonZeros();
      state.counters["matrix_norm"] = system.getOperator().norm();
      state.counters["rhs_norm"] = system.getVector().norm();
      state.counters["solver_reported_iterations"] = iterations;
      state.counters["relative_residual"] = residual();
      state.counters["relative_coefficient_error"] = coefficientError();
      state.counters["absolute_field_error"] = error;
    }

    /// @brief Registers six operations at three mesh sizes on seven geometries.
    static void registerCases()
    {
      using G = Geometry::Polytope::Type;
      const char* operations[] = {"Interpolation", "ProjectionAssembly",
        "ProjectionSolve", "AssemblyAndSolve", "PointEvaluation",
        "IntegrationPointEvaluation"};
      for (size_t operation = 0; operation < std::size(operations); ++operation)
      {
        for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
               G::Pyramid, G::Hexahedron, G::Wedge})
        {
          const std::string name = std::string("Field/") + operations[operation] + "/" +
            (Vertex ? "P1" : "H1/P" + std::to_string(K)) + "/" +
            (std::is_same_v<Scalar, Complex> ? "Complex" : "Real") + "/" +
            (Vector ? "Vector/" : "Scalar/") +
            std::string(Convergence::UniformGrid::getGeometryName(geometry));
          benchmark::RegisterBenchmark(name.c_str(),
            [geometry, operation](auto& state) { run(state, geometry, operation); })
            ->Arg(3)
            ->Arg(5)
            ->Arg(9)
            ->UseRealTime();
        }
      }
    }
  };

  /**
   * @brief Registers the distinct P1 and degree-one/two H1 implementations.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Vector Whether to measure a vector-valued field.
   */
  template <class Scalar, bool Vector>
  void registerFieldOperations()
  {
    FieldOperations<1, Scalar, Vector, true>::registerCases();
    FieldOperations<1, Scalar, Vector>::registerCases();
    FieldOperations<2, Scalar, Vector>::registerCases();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  registerFieldOperations<Real, false>();
  registerFieldOperations<Real, true>();
  registerFieldOperations<Complex, false>();
  registerFieldOperations<Complex, true>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential");
#endif
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || fieldOperationFailed ? 1 : 0;
}
