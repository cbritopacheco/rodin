/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <type_traits>
#include "Rodin/Assembly.h"
#include "Rodin/Assembly/Sequential.h"
#include "QuadratureReference.h"
#include "convergence/Convergence.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool convectionAssemblyFailed = false;

  /**
   * @brief Measures assembly consumers of the H1 directional-gradient kernel.
   *
   * A fixed physical form is measured through complete assembly, triplet
   * generation and sparse completion, with independent untimed qualification.
   *
   * ## Mathematical foundation
   *
   * The operator is @f$a_b(u,v)=\int_\Omega(Ju)b\cdot\overline v\,dx@f$.
   * The affine field @f$u_c=\alpha_c\sum_jx_j@f$ has constant directional
   * derivative @f$f_c=\alpha_c\sum_jb_j@f$. Public interpolation and load
   * assembly qualify @f$Ax_\star=\ell_f@f$, independently of the original
   * expression-entry operator. Doubling @f$x_\star@f$, rather than adding a
   * constant in the derivative nullspace, must violate this identity.
   *
   * ## Architecture
   *
   * Complete measures public warmed bilinear-form assembly. Triplets measures
   * the existing sequential assembly backend without sparse completion.
   * SparseCompletion measures replacement from a fixed qualified triplet list.
   * These are separate operations, not timings to subtract for attribution.
   * CPU time covers the process, including OpenMP workers; elapsed time is
   * reported independently and remains the primary comparison measure.
   * Mesh/space construction, interpolation, reference assembly and all guards
   * remain untimed. There are no constraints or solves; an unconstrained
   * transport operator is not assigned a fictitious coercive PDE solve.
   * Three components are independent of the physical dimension. Seven native
   * positive-dimensional cell families have three mesh sizes; Point's zero
   * derivative remains covered by the Variational-owned component executable.
   * Benchmarks are opt-in and are not CTest or publication timing gates.
   *
   * @tparam Scalar Field scalar type.
   * @tparam Degree H1 polynomial degree.
   * @tparam ComplexDirection Whether the direction has a nonzero imaginary part.
   */
  template <class Scalar, size_t Degree, bool ComplexDirection>
    requires(!ComplexDirection || std::is_same_v<Scalar, Complex>)
  class ConvectionAssembly
  {
    public:
      /// @brief Distinct public assembly operations with independent timed scopes.
      enum class Phase
      {
        Complete,
        Triplets,
        SparseCompletion
      };

      /// @brief Component count held fixed across dimensions and geometries.
      static constexpr size_t Components = 3;
      /// @brief Same base integration order as the directional component controls.
      static constexpr size_t QuadratureOrder = 4;
      /// @brief Relative budget for differently associated expression-entry assembly.
      static constexpr Real ReferenceBudget = 1e-11;
      /// @brief Relative budget for the independently derived affine patch identity.
      static constexpr Real PatchBudget = 1e-10;
      /// @brief Rejection floor for the doubled affine field, relative to load scale.
      static constexpr Real SensitivityFloor = 1e-4;

      /**
       * @brief Measures one warmed assembly operation after complete qualification.
       * @param state Grid points per axis and benchmark iteration state.
       * @param geometry Native unit-box cell family.
       * @param phase Operation included in the timed region.
       */
      static void run(
        benchmark::State& state, Geometry::Polytope::Type geometry, Phase phase)
      {
        using namespace Variational;
        using Space = H1<Degree, Math::SpatialVector<Scalar>>;
        using DirectionScalar = std::conditional_t<ComplexDirection, Complex, Real>;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dimension = mesh.getDimension();
        Space space(std::integral_constant<size_t, Degree>{}, mesh, Components);
        TrialFunction u(space);
        TestFunction v(space);
        Math::SpatialVector<DirectionScalar> beta(dimension);
        for (size_t j = 0; j < dimension; ++j)
        {
          if constexpr (ComplexDirection)
            beta[j] = Complex(Real(j + 1), Real(j + 2));
          else
            beta[j] = Real(j + 1);
        }
        const VectorFunction direction(
          dimension, [beta](const Geometry::Point&) { return beta; });
        const auto amplitude = [](size_t component) -> Scalar {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Complex(Real(component + 1), Real(component + 2));
          else
            return Real(component + 1);
        };
        u.getSolution() = VectorFunction(
          Components, [dimension, amplitude](const Geometry::Point& point) {
            Real coordinateSum = 0;
            for (size_t j = 0; j < dimension; ++j)
              coordinateSum += point(j);
            Math::SpatialVector<Scalar> value(Components);
            for (size_t c = 0; c < Components; ++c)
              value[c] = amplitude(c) * coordinateSum;
            return value;
          });
        Scalar betaSum = 0;
        for (size_t j = 0; j < dimension; ++j)
          betaSum += beta[j];
        const VectorFunction source(
          Components, [betaSum, amplitude](const Geometry::Point&) {
            Math::SpatialVector<Scalar> value(Components);
            for (size_t c = 0; c < Components; ++c)
              value[c] = amplitude(c) * betaSum;
            return value;
          });
        const auto integrand = Dot(Jacobian(u) * direction, v);
        Integral integral(integrand);
        static_assert(decltype(integral)::Specialized);
        integral.setOrder(QuadratureOrder);
        auto loadIntegral = Integral(source, v);
        loadIntegral.setOrder(QuadratureOrder);
        BilinearForm form(u, v), reference(u, v);
        LinearForm load(v);
        form = integral;
        reference = ReferenceIntegral(integral, integrand);
        load = loadIntegral;
        reference.assemble();
        load.assemble();
        form.assemble();
        const auto expected = reference.getOperator();
        const Math::Vector<Scalar> exact = u.getSolution().getData();
        const Math::Vector<Scalar> wrong = Scalar(2) * exact;
        const Real matrixScale = std::max(Real(1), expected.norm());
        const Real loadScale = std::max(Real(1), load.getVector().norm());
        const auto valid = [&](const auto& matrix) {
          const Real matrixError = (matrix - expected).norm() / matrixScale;
          const Real patchError = (matrix * exact - load.getVector()).norm() / loadScale;
          const Real wrongError = (matrix * wrong - load.getVector()).norm() / loadScale;
          return std::isfinite(matrixError) && matrixError <= ReferenceBudget &&
            std::isfinite(patchError) && patchError <= PatchBudget &&
            std::isfinite(wrongError) && wrongError > SensitivityFloor;
        };
        auto perturbed = expected;
        perturbed.coeffRef(0, 0) += Scalar(matrixScale);
        if (!valid(form.getOperator()) || valid(perturbed))
        {
          convectionAssemblyFailed = true;
          state.SkipWithError(
            "Directional operator or affine patch qualification failed");
          return;
        }

        using Triplets = std::vector<Math::SparseTriplet<Scalar>>;
        using Form = decltype(form);
        using TripletForm =
          BilinearForm<typename Form::SolutionType, Space, Space, Triplets>;
        Assembly::Sequential<Triplets, TripletForm> tripletAssembly;
        Triplets triplets;
        const auto input = [&]() {
          return typename decltype(tripletAssembly)::InputType{
            space, space, form.getLocalIntegrators(), form.getGlobalIntegrators()};
        };
        Math::SparseMatrix<Scalar> completed(space.getSize(), space.getSize());
        tripletAssembly.execute(triplets, input());
        completed.setFromTriplets(triplets.begin(), triplets.end());
        if (!valid(completed))
        {
          convectionAssemblyFailed = true;
          state.SkipWithError("Directional triplets differ from the qualified operator");
          return;
        }
        for (auto _ : state)
        {
          if (phase == Phase::Complete)
          {
            form.assemble();
            benchmark::DoNotOptimize(form.getOperator().valuePtr());
          }
          else if (phase == Phase::Triplets)
          {
            tripletAssembly.execute(triplets, input());
            benchmark::DoNotOptimize(triplets.data());
          }
          else
          {
            completed.setFromTriplets(triplets.begin(), triplets.end());
            benchmark::DoNotOptimize(completed.valuePtr());
          }
          benchmark::ClobberMemory();
        }
        if (phase == Phase::Triplets)
          completed.setFromTriplets(triplets.begin(), triplets.end());
        const auto& result = phase == Phase::Complete ? form.getOperator() : completed;
        if (!valid(result))
        {
          convectionAssemblyFailed = true;
          state.SkipWithError("Repeated directional assembly changed the operator");
        }
        size_t localEntries = 0;
        const size_t cells = mesh.getPolytopeCount(dimension);
        for (size_t cell = 0; cell < cells; ++cell)
        {
          const size_t count = space.getDOFs(dimension, cell).size();
          localEntries += count * count;
        }
        const size_t samples =
          QF::PolytopeQuadratureFormula::get(QuadratureOrder, geometry).getSize();
        state.counters["dimension"] = dimension;
        state.counters["degree"] = Degree;
        state.counters["components"] = Components;
        state.counters["field_complex"] = std::is_same_v<Scalar, Complex>;
        state.counters["direction_complex"] = ComplexDirection;
        state.counters["cells"] = cells;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = result.nonZeros();
        state.counters["triplets"] = triplets.size();
        state.counters["local_entries"] = localEntries;
        state.counters["quadrature_order"] = QuadratureOrder;
        state.counters["samples_per_cell"] = samples;
        state.counters["potential_component_updates"] =
          localEntries / Components * samples;
        state.counters["kernel_work_timed"] = phase != Phase::SparseCompletion;
        state.counters["sparse_completion_timed"] = phase != Phase::Triplets;
        state.counters["relative_matrix_error"] =
          (result - expected).norm() / matrixScale;
        state.counters["relative_patch_error"] =
          (result * exact - load.getVector()).norm() / loadScale;
      }

      /// @brief Registers matching scalar fields, geometries, phases and three sizes.
      static void registerCases()
      {
        using Geometry::Polytope;
        for (const auto geometry :
          {Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (const auto phase :
            {Phase::Complete, Phase::Triplets, Phase::SparseCompletion})
          {
            const char* phaseName = phase == Phase::Complete ? "Complete"
              : phase == Phase::Triplets                     ? "Triplets"
                                                             : "SparseCompletion";
            const std::string name = "Convection/Assembly/H1/P" + std::to_string(Degree) +
              (std::is_same_v<Scalar, Real> ? "/RealField/" : "/ComplexField/") +
              (ComplexDirection ? "ComplexDirection/" : "RealDirection/") +
              std::string(Convergence::UniformGrid::getGeometryName(geometry)) + "/" +
              phaseName;
            benchmark::RegisterBenchmark(name.c_str(),
              [geometry, phase](auto& state) { run(state, geometry, phase); })
              ->Arg(2)
              ->Arg(3)
              ->Arg(4)
              ->MeasureProcessCPUTime()
              ->UseRealTime();
          }
        }
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  const auto degree = []<size_t Degree>() {
    ConvectionAssembly<Real, Degree, false>::registerCases();
    ConvectionAssembly<Complex, Degree, false>::registerCases();
    ConvectionAssembly<Complex, Degree, true>::registerCases();
  };
  degree.template operator()<1>();
  degree.template operator()<2>();
  degree.template operator()<3>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("complete_assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("complete_assembly_backend", "Eigen/Sequential");
#endif
  benchmark::AddCustomContext("triplet_backend", "Eigen/Sequential");
  benchmark::AddCustomContext(
    "scope", "Directional assembly phases; no constraints or solve");
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || convectionAssemblyFailed ? 1 : 0;
}
