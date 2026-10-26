/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Separates field coefficient lookup from reference-gradient contraction.
 */
#include <benchmark/benchmark.h>
#include <string>
#include <type_traits>
#include <vector>
#include "convergence/Convergence.h"

namespace Rodin::Tests::Benchmarks
{
  bool fieldContractionFailed = false;

  enum class FieldContractionOperation
  {
    CoefficientLookup,
    ReferenceContraction,
    PreparedIndexContraction
  };

  /**
   * @brief Measures the consumer of already evaluated reference gradients.
   *
   * ## Mathematical foundation
   *
   * For the cell closure @f$I_K(a)@f$, the reference derivative is
   * @f$\widehat\nabla u_h(\widehat x_q)=
   * \sum_a c_{I_K(a)}\widehat\nabla\phi_a(\widehat x_q)@f$.
   * Coefficients are distinct deterministic values, not interpreted as nodal
   * samples. Real and complex fields use @f$c_i=(i+1)\eta@f$, with
   * @f$\eta=1@f$ and @f$\eta=1+i/2@f$, respectively.
   *
   * ## Architecture
   *
   * One cell, its space, quadrature, public basis gradients and output storage
   * are prepared before timing. CoefficientLookup writes one coefficient per
   * sample and basis through the public local-to-global lookup. ReferenceContraction
   * additionally performs the ordered gradient sum. PreparedIndexContraction
   * uses an owned closure copied before timing to isolate lookup from contraction.
   * It is a caller-side measurement control, not a new library implementation.
   * Physical mapping, derivative generation, assembly and solving are absent.
   * These local scopes do not establish distributed lookup costs or borrow rules.
   *
   * ## Numerical qualification
   *
   * Every prepared closure index must equal the public local-to-global index.
   * Each coefficient and contraction entry is compared exactly with an untimed
   * closure-based ordered expansion before and after timing. Coefficients change
   * through @f$1\to2\to1@f$ without changing the space. Deliberately changed
   * output must be rejected. No timing or approximation tolerance is used.
   * This validates the lookup/contraction decomposition, not independent basis
   * accuracy; that remains covered by element and convergence tests.
   *
   * @tparam K H1 degree, or one for native P1.
   * @tparam Scalar Real or complex scalar field.
   * @tparam NativeP1 Whether to select the distinct native P1 implementation.
   */
  template <size_t K, class Scalar, bool NativeP1 = false>
  class FieldContraction
  {
    public:
      /**
       * @brief Measures a coefficient or reference-contraction sweep.
       * @param state Quadrature order and benchmark iteration state.
       * @param geometry Positive-dimensional reference-cell family.
       * @param operation Measured stage and lookup protocol.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        FieldContractionOperation operation)
      {
        using namespace Variational;
        using Space = std::conditional_t<NativeP1, P1<Scalar>, H1<K, Scalar>>;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
        const auto cell = *mesh.getCell();
        const size_t dim = mesh.getDimension();
        const Index cellIndex = cell.getIndex();
        auto space = [&]() {
          if constexpr (NativeP1)
            return Space(mesh);
          else
            return Space(std::integral_constant<size_t, K>{}, mesh);
        }();
        GridFunction field(space);
        const Scalar eta = [] {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 0.5);
          else
            return Scalar(1);
        }();
        const auto setCoefficients = [&](Real scale) {
          for (Index i = 0; i < space.getSize(); ++i)
            field.getData()(i) = scale * Real(i + 1) * eta;
        };
        setCoefficients(1);
        const auto& element = space.getFiniteElement(dim, cellIndex);
        const size_t count = element.getCount();
        const IndexArray indices = space.getDOFs(dim, cellIndex);
        if (static_cast<size_t>(indices.size()) != count)
        {
          fieldContractionFailed = true;
          state.SkipWithError("Cell closure and element basis counts differ");
          return;
        }
        const auto& formula = QF::QuadratureFormula::get(state.range(0), geometry);
        const size_t nqp = formula.getSize();
        using Gradient = Math::SpatialVector<Scalar>;
        std::vector<Gradient> gradients;
        gradients.reserve(nqp * count);
        for (size_t qp = 0; qp < nqp; ++qp)
        {
          for (size_t a = 0; a < count; ++a)
          {
            gradients.emplace_back(
              element.getBasis(a).getGradient()(formula.getPoint(qp)));
          }
        }
        std::vector<Scalar> coefficients(nqp * count),
          expectedCoefficients(coefficients.size());
        std::vector<Gradient> output(nqp, Gradient(dim)), expected(nqp, Gradient(dim));
        const auto prepareExpected = [&]() {
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            expected[qp].setZero();
            for (size_t a = 0; a < count; ++a)
            {
              const Scalar coefficient = field[indices(a)];
              expectedCoefficients[qp * count + a] = coefficient;
              expected[qp] += gradients[qp * count + a] * coefficient;
            }
          }
        };
        const auto sweep = [&]() {
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            if (operation == FieldContractionOperation::CoefficientLookup)
            {
              for (size_t a = 0; a < count; ++a)
              {
                coefficients[qp * count + a] =
                  field[space.getGlobalIndex({dim, cellIndex}, a)];
              }
            }
            else
            {
              output[qp].setZero();
              for (size_t a = 0; a < count; ++a)
              {
                const Index index =
                  operation == FieldContractionOperation::PreparedIndexContraction
                  ? indices(a)
                  : space.getGlobalIndex({dim, cellIndex}, a);
                output[qp] += gradients[qp * count + a] * field[index];
              }
            }
          }
        };
        const auto matches = [&]() {
          if (static_cast<size_t>(indices.size()) != count)
            return false;
          for (size_t a = 0; a < count; ++a)
          {
            if (indices(a) != space.getGlobalIndex({dim, cellIndex}, a))
              return false;
          }
          if (operation == FieldContractionOperation::CoefficientLookup)
            return coefficients == expectedCoefficients;
          for (size_t qp = 0; qp < nqp; ++qp)
          {
            for (size_t j = 0; j < dim; ++j)
            {
              if (output[qp](j) != expected[qp](j))
                return false;
            }
          }
          return true;
        };
        const auto requireMatch = [&]() {
          if (matches())
            return true;
          fieldContractionFailed = true;
          state.SkipWithError(
            "Coefficient lookup or ordered reference contraction changed");
          return false;
        };
        for (const Real scale : {Real(1), Real(2), Real(1)})
        {
          setCoefficients(scale);
          prepareExpected();
          sweep();
          if (!requireMatch())
            return;
        }
        if (operation == FieldContractionOperation::CoefficientLookup)
          coefficients.front() += Scalar(1);
        else
          output.front()(0) += Scalar(1);
        if (matches())
        {
          fieldContractionFailed = true;
          state.SkipWithError("Perturbed contraction output was accepted");
          return;
        }
        sweep();
        for (auto _ : state)
        {
          sweep();
          benchmark::DoNotOptimize(coefficients.data());
          benchmark::DoNotOptimize(output.data());
          benchmark::ClobberMemory();
        }
        if (!requireMatch())
          return;
        Scalar checksum = 0;
        if (operation == FieldContractionOperation::CoefficientLookup)
        {
          for (const auto& value : coefficients)
            checksum += value;
        }
        else
        {
          for (const auto& gradient : output)
          {
            for (size_t j = 0; j < dim; ++j)
              checksum += gradient(j);
          }
        }
        state.counters["cells"] = 1;
        state.counters["dimension"] = dim;
        state.counters["degree"] = K;
        state.counters["global_dofs"] = space.getSize();
        state.counters["scalar_basis"] = count;
        state.counters["quadrature_order"] = state.range(0);
        state.counters["point_evaluations"] = nqp;
        state.counters["coefficient_reads"] = nqp * count;
        state.counters["scalar_contractions"] =
          operation == FieldContractionOperation::CoefficientLookup ? 0
                                                                    : nqp * count * dim;
        if constexpr (std::is_same_v<Scalar, Complex>)
        {
          state.counters["checksum_real"] = checksum.real();
          state.counters["checksum_imag"] = checksum.imag();
        }
        else
        {
          state.counters["checksum_real"] = checksum;
          state.counters["checksum_imag"] = 0;
        }
      }

      /// @brief Registers each lookup/contraction stage on seven geometries and three orders.
      static void registerCases()
      {
        using G = Geometry::Polytope::Type;
        for (const auto geometry : {G::Segment, G::Triangle, G::Quadrilateral,
               G::Tetrahedron, G::Pyramid, G::Hexahedron, G::Wedge})
        {
          for (const auto operation : {FieldContractionOperation::CoefficientLookup,
                 FieldContractionOperation::ReferenceContraction,
                 FieldContractionOperation::PreparedIndexContraction})
          {
            const std::string stage =
              operation == FieldContractionOperation::CoefficientLookup
              ? "CoefficientLookup"
              : operation == FieldContractionOperation::ReferenceContraction
              ? "ReferenceContraction"
              : "PreparedIndexContraction";
            const std::string name = "FieldContraction/" + stage +
              (NativeP1 ? "/P1/" : "/H1/P" + std::to_string(K) + "/") +
              (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
              std::string(Convergence::UniformGrid::getGeometryName(geometry));
            benchmark::RegisterBenchmark(name.c_str(),
              [geometry, operation](auto& state) { run(state, geometry, operation); })
              ->Arg(4)
              ->Arg(8)
              ->Arg(16)
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
  const auto registerScalar = []<class Scalar>() {
    FieldContraction<1, Scalar, true>::registerCases();
    FieldContraction<1, Scalar>::registerCases();
    FieldContraction<2, Scalar>::registerCases();
    FieldContraction<3, Scalar>::registerCases();
  };
  registerScalar.template operator()<Real>();
  registerScalar.template operator()<Complex>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Local field lookup/contraction; prepared reference gradients");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || fieldContractionFailed ? 1 : 0;
}
