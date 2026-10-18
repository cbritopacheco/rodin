/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Public H1 scalar mass assembly with both coefficient operand orders.
 */
#include "MassAssembly.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Selects scalar coefficient placement without changing the public workload.
   *
   * The forms are @f$a_L(u,v)=\int_\Omega(\rho u)\overline v\,dx@f$
   * and @f$a_R(u,v)=\int_\Omega(u\rho)\overline v\,dx@f$. The inherited
   * workload supplies the same scalar coefficient, constant solution,
   * projection load and independent moments. Only the left operator is
   * replaced; the existing right-weighted workload remains unchanged.
   *
   * Shared public-assembly qualification checks the complete original-entry
   * operator, constant-field energy, projection identity, repeated replacement
   * and doubled coefficient. Deliberately doubled operators and loads must
   * be rejected. Geometry and oracle construction remain untimed.
   *
   * @tparam K H1 polynomial degree.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Left Whether the scalar coefficient precedes the trial field.
   */
  template <size_t K, class Scalar, bool Left>
  class MassOperandWorkload : public MassWorkload<K, Scalar, false, MassFamily::H1>
  {
    public:
      using Base = MassWorkload<K, Scalar, false, MassFamily::H1>;

      /**
       * @brief Configures the selected expression and the inherited projection load.
       * @tparam Reference Whether the original entry loop supplies the operator.
       * @param u Trial function and interpolated constant solution.
       * @param v Test function.
       * @param matrix Public operator form.
       * @param load Public projection-load form.
       * @param dim Domain dimension.
       * @param reaction Whether a variable scalar coefficient is requested.
       * @param scale Borrowed coefficient multiplier.
       * @param curved Whether the exact quadratic physical map is installed.
       * @returns Independent constant-field action at unit coefficient scale.
       */
      template <bool Reference = false, class U, class V, class A, class B>
      static Scalar configure(U& u, V& v, A& matrix, B& load, size_t dim, bool reaction,
        const Real& scale, bool curved)
      {
        using namespace Variational;
        const Scalar expected = Base::template configure<Reference>(
          u, v, matrix, load, dim, reaction, scale, curved);
        if constexpr (Left)
        {
          if (reaction)
          {
            const auto evaluate = [dim, &scale](const Geometry::Point& point) {
              Real rho = 1;
              for (size_t j = 0; j < dim; ++j)
                rho += point(j);
              if constexpr (std::is_same_v<Scalar, Complex>)
                return scale * Scalar(1, 0.5) * rho;
              else
                return scale * rho;
            };
            const auto coefficient = [&]() {
              if constexpr (std::is_same_v<Scalar, Complex>)
                return ComplexFunction(evaluate);
              else
                return RealFunction(evaluate);
            }();
            auto integral = Integral(coefficient * u, v);
            integral.setOrder(8);
            if constexpr (Reference)
              matrix = ReferenceIntegral(integral);
            else
              matrix = integral;
          }
        }
        return expected;
      }

      /// @brief Registers operator/load measurements at three sizes on each family.
      static void registerCases()
      {
        using Geometry::Polytope;
        using Assembly =
          MassAssembly<K, Scalar, false, MassFamily::H1, MassOperandWorkload>;
        for (auto geometry :
          {Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (bool curved : {false, true})
          {
            for (bool rhs : {false, true})
            {
              const std::string name = std::string("Assembly/MassOperandOrder/") +
                (Left ? "LeftScalarWeighted/" : "RightScalarWeighted/") + "H1/P" +
                std::to_string(K) + "/" +
                (std::is_same_v<Scalar, Complex> ? "ComplexScalar/" : "RealScalar/") +
                std::string(Convergence::UniformGrid::getGeometryName(geometry)) +
                (curved ? "/CurvedQ2/" : "/Affine/") + (rhs ? "Load" : "Operator");
              benchmark::RegisterBenchmark(name.c_str(),
                [geometry, curved, rhs](
                  auto& state) { Assembly::run(state, geometry, true, rhs, curved); })
                ->Arg(3)
                ->Arg(5)
                ->Arg(9)
                ->UseRealTime();
            }
          }
        }
      }
  };

  template <size_t K>
  void registerMassOperandOrder()
  {
    MassOperandWorkload<K, Real, true>::registerCases();
    MassOperandWorkload<K, Real, false>::registerCases();
    MassOperandWorkload<K, Complex, true>::registerCases();
    MassOperandWorkload<K, Complex, false>::registerCases();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerMassOperandOrder<1>();
  registerMassOperandOrder<2>();
  registerMassOperandOrder<3>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential");
#endif
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "scope", "Warmed scalar-weighted operator/load assembly; no solve");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !cases || massCheckFailed ? 1 : 0;
}
