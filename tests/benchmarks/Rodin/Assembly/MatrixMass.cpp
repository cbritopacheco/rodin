/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Public matrix-field reaction operator and projection-load assembly.
 */
#include "MassAssembly.h"
#include "Workloads/MatrixMassWorkload.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Registers actual matrix-basis assembly consumers at three mesh sizes.
   *
   * The rectangular field has six components, independent of domain dimension.
   * Operator/load scopes share the existing public assembly qualification.
   * Setup and numerical guards are untimed. Sizes three, four and five expose
   * scaling without silently launching degree-three matrix problems at the
   * scalar workload's largest registered grid.
   *
   * @tparam K H1 degree.
   * @tparam Scalar Real or complex field scalar.
   * @tparam Left Scalar coefficient operand order.
   */
  template <size_t K, class Scalar, bool Left>
  void registerMatrixMass()
  {
    using Geometry::Polytope;
    using Workload = MatrixMassWorkload<K, Scalar, 2, 3, Left>;
    using Assembly = MassAssembly<K, Scalar, false, MassFamily::H1, Workload>;
    for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      for (bool curved : {false, true})
      {
        for (bool rhs : {false, true})
        {
          const std::string name = std::string("Assembly/MatrixMass/") +
            (Left ? "LeftScalarWeighted/" : "RightScalarWeighted/") + "H1/P" +
            std::to_string(K) + "/" +
            (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
            std::string(Convergence::UniformGrid::getGeometryName(geometry)) +
            (curved ? "/CurvedQ2/" : "/Affine/") + (rhs ? "Load" : "Operator");
          benchmark::RegisterBenchmark(name.c_str(),
            [geometry, curved, rhs](auto& state) {
              Assembly::run(state, geometry, true, rhs, curved);
              state.counters["field_rows"] = 2;
              state.counters["field_columns"] = 3;
            })
            ->Arg(3)
            ->Arg(4)
            ->Arg(5)
            ->UseRealTime();
        }
      }
    }
  }

  template <size_t K>
  void registerMatrixMassDegree()
  {
    registerMatrixMass<K, Real, true>();
    registerMatrixMass<K, Real, false>();
    registerMatrixMass<K, Complex, true>();
    registerMatrixMass<K, Complex, false>();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerMatrixMassDegree<1>();
  registerMatrixMassDegree<2>();
  registerMatrixMassDegree<3>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return massCheckFailed ? 1 : 0;
}
