/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Matched SpatialVector arithmetic across spatial and Eigen representations.
 */
#include "SpatialArithmetic.h"

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    bool failed = false;

    /**
     * @brief Registers identical mathematical operations for one representation.
     * @tparam Scalar Real or complex scalar type.
     * @tparam Storage Runtime, bounded, dynamic or compile-time-sized storage.
     * @tparam Rows Compile-time row extent for Fixed comparisons.
     * @tparam Cols Compile-time column extent for Fixed comparisons.
     * @param label Representation name in benchmark output.
     */
    template <class Scalar, SpatialStorage Storage, int Rows = Eigen::Dynamic,
      int Cols = Eigen::Dynamic>
    void registerCases(const char* label)
    {
      using Fixture = SpatialArithmetic<Scalar, false, Storage, Rows, Cols>;
      using Op = SpatialOperation;
      const std::string prefix = "SpatialVector/" +
        std::string(std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/");
      const auto add = [&]<Op Operation, bool Destination = false>(
                         const char* operation) {
        const std::string name =
          prefix + operation + "/" + label + (Destination ? "Destination" : "");
        auto* item = benchmark::RegisterBenchmark(name.c_str(), [](auto& state) {
          Fixture::template run<Operation, Destination>(state);
          failed = failed || Fixture::failed;
        });
        for (int rows = 0; rows <= 3; ++rows)
        {
          for (int cols = 1; cols <= 1; ++cols)
          {
            if constexpr (Storage == SpatialStorage::Fixed)
            {
              if (rows != Rows || cols != 1)
                continue;
            }
            if constexpr (Operation == Op::Cross)
            {
              if (rows != 3)
                continue;
            }
            if constexpr (Operation == Op::Normalized || Operation == Op::Normalize)
            {
              if (!rows)
                continue;
            }
            for (int count : {1, 32})
              item->Args({rows, count});
          }
        }
        item->ArgNames({"rows", "count"});
      };
      add.template operator()<Op::Zero>("Zero");
      add.template operator()<Op::Copy>("Copy");
      add.template operator()<Op::Add>("Add");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Subtract>("Subtract");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Negate>("Negate");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::ScaleLeft>("ScaleLeft");
      add.template operator()<Op::ScaleRight>("ScaleRight");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Divide>("Divide");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::AddAssign>("AddAssign");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::ScaleAssign>("ScaleAssign");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathSum>("MathSum");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathSubtract>("MathSubtract");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathNegate>("MathNegate");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathScale>("MathScale");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathDivide>("MathDivide");
      add.template operator()<Op::ScalarDot>("ScalarDot");
      add.template operator()<Op::MathDot>("MathDot");
      add.template operator()<Op::RightDot>("RightDot");
      add.template operator()<Op::SquaredNorm>("SquaredNorm");
      add.template operator()<Op::Norm>("Norm");
      add.template operator()<Op::Conjugate>("Conjugate");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathConjugate>("MathConjugate");
      add.template operator()<Op::Transpose>("Transpose");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::SubtractAssign>("SubtractAssign");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::DivideAssign>("DivideAssign");
      add.template operator()<Op::BilinearDot>("BilinearDot");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::StableNorm>("StableNorm");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::BlueNorm>("BlueNorm");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::LpNorm1>("LpNorm1");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::LpNorm2>("LpNorm2");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::LpNorm3>("LpNorm3");
      add.template operator()<Op::Normalized>("Normalized");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Normalize>("Normalize");
      add.template operator()<Op::Cross>("Cross");

      if constexpr (Storage != SpatialStorage::Spatial)
      {
        add.template operator()<Op::Add, true>("Add");
        add.template operator()<Op::ScaleRight, true>("ScaleRight");
        add.template operator()<Op::Conjugate, true>("Conjugate");
        add.template operator()<Op::Transpose, true>("Transpose");
        add.template operator()<Op::MathSum, true>("MathSum");
        add.template operator()<Op::MathScale, true>("MathScale");
      }
    }
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  registerCases<Real, SpatialStorage::Spatial>("Spatial");
  registerCases<Real, SpatialStorage::Bounded>("Bounded");
  registerCases<Real, SpatialStorage::Dynamic>("Dynamic");
  registerCases<Real, SpatialStorage::Fixed, 1, 1>("Fixed");
  registerCases<Real, SpatialStorage::Fixed, 2, 1>("Fixed");
  registerCases<Real, SpatialStorage::Fixed, 3, 1>("Fixed");
  registerCases<Complex, SpatialStorage::Spatial>("Spatial");
  registerCases<Complex, SpatialStorage::Bounded>("Bounded");
  registerCases<Complex, SpatialStorage::Dynamic>("Dynamic");
  registerCases<Complex, SpatialStorage::Fixed, 1, 1>("Fixed");
  registerCases<Complex, SpatialStorage::Fixed, 2, 1>("Fixed");
  registerCases<Complex, SpatialStorage::Fixed, 3, 1>("Fixed");
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext(
    "timing", "Calling-thread CPU; input setup and numerical guards excluded");
  benchmark::AddCustomContext(
    "result", "Owned result assigned to retained output; dynamic allocations included");
  const auto count = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return count == 0 || failed ? 1 : 0;
}
