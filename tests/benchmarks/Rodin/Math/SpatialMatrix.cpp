/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Matched SpatialMatrix arithmetic across spatial and Eigen representations.
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
      using Fixture = SpatialArithmetic<Scalar, true, Storage, Rows, Cols>;
      using Op = SpatialOperation;
      const std::string prefix = "SpatialMatrix/" +
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
          for (int cols = 0; cols <= 3; ++cols)
          {
            if constexpr (Storage == SpatialStorage::Fixed)
            {
              if (rows != Rows || cols != Cols)
                continue;
            }
            if constexpr (Operation == Op::Trace || Operation == Op::Determinant ||
              Operation == Op::Inverse || Operation == Op::Solve)
            {
              if (rows != cols)
                continue;
            }
            if constexpr (Operation == Op::Inverse || Operation == Op::Solve ||
              Operation == Op::PseudoInverse || Operation == Op::Row ||
              Operation == Op::Column)
            {
              if (!rows || !cols)
                continue;
            }
            for (int count : {1, 32})
              item->Args({rows, cols, count});
          }
        }
        item->ArgNames({"rows", "cols", "count"});
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
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::MathProduct>("MathProduct");
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
        add.template operator()<Op::Adjoint>("Adjoint");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Trace>("Trace");
      add.template operator()<Op::Determinant>("Determinant");
      add.template operator()<Op::Inverse>("Inverse");
      add.template operator()<Op::Solve>("Solve");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::PseudoInverse>("PseudoInverse");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Row>("Row");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::Column>("Column");
      add.template operator()<Op::MatrixVector>("MatrixVector");
      if constexpr (Storage != SpatialStorage::Fixed)
        add.template operator()<Op::VectorMatrix>("VectorMatrix");
      add.template operator()<Op::MatrixProduct>("MatrixProduct");

      if constexpr (Storage != SpatialStorage::Spatial)
      {
        add.template operator()<Op::Add, true>("Add");
        add.template operator()<Op::ScaleRight, true>("ScaleRight");
        add.template operator()<Op::Conjugate, true>("Conjugate");
        add.template operator()<Op::Transpose, true>("Transpose");
        add.template operator()<Op::MathSum, true>("MathSum");
        add.template operator()<Op::MathScale, true>("MathScale");
        add.template operator()<Op::MatrixProduct, true>("MatrixProduct");
        add.template operator()<Op::MatrixVector, true>("MatrixVector");
        add.template operator()<Op::VectorMatrix, true>("VectorMatrix");
        add.template operator()<Op::MathProduct, true>("MathProduct");
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
  registerCases<Real, SpatialStorage::Fixed, 2, 2>("Fixed");
  registerCases<Real, SpatialStorage::Fixed, 3, 3>("Fixed");
  registerCases<Complex, SpatialStorage::Spatial>("Spatial");
  registerCases<Complex, SpatialStorage::Bounded>("Bounded");
  registerCases<Complex, SpatialStorage::Dynamic>("Dynamic");
  registerCases<Complex, SpatialStorage::Fixed, 1, 1>("Fixed");
  registerCases<Complex, SpatialStorage::Fixed, 2, 2>("Fixed");
  registerCases<Complex, SpatialStorage::Fixed, 3, 3>("Fixed");
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
