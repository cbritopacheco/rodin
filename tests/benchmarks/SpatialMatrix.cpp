/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief CPU-time attribution of matrix assembly versus nine scalar assemblies. */
#include <benchmark/benchmark.h>
#include "Rodin/Assembly.h"
#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Variational.h"
using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
namespace
{
  // Stage: 0 complete assembly, 1 triplets, 2 element binding, 3 cached entry
  // scan on one cell, 4 sparse conversion from precomputed triplets.
  template <size_t K, size_t Stage, bool Threaded = false>
  void matrixAssembly(benchmark::State& state)
  {
    const bool threeD = state.range(0) == 3;
    auto mesh = threeD ? LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3})
                       : LocalMesh::UniformGrid(Polytope::Type::Triangle, {8, 8});
    for (size_t d = 1; d <= mesh.getDimension(); ++d)
      for (size_t l = 0; l < d; ++l)
        mesh.getConnectivity().compute(d, l);
    auto run = [&](const auto& fes, size_t repetitions) {
      using Space = std::remove_cvref_t<decltype(fes)>;
      using Triplets = std::vector<Math::SparseTriplet<Real>>;
      TrialFunction u(fes);
      TestFunction v(fes);
      BilinearForm form(u, v);
      if (state.range(1))
        form = Integral(Grad(u), Grad(v));
      else
        form = Integral(u, v);
      using Form = decltype(form);
      using TripletForm =
        BilinearForm<typename Form::SolutionType, Space, Space, Triplets>;
      Assembly::Sequential<Triplets, TripletForm> tripletAssembly;
      Assembly::Sequential<typename Form::OperatorType, Form> sequential;
      Triplets triplets;
      auto input = [&] {
        return typename decltype(tripletAssembly)::InputType{
          fes, fes, form.getLocalIntegrators(), form.getGlobalIntegrators()};
      };
      // Warm geometry/tabulation caches and retain the numerical reference.
      sequential.execute(form.getOperator(),
        {fes, fes, form.getLocalIntegrators(), form.getGlobalIntegrators()});
      const Math::SparseMatrix<Real> reference = form.getOperator();
      if constexpr (FormLanguage::IsMatrixRange<
                      typename FormLanguage::Traits<Space>::RangeType>::Value)
      {
        const auto& scalarSpace = fes.getScalarSpace();
        TrialFunction scalarTrial(scalarSpace);
        TestFunction scalarTest(scalarSpace);
        BilinearForm scalarForm(scalarTrial, scalarTest);
        if (state.range(1))
          scalarForm = Integral(Grad(scalarTrial), Grad(scalarTest));
        else
          scalarForm = Integral(scalarTrial, scalarTest);
        Assembly::Sequential<typename decltype(scalarForm)::OperatorType,
          decltype(scalarForm)>
          scalarAssembly;
        scalarAssembly.execute(scalarForm.getOperator(),
          {scalarSpace, scalarSpace, scalarForm.getLocalIntegrators(),
            scalarForm.getGlobalIntegrators()});
        const auto& scalar = scalarForm.getOperator();
        // This comparison uses the documented component ordering of H1 matrix
        // spaces, not an assumed ordering for arbitrary finite element families.
        Triplets expectedTriplets;
        expectedTriplets.reserve(9 * scalar.nonZeros());
        for (Index outer = 0; outer < static_cast<Index>(scalar.outerSize()); ++outer)
          for (Math::SparseMatrix<Real>::InnerIterator entry(scalar, outer); entry;
               ++entry)
            for (Index component = 0; component < 9; ++component)
              expectedTriplets.emplace_back(
                9 * entry.row() + component, 9 * entry.col() + component, entry.value());
        Math::SparseMatrix<Real> expected(reference.rows(), reference.cols());
        expected.setFromTriplets(expectedTriplets.begin(), expectedTriplets.end());
        const Real difference = (reference - expected).norm();
        state.counters["scalar_block_relative_error"] =
          difference / std::max(Real(1), expected.norm());
        // Cancellation can retain different near-zero sparse entries. Compare
        // numerical operators rather than requiring identical storage patterns.
        if (difference > 1e-12 * std::max(Real(1), expected.norm()))
        {
          state.SkipWithError("Matrix assembly differs from nine scalar blocks.");
          return;
        }
      }
      tripletAssembly.execute(triplets, input());
      Math::SparseMatrix<Real> converted(fes.getSize(), fes.getSize());
      converted.setFromTriplets(triplets.begin(), triplets.end());
      if ((converted - reference).norm() > 1e-12 * std::max(Real(1), reference.norm()))
      {
        state.SkipWithError("Triplet assembly differs from sparse reference.");
        return;
      }
      const size_t dimension = mesh.getDimension();
      const auto cell = mesh.getPolytope(dimension, 0);
      const size_t localDOFs = fes.getDOFs(dimension, 0).size();
      // Numerical work counters refer to one full nine-component problem.
      state.counters["dofs"] = fes.getSize() * repetitions;
      state.counters["local_entries"] = localDOFs * localDOFs * repetitions;
      state.counters["triplets"] = triplets.size() * repetitions;
      state.counters["nonzeros"] = reference.nonZeros() * repetitions;
      state.counters["cells"] = mesh.getPolytopeCount(dimension);
      for (auto& integral : form.getLocalIntegrators())
        integral.setPolytope(*cell);
      for (auto iteration : state)
      {
        for (size_t repeat = 0; repeat < repetitions; ++repeat)
        {
          if constexpr (Stage == 0)
          {
            if constexpr (Threaded)
              form.assemble();
            else
              sequential.execute(form.getOperator(),
                {fes, fes, form.getLocalIntegrators(), form.getGlobalIntegrators()});
            benchmark::DoNotOptimize(form.getOperator().valuePtr());
          }
          else if constexpr (Stage == 1)
          {
            tripletAssembly.execute(triplets, input());
            benchmark::DoNotOptimize(triplets.data());
          }
          else if constexpr (Stage == 2)
          {
            for (auto& integral : form.getLocalIntegrators())
              for (Index i = 0; i < mesh.getPolytopeCount(dimension); ++i)
                integral.setPolytope(*mesh.getPolytope(dimension, i));
          }
          else if constexpr (Stage == 3)
          {
            Real checksum = 0;
            for (auto& integral : form.getLocalIntegrators())
              for (size_t i = 0; i < localDOFs; ++i)
                for (size_t j = 0; j < localDOFs; ++j)
                  checksum += integral.integrate(j, i);
            benchmark::DoNotOptimize(checksum);
          }
          else
          {
            converted.setFromTriplets(triplets.begin(), triplets.end());
            benchmark::DoNotOptimize(converted.valuePtr());
          }
        }
        benchmark::ClobberMemory();
      }
      if constexpr (Stage == 0)
        if ((form.getOperator() - reference).norm() >
          1e-12 * std::max(Real(1), reference.norm()))
          state.SkipWithError("Timed assembly differs from sequential reference.");
    };
    if (state.range(2))
      run(H1(std::integral_constant<size_t, K>{}, mesh), 9);
    else
      run(H1(std::integral_constant<size_t, K>{}, mesh, 3, 3), 1);
  }

  template <size_t K>
  void registerAssemblyCases()
  {
    const std::string prefix = "SpatialAssembly/P" + std::to_string(K) + "/";
    // Arguments: dimension, stiffness (otherwise mass), nine scalar assemblies
    // (otherwise one native matrix assembly). CPU time includes all worker threads.
    benchmark::RegisterBenchmark(
      (prefix + "Full/Sequential").c_str(), &matrixAssembly<K, 0>)
      ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
      ->MeasureProcessCPUTime();
    benchmark::RegisterBenchmark(
      (prefix + "Triplets/Sequential").c_str(), &matrixAssembly<K, 1>)
      ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
      ->MeasureProcessCPUTime();
    benchmark::RegisterBenchmark(
      (prefix + "Bind/Sequential").c_str(), &matrixAssembly<K, 2>)
      ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
      ->MeasureProcessCPUTime();
    benchmark::RegisterBenchmark(
      (prefix + "Scan/Sequential").c_str(), &matrixAssembly<K, 3>)
      ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
      ->MeasureProcessCPUTime();
    benchmark::RegisterBenchmark(
      (prefix + "Sparse/Sequential").c_str(), &matrixAssembly<K, 4>)
      ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
      ->MeasureProcessCPUTime();
#ifdef RODIN_USE_OPENMP
    benchmark::RegisterBenchmark(
      (prefix + "Full/OpenMP").c_str(), &matrixAssembly<K, 0, true>)
      ->ArgsProduct({{2, 3}, {0, 1}, {0, 1}})
      ->MeasureProcessCPUTime();
#endif
  }
  const bool Registered = [] {
    registerAssemblyCases<1>();
    registerAssemblyCases<2>();
    registerAssemblyCases<3>();
    return true;
  }();
}
