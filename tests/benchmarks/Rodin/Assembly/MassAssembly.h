/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Full Eigen mass/reaction operator and projection-load assembly.
 */
#ifndef RODIN_RODIN_ASSEMBLY_MASSASSEMBLY_H
#define RODIN_RODIN_ASSEMBLY_MASSASSEMBLY_H

#include "Rodin/Assembly.h"
#include "Workloads/MassWorkload.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  inline bool massCheckFailed = false;

  /**
   * @brief Times one public operator or load assembly, with independent oracles.
   *
   * ## Architecture
   *
   * Setup, independent geometry, original-loop assembly, constant-field energy and
   * projection identities are untimed. Repeated assembly must replace output.
   * A doubled reaction coefficient is checked after timing to detect stale
   * expression values. No projection solve or constraint elimination is timed.
   * The workload parameter supplies the form without duplicating assembly or
   * qualification. Its default retains the original right-weighted reaction.
   *
   * @tparam K Field polynomial degree.
   * @tparam Scalar Real or complex scalar field.
   * @tparam Vector Whether the field has one component per spatial dimension.
   * @tparam Family Finite-element family.
   * @tparam Workload Space, geometry, form and independent action definition.
   */
  template <size_t K, class Scalar, bool Vector, MassFamily Family,
    class Workload = MassWorkload<K, Scalar, Vector, Family>>
  class MassAssembly
  {
    public:
      /// @brief Relative budget for the independently integrated constant-field energy.
      static constexpr Real EnergyBudget = 1e-10;

      /// @brief Relative budget for differently associated original-entry assembly.
      static constexpr Real OperatorBudget = 1e-11;

      /**
       * @brief Measures one affine or curved operator/load workload.
       * @param state Benchmark iterations and grid size.
       * @param geometry Reference-cell family.
       * @param reaction Whether to use variable reaction rather than unit mass.
       * @param rhs Whether to time the projection load rather than the operator.
       * @param curved Whether to use the prescribed quadratic map.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        bool reaction, bool rhs, bool curved)
      {
        using namespace Variational;
        auto mesh = Workload::makeMesh(geometry, state.range(0), curved);
        const size_t dim = mesh.getDimension();
        const bool geometryValid = Workload::checkGeometry(mesh, geometry, curved);
        if (!geometryValid)
        {
          massCheckFailed = true;
          state.SkipWithError("Mapped quadrature differs from pointwise geometry or "
                              "loses quadratic degree");
          return;
        }
        auto space = Workload::makeSpace(mesh, dim);
        TrialFunction u(space);
        TestFunction v(space);
        BilinearForm matrix(u, v);
        LinearForm load(v);
        Real scale = 1;
        const Scalar expected =
          Workload::configure(u, v, matrix, load, dim, reaction, scale, curved);
        matrix.assemble();
        load.assemble();
        BilinearForm reference(u, v);
        LinearForm referenceLoad(v);
        Workload::template configure<true>(
          u, v, reference, referenceLoad, dim, reaction, scale, curved);
        reference.assemble();
        const auto baseline = reference.getOperator();
        const auto baselineLoad = load.getVector();
        const auto& c = u.getSolution().getData();
        const auto valid = [&](const auto& a, const auto& b) {
          const Scalar energy = Math::dot(a * c, c);
          return std::isfinite(std::abs(energy)) &&
            std::abs(energy / (scale * expected) - Scalar(1)) < EnergyBudget &&
            (a - scale * baseline).norm() <
            OperatorBudget * std::max(Real(1), baseline.norm() * scale) &&
            (a * c - b).norm() < OperatorBudget * std::max(Real(1), b.norm()) &&
            (b - scale * baselineLoad).norm() <
            OperatorBudget * std::max(Real(1), b.norm());
        };
        matrix.assemble();
        load.assemble();
        if (!valid(matrix.getOperator(), load.getVector()))
        {
          massCheckFailed = true;
          state.SkipWithError("Mass/reaction operator, load, or analytic oracle failed");
          return;
        }
        {
          auto perturbed = matrix.getOperator();
          perturbed *= Scalar(2);
          if (valid(perturbed, load.getVector()))
          {
            massCheckFailed = true;
            state.SkipWithError("Qualification accepted a doubled operator");
            return;
          }
        }
        {
          const auto perturbed = (Scalar(2) * load.getVector()).eval();
          if (valid(matrix.getOperator(), perturbed))
          {
            massCheckFailed = true;
            state.SkipWithError("Qualification accepted a doubled projection load");
            return;
          }
        }
        for (auto _ : state)
        {
          if (rhs)
          {
            load.assemble();
            benchmark::DoNotOptimize(load.getVector().data());
          }
          else
          {
            matrix.assemble();
            benchmark::DoNotOptimize(matrix.getOperator().nonZeros());
          }
          benchmark::ClobberMemory();
        }
        if (!valid(matrix.getOperator(), load.getVector()))
        {
          massCheckFailed = true;
          state.SkipWithError("Timed assembly changed operator or load");
        }
        if (reaction)
        {
          scale = 2;
          matrix.assemble();
          load.assemble();
          if (!valid(matrix.getOperator(), load.getVector()))
          {
            massCheckFailed = true;
            state.SkipWithError("Coefficient update did not replace operator and load");
          }
        }
        state.counters["cells"] = mesh.getPolytopeCount(dim);
        state.counters["dimension"] = dim;
        state.counters["ambient_dimension"] = mesh.getSpaceDimension();
        state.counters["points_per_axis"] = state.range(0);
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = matrix.getOperator().nonZeros();
        state.counters["degree"] = K;
        state.counters["geometry_degree"] = curved ? 2 : 1;
        state.counters["geometry_amplitude"] = curved ? Workload::GeometryAmplitude : 0;
        state.counters["domain_measure"] =
          dim == 1 && curved ? 1 + Workload::GeometryAmplitude : 1;
        state.counters["components"] = space.getVectorDimension();
        const auto firstCell = mesh.getCell();
        assert(firstCell);
        const size_t localBasis = u.getDOFs(*firstCell);
        state.counters["local_basis"] = localBasis;
        state.counters["local_entries"] = localBasis * localBasis;
        state.counters["quadrature_order"] = 8;
        state.counters["quadrature_points_per_cell"] =
          QF::QuadratureFormula::get(8, geometry).getSize();
      }

      static void registerCases()
      {
        Workload::registerCases(run, false);
      }
  };

  template <class Scalar, bool Vector>
  void registerMass()
  {
    MassAssembly<1, Scalar, Vector, MassFamily::H1>::registerCases();
    MassAssembly<2, Scalar, Vector, MassFamily::H1>::registerCases();
    MassAssembly<3, Scalar, Vector, MassFamily::H1>::registerCases();
    MassAssembly<1, Scalar, Vector, MassFamily::P1>::registerCases();
    MassAssembly<0, Scalar, Vector, MassFamily::P0>::registerCases();
    MassAssembly<0, Scalar, Vector, MassFamily::P0g>::registerCases();
  }
}

#endif
