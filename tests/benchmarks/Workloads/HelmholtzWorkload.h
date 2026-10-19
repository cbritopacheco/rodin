/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_WORKLOADS_HELMHOLTZWORKLOAD_H
#define RODIN_WORKLOADS_HELMHOLTZWORKLOAD_H

#include <benchmark/benchmark.h>
#include "convergence/Helmholtz.h"
#include "QuadratureReference.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Shared complex Helmholtz forms for assembly and solver measurements.
   *
   * For @f$-\Delta u+\alpha u=f@f$ on the unit box, full Dirichlet data
   * and the affine/quadratic exact field come from
   * @ref Convergence::HelmholtzData. The load is shifted analytically from
   * its @f$\alpha=-1/4@f$ source. The second case uses
   * @f$\alpha=-64+i@f$, whose continuous real part is indefinite in
   * dimensions one through three. Coarse discrete interior spectra need not
   * have both signs. No assembly backend or solver policy is selected here.
   *
   * @tparam K Degree of the complex H1 space.
   */
  template <size_t K>
  class HelmholtzWorkload
  {
    public:
      /// @brief Integrates polynomial patches and rational pyramid basis products.
      static constexpr size_t QuadratureOrder = 12;

      /// @brief Above the first unit-box Dirichlet eigenvalue in dimensions 1--3.
      static constexpr Real ShiftedWaveNumberSquared = 64;

      /// @brief Positive imaginary physical mass shift, not a solver parameter.
      static constexpr Real Absorption = 1;

      /**
       * @brief Returns the mass coefficient for the selected physical operator.
       * @param shifted Whether the complex-shifted case is selected.
       * @returns Complex mass coefficient.
       */
      static Complex getMassCoefficient(bool shifted)
      {
        return shifted ? Complex(-ShiftedWaveNumberSquared, Absorption)
                       : Complex(-Convergence::HelmholtzData::WaveNumberSquared);
      }

      /**
       * @brief Defines either the optimized or original local-kernel formulation.
       * @tparam Reference Whether original integral kernels are requested.
       * @tparam U Trial function type.
       * @tparam V Test function type.
       * @tparam ProblemType Native complex problem type.
       * @param u Trial function; its coefficients receive the exact interpolant.
       * @param v Test function.
       * @param problem Destination formulation.
       * @param data Manufactured physical field and low-wave-number source.
       * @param alpha Complex mass coefficient, retained by value in the source.
       */
      template <bool Reference = false, class U, class V, class ProblemType>
      static void configure(U& u, V& v, ProblemType& problem,
        const Convergence::HelmholtzData& data, Complex alpha)
      {
        using namespace Variational;
        const auto exact = data.getSolution();
        const auto lowSource = data.getSource();
        ComplexFunction source([exact, lowSource, alpha](const Geometry::Point& p) {
          return lowSource(p) +
            (alpha + Convergence::HelmholtzData::WaveNumberSquared) * exact(p);
        });
        auto diffusion = Integral(Grad(u), Grad(v));
        auto mass = Integral(alpha * u, v);
        auto load = Integral(source, v);
        diffusion.setOrder(QuadratureOrder);
        mass.setOrder(QuadratureOrder);
        load.setOrder(QuadratureOrder);
        if constexpr (Reference)
          problem = ReferenceIntegral(diffusion) + ReferenceIntegral(mass) - load +
            DirichletBC(u, exact);
        else
          problem = diffusion + mass - load + DirichletBC(u, exact);
        u.getSolution() = exact;
      }

      /**
       * @brief Registers matched physics, degrees, geometries and grid sizes.
       * @tparam Callback Callable accepting state, geometry and shifted-case flag.
       * @param scope Owning operation name, such as Assembly or GMRES.
       * @param callback Operation to measure.
       */
      template <class Callback>
      static void registerCases(const char* scope, Callback callback)
      {
        using Geometry::Polytope;
        for (const auto geometry :
          {Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (const bool shifted : {false, true})
          {
            const std::string name = std::string("Helmholtz/") + scope + "/" +
              (shifted ? "Shifted/" : "LowWave/") + "H1/P" + std::to_string(K) +
              "/Complex/" +
              std::string(Convergence::UniformGrid::getGeometryName(geometry));
            benchmark::RegisterBenchmark(name.c_str(),
              [callback, geometry, shifted](
                auto& state) { callback(state, geometry, shifted); })
              ->Arg(3)
              ->Arg(5)
              ->Arg(9)
              ->UseRealTime();
          }
        }
      }
  };
}
#endif
