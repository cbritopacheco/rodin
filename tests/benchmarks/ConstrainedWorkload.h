/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_CONSTRAINEDWORKLOAD_H
#define RODIN_CONSTRAINEDWORKLOAD_H

#include <benchmark/benchmark.h>
#include "Rodin/Variational/H1.h"
#include "../convergence/Convergence.h"
#include "../QuadratureReference.h"

namespace Rodin::Tests::Benchmarks
{
  /**
   * @brief Manufactured scalar problems with nonhomogeneous full Dirichlet data.
   *
   * On the unit box, @f$u=\eta t(1+\sum_j(x_j+\chi x_j^2))@f$,
   * where @f$\chi=0@f$ at degree one and @f$\chi=1@f$ otherwise.
   * Poisson uses @f$\gamma=s@f$, conductivity uses
   * @f$\gamma=s(1+\sum_jx_j)@f$, and reaction--diffusion additionally
   * uses @f$\rho=s@f$. The source is derived from
   * @f$f=-\nabla\cdot(\gamma\nabla u)+\rho u@f$.
   * The exact field is interpolated through the space's functionals.
   * Referenced scales permit coefficient and boundary state transitions.
   * No solver or backend policy is defined here.
   */
  template <size_t K, class Scalar>
  class ConstrainedWorkload
  {
    public:
      static constexpr size_t quadratureOrder = K == 1 ? 8 : 12;

      static auto field(size_t dim, const Real& boundaryScale)
      {
        const auto evaluate = [dim, &boundaryScale](const Geometry::Point& p) {
          Real value = 1;
          for (size_t j = 0; j < dim; ++j)
            value += p(j) + (K > 1 ? p(j) * p(j) : 0);
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 2) * boundaryScale * value;
          else
            return boundaryScale * value;
        };
        if constexpr (std::is_same_v<Scalar, Complex>)
          return Variational::ComplexFunction(evaluate);
        else
          return Variational::RealFunction(evaluate);
      }

      template <bool Reference = false, class U, class V, class ProblemType>
      static void configure(U& u, V& v, ProblemType& problem, size_t dim,
        unsigned physics, const Real& operatorScale, const Real& boundaryScale)
      {
        using namespace Variational;
        auto exact = field(dim, boundaryScale);
        RealFunction gamma([dim, physics, &operatorScale](const Geometry::Point& p) {
          Real value = 1;
          if (physics != 0)
            for (size_t j = 0; j < dim; ++j)
              value += p(j);
          return operatorScale * value;
        });
        const auto source = [dim, physics, &operatorScale, &boundaryScale](
                              const Geometry::Point& p) {
          Real sum = 0, derivative = 0, value = 1;
          for (size_t j = 0; j < dim; ++j)
          {
            sum += p(j);
            derivative += 1 + (K > 1 ? 2 * p(j) : 0);
            value += p(j) + (K > 1 ? p(j) * p(j) : 0);
          }
          const Real gammaValue = physics != 0 ? 1 + sum : 1;
          const Real forcing = -(K > 1 ? 2 * Real(dim) * gammaValue : 0) -
            (physics != 0 ? derivative : 0) + (physics == 2 ? value : 0);
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 2) * operatorScale * boundaryScale * forcing;
          else
            return operatorScale * boundaryScale * forcing;
        };
        auto forcing = [&]() {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(source);
          else
            return RealFunction(source);
        }();
        RealFunction rho([physics, &operatorScale](const Geometry::Point&) {
          return physics == 2 ? operatorScale : 0;
        });
        auto diffusion = Integral(gamma * Grad(u), Grad(v));
        auto reaction = Integral(rho * u, v);
        auto load = Integral(forcing, v);
        diffusion.setOrder(quadratureOrder);
        reaction.setOrder(quadratureOrder);
        load.setOrder(quadratureOrder);
        if constexpr (Reference)
        {
          if (physics == 2)
            problem = ReferenceIntegral(diffusion) + ReferenceIntegral(reaction) - load +
              DirichletBC(u, exact);
          else
            problem = ReferenceIntegral(diffusion) - load + DirichletBC(u, exact);
        }
        else
        {
          if (physics == 2)
            problem = diffusion + reaction - load + DirichletBC(u, exact);
          else
            problem = diffusion - load + DirichletBC(u, exact);
        }
        u.getSolution() = exact;
      }

      template <class Callback>
      static void registerCases(Callback callback, bool manual)
      {
        using Geometry::Polytope;
        for (auto geometry :
          {Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
          for (unsigned physics : {0u, 1u, 2u})
          {
            const std::string name = std::string(physics == 0 ? "Poisson"
                                         : physics == 1       ? "Conductivity"
                                                              : "ReactionDiffusion") +
              "/Constrained/H1/P" + std::to_string(K) + "/" +
              (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
              std::string(Convergence::UniformGrid::getGeometryName(geometry));
            auto* entry = benchmark::RegisterBenchmark(name.c_str(),
              [callback, geometry, physics](
                auto& state) { callback(state, geometry, physics); })
                            ->Arg(3)
                            ->Arg(5)
                            ->Arg(9);
            if (manual)
              entry->Iterations(3)->UseManualTime();
            else
              entry->UseRealTime();
          }
      }
  };
}
#endif
