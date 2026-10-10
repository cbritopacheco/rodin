/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MASSWORKLOAD_H
#define RODIN_MASSWORKLOAD_H

/** @file @brief Backend-independent mass/reaction workloads and constant-field oracles. */
#include <benchmark/benchmark.h>
#include "Rodin/Variational/H1.h"
#include "Rodin/Variational/P0g.h"
#include "../convergence/Convergence.h"
#include "../QuadratureReference.h"

namespace Rodin::Tests::Benchmarks
{
  enum class MassFamily
  {
    H1,
    P1,
    P0,
    P0g
  };

  /**
   * @brief Defines mass and reaction operators and projection right-hand sides.
   *
   * Architecture: a space is constructed through its public constructors;
   * constants are interpolated through DOF functionals, not coefficient layouts.
   * The forms are @f$a(u,v)=\int_\Omega\rho u\cdot\overline v\,dx@f$ and
   * @f$L(v)=a(c,v)@f$. On the unit box, @f$Mc=b@f$ and
   * @f$c^*Mc=\alpha(1+d/2)|c|^2@f$ for reaction; mass uses @f$\rho=1@f$.
   * Real reaction uses @f$\alpha=1@f$, complex reaction uses
   * @f$\alpha=1+i/2@f$. A referenced scale permits untimed coefficient rebinds.
   * Registration records three sizes and separates operator from load timing.
   */
  template <size_t K, class Scalar, bool Vector, MassFamily Family>
  class MassWorkload
  {
    public:
      template <class Mesh>
      static auto makeSpace(const Mesh& mesh, size_t dim)
      {
        using Range = std::conditional_t<Vector, Math::SpatialVector<Scalar>, Scalar>;
        using FES =
          std::conditional_t<Family == MassFamily::H1, Variational::H1<K, Range, Mesh>,
            std::conditional_t<Family == MassFamily::P1, Variational::P1<Range, Mesh>,
              std::conditional_t<Family == MassFamily::P0, Variational::P0<Range, Mesh>,
                Variational::P0g<Range, Mesh>>>>;
        if constexpr (Family == MassFamily::H1)
        {
          if constexpr (Vector)
            return FES(std::integral_constant<size_t, K>{}, mesh, dim);
          else
            return FES(std::integral_constant<size_t, K>{}, mesh);
        }
        else if constexpr (Vector)
          return FES(mesh, dim);
        else
          return FES(mesh);
      }

      template <bool Reference = false, class U, class V, class A, class B>
      static Scalar configure(
        U& u, V& v, A& matrix, B& load, size_t dim, bool reaction, const Real& scale)
      {
        using namespace Variational;
        const Scalar phase = []() {
          if constexpr (std::is_same_v<Scalar, Complex>)
            return Scalar(1, 2);
          else
            return Scalar(1);
        }();
        Scalar alpha = 1;
        if constexpr (std::is_same_v<Scalar, Complex>)
          if (reaction)
            alpha = Scalar(1, 0.5);
        auto field = [&]() {
          if constexpr (Vector)
            return VectorFunction(dim, [dim, phase](const Geometry::Point&) {
              Math::SpatialVector<Scalar> value(static_cast<std::uint8_t>(dim));
              for (size_t j = 0; j < dim; ++j)
                value(j) = Scalar(j + 1) * phase;
              return value;
            });
          else if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(phase);
          else
            return RealFunction(phase);
        }();
        auto coefficient = [&]() {
          const auto evaluate = [dim, reaction, alpha, &scale](const Geometry::Point& p) {
            Real rho = 1;
            if (reaction)
              for (size_t j = 0; j < dim; ++j)
                rho += p(j);
            return scale * alpha * rho;
          };
          if constexpr (std::is_same_v<Scalar, Complex>)
            return ComplexFunction(evaluate);
          else
            return RealFunction(evaluate);
        }();
        const auto assign = [&matrix](auto integral) {
          integral.setOrder(8);
          if constexpr (Reference)
            matrix = ReferenceIntegral(integral);
          else
            matrix = integral;
        };
        if (reaction)
          assign(Integral(u * coefficient, v));
        else
          assign(Integral(u, v));
        auto rhs = Integral(field * coefficient, v);
        rhs.setOrder(8);
        load = rhs;
        u.getSolution() = field;
        const Real square = Vector ? Real(dim * (dim + 1) * (2 * dim + 1)) / 6 : 1;
        const Real phaseSquare = std::abs(phase) * std::abs(phase);
        return alpha * square * phaseSquare * (reaction ? 1 + Real(dim) / 2 : 1);
      }

      template <class Callback>
      static void registerCases(Callback callback, bool manual)
      {
        using Geometry::Polytope;
        const std::string family = Family == MassFamily::H1 ? "H1/P" + std::to_string(K)
          : Family == MassFamily::P1                        ? "P1"
          : Family == MassFamily::P0                        ? "P0"
                                                            : "P0g";
        for (auto geometry :
          {Polytope::Type::Segment, Polytope::Type::Triangle,
            Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
            Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
        {
          for (bool reaction : {false, true})
          {
            for (bool rhs : {false, true})
            {
              const std::string name = std::string(reaction ? "Reaction/" : "Mass/") +
                (rhs ? "Load/" : "Operator/") + family + "/" +
                (std::is_same_v<Scalar, Complex> ? "Complex" : "Real") +
                (Vector ? "Vector/" : "Scalar/") +
                std::string(Convergence::UniformGrid::getGeometryName(geometry));
              auto* entry = benchmark::RegisterBenchmark(name.c_str(),
                [callback, geometry, reaction, rhs](
                  auto& state) { callback(state, geometry, reaction, rhs); })
                              ->Arg(3)
                              ->Arg(5)
                              ->Arg(9);
              if (manual)
                entry->Iterations(3)->UseManualTime();
              else
                entry->UseRealTime();
            }
          }
        }
      }
  };
}
#endif
