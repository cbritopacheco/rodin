/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_BENCHMARKS_PHYSICSFORM_H
#define RODIN_TESTS_BENCHMARKS_PHYSICSFORM_H

#include "Rodin/Variational.h"
#include "../QuadratureReference.h"

namespace Rodin::Tests::Benchmarks
{
  /** Defines the same physical form and analytic affine-field energy for every backend. */
  template <bool Elasticity>
  class PhysicsForm
  {
    public:
      template <bool Reference = false, class U, class V, class Form>
      static Real configure(U& u, V& v, Form& form, size_t dim, bool conductivity)
      {
        using namespace Variational;
        const auto wrap = [](const auto& integral) {
          if constexpr (Reference)
            return Tests::ReferenceIntegral(integral);
          else
            return integral;
        };
        if constexpr (Elasticity)
        {
          auto volumetric = Integral(1.5 * Div(u), Div(v));
          auto shear = Integral(
            0.5 * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
          volumetric.setOrder(6);
          shear.setOrder(6);
          form = wrap(volumetric) + wrap(shear);
          u.getSolution() = VectorFunction(dim, [dim](const Geometry::Point& p) {
            Real s = 0;
            for (size_t j = 0; j < dim; ++j)
              s += p(j);
            Math::SpatialVector<Real> value(static_cast<std::uint8_t>(dim));
            for (size_t i = 0; i < dim; ++i)
              value(i) = Real(i + 1) * s;
            return value;
          });
          const Real c = Real(dim * (dim + 1)) / 2;
          const Real squares = Real(dim * (dim + 1) * (2 * dim + 1)) / 6;
          return 1.5 * c * c + 0.5 * (Real(dim) * squares + c * c);
        }
        else
        {
          const RealFunction gamma([dim, conductivity](const Geometry::Point& p) {
            Real value = 1;
            if (conductivity)
              for (size_t j = 0; j < dim; ++j)
                value += p(j);
            return value;
          });
          auto diffusion = Integral(gamma * Grad(u), Grad(v));
          diffusion.setOrder(6);
          form = wrap(diffusion);
          u.getSolution() = RealFunction([dim](const Geometry::Point& p) {
            Real value = 0;
            for (size_t j = 0; j < dim; ++j)
              value += p(j);
            return value;
          });
          return Real(dim) * (conductivity ? 1 + Real(dim) / 2 : 1);
        }
      }
  };
}
#endif
