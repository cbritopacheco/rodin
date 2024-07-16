/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_SIGNEDPOISSON_H
#define RODIN_MODELS_DISTANCE_SIGNEDPOISSON_H

#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Poisson approximation to the signed distance function.
   */
  class SignedPoisson
  {
    public:
      template <class FES>
      auto operator()(const FES& fes) const
      {
        const auto& mesh = fes.getMesh();
        const auto& mdim = mesh.getDimension();
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, mdim - 1, mdim);
        Variational::TrialFunction u(fes);
        Variational::TestFunction  v(fes);
        Variational::RealFunction zero(0);
        Variational::Problem sp(u, v);
        sp = Variational::Integral(Variational::Grad(u), Variational::Grad(v))
           + Variational::Integral(v)
           + Variational::DirichletBC(u, zero);
        Solver::CG(sp).solve();
        return u.getSolution();
      }

      template <class FES>
      auto operator()(
          Geometry::Attribute interface,
          Geometry::Attribute region,
          const FES& fes) const
      {
        return operator()(
            FlatSet<Geometry::Attribute>{interface}, FlatSet<Geometry::Attribute>{region}, fes);
      }

      template <class FES>
      auto operator()(
          const FlatSet<Geometry::Attribute>& interface,
          const FlatSet<Geometry::Attribute>& region,
          const FES& fes) const
      {
        const auto& mesh = fes.getMesh();
        Variational::TrialFunction u(fes);
        Variational::TestFunction  v(fes);
        Variational::RealFunction zero = 0;
        Variational::Problem sp(u, v);
        sp = Variational::Integral(Variational::Grad(u), Variational::Grad(v))
           + Variational::Integral(2 * v).over(region)
           - Variational::Integral(v)
           + Variational::DirichletBC(u, zero).on(interface);
        Solver::CG(sp).solve();
        return u.getSolution();
      }
  };
}

#endif

