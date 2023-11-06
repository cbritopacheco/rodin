/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_DISTANCE_POISSON_H
#define RODIN_MODELS_DISTANCE_POISSON_H

#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Models::Distance
{
  /**
   * @brief Poisson approximation to the distance function.
   */
  class Poisson
  {
    public:
      template <class FES>
      auto operator()(const FES& fes) const
      {
        const auto& mesh = fes.getMesh();
        const auto& mdim = mesh.getDimension();
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, mdim - 1, mdim);
        Solver::CG cg;
        Variational::TrialFunction u(fes);
        Variational::TestFunction  v(fes);
        Variational::ScalarFunction zero(0);
        Variational::Problem sp(u, v);
        sp = Variational::Integral(Variational::Grad(u), Variational::Grad(v))
           - Variational::Integral(v)
           + Variational::DirichletBC(u, zero);
        sp.solve(cg);
        return u.getSolution();
      }
  };
}

#endif


