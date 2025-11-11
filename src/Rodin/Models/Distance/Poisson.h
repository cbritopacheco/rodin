/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Poisson.h
 * @brief Poisson equation-based distance function approximation.
 *
 * This file provides a distance function computation method using
 * the Poisson equation as an approximation.
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
   *
   * This class computes an approximation to the distance function by solving
   * a Poisson equation with appropriate boundary conditions. The method solves:
   * @f[
   *   -\Delta u = 1 \quad \text{in } \Omega
   * @f]
   * @f[
   *   u = 0 \quad \text{on } \partial\Omega
   * @f]
   * where @f$ u @f$ approximates the distance to the boundary.
   *
   * ## Mathematical Background
   * Unlike the Eikonal equation which gives exact distances, the Poisson
   * approach provides a smooth approximation that can be computed more
   * efficiently for certain applications.
   *
   * ## Usage Example
   * ```cpp
   * P1 fes(mesh);
   * Poisson poissonDist;
   * auto u = poissonDist(fes);
   * ```
   */
  class Poisson
  {
    public:
      /**
       * @brief Computes the Poisson-based distance approximation.
       *
       * @tparam FES Finite element space type
       * @param[in] fes Finite element space on which to solve
       * @return Grid function containing the distance approximation
       */
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
           - Variational::Integral(v)
           + Variational::DirichletBC(u, zero);
        Solver::CG(sp).solve();
        return u.getSolution();
      }
  };
}

#endif


