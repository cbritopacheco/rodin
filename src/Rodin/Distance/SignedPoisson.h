/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file SignedPoisson.h
 * @brief Signed Poisson equation-based distance function approximation.
 *
 * This file provides a signed distance function computation method using
 * the Poisson equation with region-specific forcing terms.
 */
#ifndef RODIN_MODELS_DISTANCE_SIGNEDPOISSON_H
#define RODIN_MODELS_DISTANCE_SIGNEDPOISSON_H

#include "Rodin/Solver/CG.h"

#include "Rodin/Variational/Problem.h"
#include "Rodin/Variational/Integral.h"
#include "Rodin/Variational/DirichletBC.h"
#include "Rodin/Variational/GridFunction.h"

#include "Rodin/Alert/MemberFunctionException.h"

namespace Rodin::Distance
{
  /**
   * @brief Poisson approximation to the signed distance function.
   *
   * This class computes an approximation to the signed distance function by
   * solving a Poisson equation with different forcing terms in interior and
   * exterior regions:
   * @f[
   *   -\Delta u = f \quad \text{in } \Omega
   * @f]
   * @f[
   *   u = 0 \quad \text{on } \Gamma
   * @f]
   * where @f$ f @f$ changes sign based on the region, and @f$ \Gamma @f$ is
   * the interface.
   *
   * ## Variants
   * This class provides two operator() overloads:
   * 1. Unsigned distance (full domain)
   * 2. Signed distance (with specified interior/interface regions)
   *
   * ## Usage Example
   * ```cpp
   * P1 fes(mesh);
   * SignedPoisson signedDist;
   * // Signed distance with regions:
   * auto u = signedDist(interfaceAttr, interiorAttr, fes);
   * ```
   */
  class SignedPoisson
  {
    public:
      /**
       * @brief Computes the unsigned Poisson-based distance approximation.
       *
       * Solves the Poisson equation over the entire domain with uniform
       * forcing term.
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
           + Variational::Integral(v)
           + Variational::DirichletBC(u, zero);
        Solver::CG(sp).solve();
        return u.getSolution();
      }

      /**
       * @brief Computes the signed distance for specified interface and region.
       *
       * Convenience overload that accepts single attributes for interface
       * and region.
       *
       * @tparam FES Finite element space type
       * @param[in] interface Attribute marking the interface (zero level set)
       * @param[in] region Attribute marking the interior region
       * @param[in] fes Finite element space on which to solve
       * @return Grid function containing the signed distance approximation
       */
      template <class FES>
      auto operator()(
          Geometry::Attribute interface,
          Geometry::Attribute region,
          const FES& fes) const
      {
        return operator()(
            FlatSet<Geometry::Attribute>{interface}, FlatSet<Geometry::Attribute>{region}, fes);
      }

      /**
       * @brief Computes the signed distance for specified interface and region sets.
       *
       * Solves the Poisson equation with different forcing terms in the interior
       * and exterior regions, creating a signed distance approximation.
       *
       * @tparam FES Finite element space type
       * @param[in] interface Set of attributes marking the interface
       * @param[in] region Set of attributes marking the interior region
       * @param[in] fes Finite element space on which to solve
       * @return Grid function containing the signed distance approximation
       */
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

