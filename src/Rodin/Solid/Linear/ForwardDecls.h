/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the linear (small-strain) elasticity
 * integrators.
 */
#ifndef RODIN_SOLID_LINEAR_FORWARDDECLS_H
#define RODIN_SOLID_LINEAR_FORWARDDECLS_H

namespace Rodin::Variational
{
  template <class Solution, class FES, class LambdaDerived, class MuDerived>
  class LinearElasticityIntegrator;

  template <class Solution, class FES>
  class LinearElasticityIntegral;
}

#endif
