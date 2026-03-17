/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HYPERELASTICITY_FORWARDDECLS_H
#define RODIN_HYPERELASTICITY_FORWARDDECLS_H

namespace Rodin::Hyperelasticity
{
  template <class Derived>
  class ConstitutiveLawBase;

  class NeoHookean;

  template <class ConstitutiveLaw>
  class HyperelasticProblem;
}

#endif
