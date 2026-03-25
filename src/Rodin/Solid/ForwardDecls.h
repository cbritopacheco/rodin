/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ForwardDecls.h
 * @brief Forward declarations for the Solid mechanics module.
 */
#ifndef RODIN_SOLID_FORWARDDECLS_H
#define RODIN_SOLID_FORWARDDECLS_H

namespace Rodin::Solid
{
  class KinematicState;

  class ConstitutivePoint;

  namespace Tags
  {
    struct FiberDirection;
    struct SheetDirection;
    struct SheetNormalDirection;
    struct Activation;
  }

  template <class Derived>
  class Input;

  class IsotropicInvariants;

  template <class Derived>
  class HyperElasticLaw;

  class Hooke;

  template <class LawDerived, class FES>
  class InternalForce;

  template <class LawDerived, class Solution, class FES>
  class MaterialTangent;
}

#endif
