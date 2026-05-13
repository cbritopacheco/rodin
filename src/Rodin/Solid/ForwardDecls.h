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
    struct ActiveExtension;
    struct PreviousActiveExtension;
    struct TimeStep;
    struct ElectricalActivation;
    struct PreviousActiveGamma;
    struct PreviousActiveBeta;
  }

  template <class Derived>
  class Input;

  class IsotropicInvariants;

  class FiberKinematics;

  template <class Derived>
  class HyperElasticLaw;

  class Hooke;

  class HolzapfelOgden;

  class ActiveFiberLaw;

  template <class PassiveLaw, class ActiveLaw>
  class ActiveContraction;

  template <class LawDerived, class TestFunctionType, class DisplacementType>
  class InternalForce;

  template <class LawDerived, class TrialFunctionType,
            class TestFunctionType, class DisplacementType>
  class MaterialTangent;
}

#endif
