/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_HEART_CCMLC2014_PASSIVEENERGY_H
#define RODIN_HEART_CCMLC2014_PASSIVEENERGY_H

namespace Rodin::Heart
{
  template <class Scalar>
  struct ReducedInvariants
  {
    Scalar J1 = Scalar(0);
    Scalar J2 = Scalar(0);
    Scalar J4 = Scalar(0);
  };

  template <class Scalar>
  struct ReducedInvariantGradient
  {
    Scalar dW_dJ1 = Scalar(0);
    Scalar dW_dJ2 = Scalar(0);
    Scalar dW_dJ4 = Scalar(0);
  };

  template <class Scalar>
  struct ReducedInvariantHessian
  {
    Scalar d2W_dJ1dJ1 = Scalar(0);
    Scalar d2W_dJ1dJ2 = Scalar(0);
    Scalar d2W_dJ1dJ4 = Scalar(0);
    Scalar d2W_dJ2dJ2 = Scalar(0);
    Scalar d2W_dJ2dJ4 = Scalar(0);
    Scalar d2W_dJ4dJ4 = Scalar(0);
  };

  template <class Scalar>
  struct PassiveEnergyDerivatives
  {
    ReducedInvariantGradient<Scalar> grad;
    ReducedInvariantHessian<Scalar> hess;
  };
}

#endif
