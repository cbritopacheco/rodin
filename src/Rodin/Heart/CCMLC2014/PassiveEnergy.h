/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file PassiveEnergy.h
 * @brief Reduced invariant containers and passive-energy derivatives for CCMLC2014.
 */
#ifndef RODIN_HEART_CCMLC2014_PASSIVEENERGY_H
#define RODIN_HEART_CCMLC2014_PASSIVEENERGY_H

namespace Rodin::Heart
{
  /**
   * @brief Reduced invariants used by the passive-energy law.
   *
   * For the 0D reduction, these are scalarized forms of the 3D invariants.
   */
  template <class Scalar>
  struct ReducedInvariants
  {
    Scalar J1 = Scalar(0); ///< First reduced invariant.
    Scalar J2 = Scalar(0); ///< Second reduced invariant.
    Scalar J4 = Scalar(0); ///< Fiber reduced invariant.
  };

  /**
   * @brief First derivatives of passive-energy density with respect to reduced invariants.
   */
  template <class Scalar>
  struct ReducedInvariantGradient
  {
    Scalar dW_dJ1 = Scalar(0); ///< @f$ \partial W / \partial J_1 @f$.
    Scalar dW_dJ2 = Scalar(0); ///< @f$ \partial W / \partial J_2 @f$.
    Scalar dW_dJ4 = Scalar(0); ///< @f$ \partial W / \partial J_4 @f$.
  };

  /**
   * @brief Second derivatives of passive-energy density with respect to reduced invariants.
   */
  template <class Scalar>
  struct ReducedInvariantHessian
  {
    Scalar d2W_dJ1dJ1 = Scalar(0); ///< @f$ \partial^2 W / \partial J_1^2 @f$.
    Scalar d2W_dJ1dJ2 = Scalar(0); ///< @f$ \partial^2 W / \partial J_1 \partial J_2 @f$.
    Scalar d2W_dJ1dJ4 = Scalar(0); ///< @f$ \partial^2 W / \partial J_1 \partial J_4 @f$.
    Scalar d2W_dJ2dJ2 = Scalar(0); ///< @f$ \partial^2 W / \partial J_2^2 @f$.
    Scalar d2W_dJ2dJ4 = Scalar(0); ///< @f$ \partial^2 W / \partial J_2 \partial J_4 @f$.
    Scalar d2W_dJ4dJ4 = Scalar(0); ///< @f$ \partial^2 W / \partial J_4^2 @f$.
  };

  /**
   * @brief Bundle of gradient and Hessian for passive-energy derivatives.
   */
  template <class Scalar>
  struct PassiveEnergyDerivatives
  {
    ReducedInvariantGradient<Scalar> grad; ///< First-order derivatives.
    ReducedInvariantHessian<Scalar> hess;  ///< Second-order derivatives.
  };
}

#endif
