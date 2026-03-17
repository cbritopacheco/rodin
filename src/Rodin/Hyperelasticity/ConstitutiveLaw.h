/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ConstitutiveLaw.h
 * @brief Abstract base class for hyperelastic constitutive laws.
 *
 * Provides a CRTP interface for hyperelastic material models. Each concrete
 * law must supply:
 * - The first Piola-Kirchhoff stress @f$ \mathbf{P}(\mathbf{F}) @f$
 * - The tangent modulus @f$ \frac{\partial \mathbf{P}}{\partial \mathbf{F}} @f$
 * - The strain energy density @f$ W(\mathbf{F}) @f$
 *
 * The constitutive law is decoupled from mesh and solver internals: it
 * operates solely on the deformation gradient.
 */
#ifndef RODIN_HYPERELASTICITY_CONSTITUTIVELAW_H
#define RODIN_HYPERELASTICITY_CONSTITUTIVELAW_H

#include <Eigen/Dense>

#include "Rodin/Types.h"
#include "Rodin/Copyable.h"

#include "ForwardDecls.h"

namespace Rodin::Hyperelasticity
{
  /**
   * @brief CRTP base for hyperelastic constitutive laws.
   *
   * A constitutive law maps the deformation gradient
   * @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$
   * to:
   * - the first Piola-Kirchhoff stress @f$ \mathbf{P} @f$
   * - the material tangent @f$ \mathbb{A} = \partial \mathbf{P} / \partial \mathbf{F} @f$
   * - the stored energy density @f$ W @f$
   *
   * @tparam Derived Concrete law type (CRTP)
   */
  template <class Derived>
  class ConstitutiveLawBase : public Copyable
  {
    public:
      virtual ~ConstitutiveLawBase() = default;

      /**
       * @brief Computes the first Piola-Kirchhoff stress.
       * @param[in]  F  Deformation gradient (d x d)
       * @param[out] P  First Piola-Kirchhoff stress (d x d)
       */
      void stress(const Eigen::Ref<const Eigen::MatrixXd>& F,
                  Eigen::Ref<Eigen::MatrixXd> P) const
      {
        static_cast<const Derived*>(this)->stressImpl(F, P);
      }

      /**
       * @brief Computes the material tangent modulus.
       *
       * The tangent is stored as a @f$ d^2 \times d^2 @f$ matrix with the
       * Voigt-like indexing @f$ C(i d + J,\, k d + L) =
       * \partial P_{iJ} / \partial F_{kL} @f$.
       *
       * @param[in]  F  Deformation gradient (d x d)
       * @param[out] C  Tangent modulus (d*d x d*d)
       */
      void tangent(const Eigen::Ref<const Eigen::MatrixXd>& F,
                   Eigen::Ref<Eigen::MatrixXd> C) const
      {
        static_cast<const Derived*>(this)->tangentImpl(F, C);
      }

      /**
       * @brief Computes the strain energy density.
       * @param[in] F Deformation gradient (d x d)
       * @returns @f$ W(\mathbf{F}) @f$
       */
      Real energy(const Eigen::Ref<const Eigen::MatrixXd>& F) const
      {
        return static_cast<const Derived*>(this)->energyImpl(F);
      }
  };
}

#endif
