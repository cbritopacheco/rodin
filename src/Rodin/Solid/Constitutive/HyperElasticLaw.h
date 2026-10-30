/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HyperElasticLaw.h
 * @brief Base class for hyperelastic constitutive laws.
 *
 * Provides the CRTP base class for all hyperelastic constitutive laws.
 * Each derived law must provide:
 * - Non-const `setConstitutivePoint(const ConstitutivePoint&)` to prepare owned data
 * - Const `getFirstPiolaKirchhoffStress(SpatialMatrix&)`
 *   to compute @f$ \mathbf{P} = \partial W / \partial \mathbf{F} @f$
 * - `getMaterialTangent(SpatialMatrix&, const SpatialMatrix&)`
 *   to compute @f$ D\mathbf{P}[\delta\mathbf{F}] @f$
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_HYPERELASTICLAW_H
#define RODIN_SOLID_CONSTITUTIVE_HYPERELASTICLAW_H

#include "Rodin/Types.h"
#include "Rodin/Alert/MemberFunctionException.h"
#include "Rodin/Alert/Raise.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"

namespace Rodin::Solid
{
  /**
   * @brief CRTP base class for hyperelastic constitutive laws.
   *
   * A hyperelastic material is characterized by a stored energy density
   * function @f$ W(\mathbf{F}) @f$ from which all stress measures derive.
   *
   * Each derived class implements:
   * - **Binding**: Captures inputs and prepares private law-specific data.
   * - **First Piola-Kirchhoff stress**: @f$ \mathbf{P} = \frac{\partial W}{\partial \mathbf{F}} @f$
   * - **Material tangent**: Directional derivative
   *   @f$ D\mathbf{P}[\delta\mathbf{F}] = \frac{\partial^2 W}{\partial \mathbf{F}^2} : \delta\mathbf{F} @f$
   *
   * Binding is non-const. Evaluation methods are const and use the captured
   * inputs until the next binding, without retaining references to the point.
   * Copies own independent prepared state; each assembly worker uses its own
   * law. Concurrent binding of the same law requires external synchronization.
   *
   * @tparam Derived The derived constitutive law type (CRTP)
   */
  template <class Derived>
  class HyperElasticLaw
  {
    public:
      /**
       * @brief Binds the law to constitutive inputs through the derived implementation.
       * @param cp Constitutive point whose consumed values are captured by the law.
       * @returns The prepared derived law.
       */
      Derived& setConstitutivePoint(const ConstitutivePoint& cp)
      {
        return static_cast<Derived&>(*this).setConstitutivePoint(cp);
      }

      /**
       * @brief Returns the prepared kinematic snapshot.
       * @returns Owned kinematics of the prepared law.
       */
      const KinematicState& getKinematicState() const
      {
        return static_cast<const Derived&>(*this).getKinematicState();
      }

      /**
       * @brief Computes the prepared first Piola-Kirchhoff stress.
       * @param P Output stress tensor.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P) const
      {
        static_cast<const Derived&>(*this).getFirstPiolaKirchhoffStress(P);
      }

      /**
       * @brief Computes the prepared material tangent action.
       * @param dP Output stress increment.
       * @param dF Deformation-gradient increment.
       */
      void getMaterialTangent(
        Math::SpatialMatrix<Real>& dP, const Math::SpatialMatrix<Real>& dF) const
      {
        static_cast<const Derived&>(*this).getMaterialTangent(dP, dF);
      }

      /**
       * @brief Returns the prepared strain-energy density.
       * @returns Strain-energy density.
       */
      Real getStrainEnergyDensity() const
      {
        return static_cast<const Derived&>(*this).getStrainEnergyDensity();
      }

    protected:
      /**
       * @brief Requires successful preparation before reading constitutive results.
       * @param condition Whether the derived evaluation is prepared.
       * @param function Calling member function name.
       */
      void require(bool condition, const char* function) const
      {
        if (!condition) [[unlikely]]
        {
          raiseUnbound(function);
        }
      }

      /**
       * @brief Compares the independent kinematic inputs exactly.
       * @param lhs Previously captured kinematics.
       * @param rhs Incoming kinematics.
       * @returns Whether the dimensions and displacement gradients agree.
       */
      static bool hasSameKinematics(const KinematicState& lhs, const KinematicState& rhs)
      {
        if (lhs.getDimension() != rhs.getDimension())
        {
          return false;
        }
        const auto& a = lhs.getDisplacementGradient();
        const auto& b = rhs.getDisplacementGradient();
        return (a.getData().topLeftCorner(a.rows(), a.cols()).array() ==
          b.getData().topLeftCorner(b.rows(), b.cols()).array())
          .all();
      }

      /**
       * @brief Compares consumed preferred directions exactly.
       * @param lhs Previously captured direction.
       * @param rhs Incoming direction.
       * @returns Whether the dimensions and coordinates agree.
       */
      static bool hasSameDirection(
        const Math::SpatialVector<Real>& lhs, const Math::SpatialVector<Real>& rhs)
      {
        if (lhs.size() != rhs.size())
        {
          return false;
        }
        for (size_t i = 0; i < lhs.size(); ++i)
        {
          if (lhs(i) != rhs(i))
          {
            return false;
          }
        }
        return true;
      }

      /// @brief Default constructor for derived CRTP laws.
      HyperElasticLaw() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      HyperElasticLaw(const HyperElasticLaw& other) = default;

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      HyperElasticLaw(HyperElasticLaw&& other) = default;

      /**
       * @brief Copy assignment operator.
       * @returns Reference to this object after the operation.
       * @param other Object to copy from.
       */
      HyperElasticLaw& operator=(const HyperElasticLaw& other) = default;

      /**
       * @brief Move assignment operator.
       * @returns Reference to this object after the operation.
       * @param other Object to move from.
       */
      HyperElasticLaw& operator=(HyperElasticLaw&& other) = default;

    private:
      /**
       * @brief Reports a read before successful constitutive binding.
       * @param function Calling member function name.
       */
      void raiseUnbound(const char* function) const
      {
        Alert::MemberFunctionException(static_cast<const Derived&>(*this), function)
          << "Bind a constitutive point before evaluating the material law."
          << Alert::Raise;
      }
  };
}

#endif
