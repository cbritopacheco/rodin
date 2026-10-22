/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Output.h
 * @brief Extensible extraction mechanism for constitutive evaluation.
 *
 * An Output is the dual of Input: where an Input populates a
 * ConstitutivePoint with auxiliary data before the constitutive law is
 * evaluated, an Output reads the prepared law at each quadrature point
 * afterwards. This is the extraction site for quantities derived from the
 * evaluation -- committed internal variables, recovered stresses, strain
 * energy, diagnostics -- without hard-coding them into the integrator or the
 * constitutive law, and without recomputing the kinematics at the call site.
 *
 * An Output receives a const prepared law, so it is written against a specific
 * law. The ConstitutivePoint carries @c Tags::CellIndex and
 * @c Tags::QuadraturePointIndex, so an Output and an Input driven by the
 * same pass address the same quadrature point.
 *
 * ## Usage with CRTP
 *
 * @code
 * struct MyOutput : Solid::Output<MyOutput>
 * {
 *   MyOutput() = default;
 *
 *   template <class Law>
 *   void extract(const Solid::ConstitutivePoint& cp, const Law& law) const
 *   {
 *     // e.g., commit an internal variable
 *     state[cp.get<Solid::Tags::CellIndex>()]
 *          [cp.get<Solid::Tags::QuadraturePointIndex>()] = law.getActiveExtension();
 *   }
 * };
 *
 * ivw.setOutput([output = MyOutput{}](const Solid::ConstitutivePoint& cp,
 *                                  const auto& law) {
 *   const Solid::Output<MyOutput>& base = output;
 *   base.extract(cp, law);
 * });
 * @endcode
 *
 * ## Usage with lambda / std::function
 *
 * @code
 * ivw.setOutput([&](const Solid::ConstitutivePoint& cp, const auto& law) {
 *   energy += law.getStrainEnergyDensity();
 * });
 * ivw.commit();
 * @endcode
 */
#ifndef RODIN_SOLID_LOCAL_OUTPUT_H
#define RODIN_SOLID_LOCAL_OUTPUT_H

#include <functional>

#include "ConstitutivePoint.h"

namespace Rodin::Solid
{
  /**
   * @brief CRTP base class for outputs.
   *
   * Derived classes implement an `extract(const ConstitutivePoint&, const
   * Law&)` method that reads quantities out of a constitutive evaluation at
   * each quadrature point.
   *
   * @tparam Derived The concrete output type (CRTP)
   */
  template <class Derived>
  class Output
  {
    public:
      /**
       * @brief Extracts data from a constitutive evaluation.
       *
       * Called at each quadrature point after the constitutive law has
       * prepared its owned state.
       *
       * @param[in] cp The evaluated constitutive point
       * @param[in] law The prepared law at @p cp
       */
      template <class Law>
      void extract(const ConstitutivePoint& cp, const Law& law) const
      {
        static_cast<const Derived&>(*this).extract(cp, law);
      }

    protected:
      /// @brief Default constructor for derived CRTP outputs.
      Output() = default;

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Output(const Output& other) = default;

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Output(Output&& other) = default;

      /**
       * @brief Copy assignment operator.
       * @returns Reference to this object after the operation.
       * @param other Object to copy from.
       */
      Output& operator=(const Output& other) = default;

      /**
       * @brief Move assignment operator.
       * @returns Reference to this object after the operation.
       * @param other Object to move from.
       */
      Output& operator=(Output&& other) = default;
  };

  /**
   * @brief Type-erased callable for extraction from a constitutive evaluation.
   * @tparam Law The prepared constitutive law type.
   */
  template <class Law>
  using OutputFunction = std::function<void(const ConstitutivePoint&, const Law&)>;
}

#endif
