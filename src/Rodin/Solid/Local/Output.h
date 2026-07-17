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
 * An Output is the dual of @ref Input: where an Input populates a
 * ConstitutivePoint with auxiliary data before the constitutive law is
 * evaluated, an Output reads the resulting cache at each quadrature point
 * afterwards. This is the extraction site for quantities derived from the
 * evaluation -- committed internal variables, recovered stresses, strain
 * energy, diagnostics -- without hard-coding them into the integrator or the
 * constitutive law, and without recomputing the kinematics at the call site.
 *
 * The cache type is the law's, so an Output is written against a specific
 * law. The ConstitutivePoint carries @ref Tags::CellIndex and
 * @ref Tags::QuadraturePointIndex, so an Output and an Input driven by the
 * same pass address the same quadrature point.
 *
 * ## Usage with CRTP
 *
 * @code
 * struct MyOutput : Solid::Output<MyOutput>
 * {
 *   void extract(const Solid::ConstitutivePoint& cp, const Cache& cache) const
 *   {
 *     // e.g., commit an internal variable
 *     state[cp.get<Solid::Tags::CellIndex>()]
 *          [cp.get<Solid::Tags::QuadraturePointIndex>()] = cache.activeExtension;
 *   }
 * };
 *
 * ivw.setOutput(MyOutput{});
 * @endcode
 *
 * ## Usage with lambda / std::function
 *
 * @code
 * ivw.setOutput([&](const Solid::ConstitutivePoint& cp, const auto& cache) {
 *   energy += cache.strain;
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
   * Cache&)` method that reads quantities out of a constitutive evaluation at
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
       * populated its cache.
       *
       * @param[in] cp The evaluated constitutive point
       * @param[in] cache The law's cache at @p cp
       */
      template <class Cache>
      void extract(const ConstitutivePoint& cp, const Cache& cache) const
      {
        static_cast<const Derived&>(*this).extract(cp, cache);
      }

    protected:
      /// @brief Default constructor for derived CRTP outputs.
      Output() = default;
      /// @brief Copy constructor.
      Output(const Output&) = default;
      /// @brief Move constructor.
      Output(Output&&) = default;
      /// @brief Copy assignment operator.
      Output& operator=(const Output&) = default;
      /// @brief Move assignment operator.
      Output& operator=(Output&&) = default;
  };

  /**
   * @brief Type-erased callable for extraction from a constitutive evaluation.
   * @tparam Cache The constitutive law's cache type.
   */
  template <class Cache>
  using OutputFunction = std::function<void(const ConstitutivePoint&, const Cache&)>;
}

#endif
