/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Integrator.h
 * @brief Base class for all integrator types in variational formulations.
 *
 * This file defines the Integrator base class, which serves as the common
 * ancestor for both linear and bilinear form integrators. Integrators are
 * responsible for computing local contributions during finite element assembly.
 */
#ifndef RODIN_VARIATIONAL_INTEGRATOR_H
#define RODIN_VARIATIONAL_INTEGRATOR_H

#include "Rodin/FormLanguage/Base.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Abstract base class for integrators in variational formulations.
   *
   * The Integrator class provides a common interface for all types of
   * integrators used in finite element assembly. Integrators compute local
   * (element-level) contributions that are then assembled into global matrices
   * and vectors.
   *
   * ## Integrator Types
   * There are two main types of integrators:
   * - **Linear Integrators**: Compute contributions to the load vector (RHS)
   *   - Example: @f$ b_i = \int_K f \psi_i \, dx @f$
   * - **Bilinear Integrators**: Compute contributions to the system matrix
   *   - Example: @f$ A_{ij} = \int_K \nabla \phi_j \cdot \nabla \psi_i \, dx @f$
   *
   * ## Role in Assembly
   * During the assembly process:
   * 1. Integrators are iterated over mesh elements
   * 2. Local contributions are computed using numerical quadrature
   * 3. Local contributions are mapped to global indices
   * 4. Global matrix/vector is updated with local contributions
   *
   * @see LinearFormIntegratorBase, BilinearFormIntegratorBase
   */
  class Integrator : public FormLanguage::Base
  {
    public:
      /// @brief Parent class type
      using Parent = FormLanguage::Base;

      /**
       * @brief Enumeration of integrator types.
       */
      enum class Type
      {
        Linear,   ///< Linear form integrator (load vector)
        Bilinear  ///< Bilinear form integrator (system matrix)
      };

      /// @brief Default constructor
      Integrator() = default;

      /**
       * @brief Copy constructor.
       * @param[in] other Integrator to copy
       */
      Integrator(const Integrator& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integrator to move
       */
      Integrator(Integrator&& other)
        : Parent(std::move(other))
      {}

      /// @brief Virtual destructor
      virtual ~Integrator() = default;

      /**
       * @brief Gets the type of this integrator.
       * @returns Type indicating whether this is a linear or bilinear integrator
       */
      virtual Type getType() const = 0;

      /**
       * @brief Creates a copy of this integrator.
       * @returns Pointer to newly allocated copy
       */
      virtual Integrator* copy() const noexcept override = 0;
  };
}

#endif
