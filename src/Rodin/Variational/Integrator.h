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

#include <functional>

#include "Rodin/FormLanguage/Base.h"
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/Types.h"

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
        : Parent(other),
          m_order(other.m_order)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integrator to move
       */
      Integrator(Integrator&& other)
        : Parent(std::move(other)),
          m_order(std::move(other.m_order))
      {}

      /// @brief Virtual destructor
      /**
       * @brief Rule giving the integration order to use on a polytope.
       *
       * @see setOrder(OrderType)
       */
      using OrderType = std::function<size_t(const Geometry::Polytope&)>;

      /**
       * @brief Restores order inference from the integrand.
       *
       * This is the default: the order is derived per polytope from the
       * polynomial degrees of the finite elements living on it, and is
       * therefore right exactly when the integrand's degree is the sum of
       * those degrees --- when nothing else in the expression contributes any.
       *
       * A coefficient does contribute. Where it can report its degree that is
       * accounted for; where it cannot, inference undercounts and the form is
       * under-integrated, and this happens whether or not the integrand is
       * polynomial. A mass form weighted by a lambda returning @f$ x + y @f$
       * is cubic in a first-order basis, but the lambda reports nothing, so
       * inference sees the basis alone and asks for a quadratic rule. The
       * integrand is a polynomial and the answer is still wrong.
       *
       * State the order with @ref setOrder(size_t) whenever the degree the
       * integrand actually has is not the degree inference can see.
       */
      Integrator& setOrder(std::nullopt_t)
      {
        m_order = nullptr;
        return *this;
      }

      /**
       * @brief Sets a constant integration order.
       *
       * The order is the degree the quadrature rule must integrate exactly,
       * so this is how a caller says what its integrand actually needs: the
       * degree of a polynomial integrand that inference cannot see in full,
       * or, for an integrand that is not polynomial at all, the accuracy
       * being asked for.
       *
       * Note that the inferred order is a function of the polytope, so a
       * constant is only equivalent to inference on a mesh whose elements all
       * share one degree; prefer @ref setOrder(OrderType) otherwise.
       *
       * @param[in] order Integration order to use on every polytope
       */
      Integrator& setOrder(size_t order)
      {
        m_order = [order](const Geometry::Polytope&) { return order; };
        return *this;
      }

      /**
       * @brief Sets a rule computing the integration order per polytope.
       *
       * On an integrator the order is the degree of the quadrature rule, not
       * the polynomial degree of an expression. Set it whenever the degree
       * inference can see is not the degree the integrand has: a coefficient
       * that reports no degree undercounts a polynomial integrand, and one
       * that is not polynomial at all --- a level set composed with a mesh
       * lookup, say --- has no degree to infer, so the inferred order is
       * meaningless and typically too low in both cases.
       *
       * @param[in] order Rule invoked with the polytope being integrated
       */
      Integrator& setOrder(OrderType order)
      {
        m_order = std::move(order);
        return *this;
      }

      /**
       * @brief The integration order to use on a polytope.
       * @param[in] polytope Polytope being integrated
       * @returns The order given by the rule, or an empty optional when the
       * order is to be inferred from the integrand.
       */
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const
      {
        if (m_order)
          return m_order(polytope);
        return std::nullopt;
      }

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

    private:
        // Empty unless setOrder() was called, in which case it is the rule
        // resolving the order. Distinct from the m_order each QuadratureRule
        // specialisation declares, which caches the *resolved* order for its
        // recompute check; that member hides this one, which is private here
        // and reached only through getOrder().
      OrderType m_order;
  };
}

#endif
