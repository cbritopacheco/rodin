/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearFormIntegrator.h
 * @brief Base classes for linear form integrators.
 *
 * This file defines the abstract base classes for linear form integrators,
 * which are responsible for computing local contributions to the right-hand
 * side vector (load vector) in finite element assembly. Integrators represent
 * specific terms in the weak formulation involving only test functions.
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include "Rodin/Geometry/Region.h"

#include "ForwardDecls.h"
#include "TestFunction.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Abstract base class for linear form integrators.
   *
   * LinearFormIntegratorBase provides the foundation for implementing specific
   * linear form terms in variational formulations. A linear form integrator
   * computes the local element vector contributions:
   * @f[
   *   b^K_i = \int_K l(\psi_i) \, dx
   * @f]
   * where @f$ \psi_i @f$ are test basis functions on element @f$ K @f$ and
   * @f$ l(\cdot) @f$ is a linear functional.
   *
   * ## Role in Assembly
   * During finite element assembly:
   * 1. The integrator is called for each element (or boundary face) in the mesh
   * 2. Local element vectors are computed using quadrature
   * 3. Local vectors are assembled into the global right-hand side vector
   *
   * ## Common Linear Forms
   * - **Volume load**: @f$ \int_K f \cdot v \, dx @f$
   * - **Neumann BC**: @f$ \int_{\partial K} g \cdot v \, ds @f$
   * - **Point source**: @f$ f_i v(x_i) @f$
   * - **Initial conditions**: @f$ \int_K u_0 \cdot v \, dx @f$
   *
   * ## Domain Specification
   * Linear form integrators support the `over()` method to restrict integration
   * to specific mesh regions identified by attributes:
   * ```cpp
   * auto load = Integral(f, v).over(1, 2);  // Integrate over regions 1 and 2
   * ```
   *
   * @tparam Number Scalar type for vector entries
   *
   * @see BilinearFormIntegratorBase, LinearForm
   */
  template <class Number>
  class LinearFormIntegratorBase : public Integrator
  {
    public:
      /// @brief Scalar type for vector entries
      using ScalarType = Number;

      /// @brief Parent class type
      using Parent = Integrator;

      /**
       * @brief Constructs a linear form integrator.
       * @param[in] v Test function
       *
       * Creates an integrator that will compute contributions involving
       * the test function @f$ v @f$.
       */
      template <class FES>
      LinearFormIntegratorBase(const TestFunction<FES>& v)
        : m_v(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integrator to copy
       */
      LinearFormIntegratorBase(const LinearFormIntegratorBase& other)
        : Parent(other),
          m_v(other.m_v->copy()),
          m_attrs(other.m_attrs)
      {}

      /**
       * @brief Copy constructor from different scalar type.
       * @param[in] other Integrator to copy
       */
      template <class OtherNumber>
      LinearFormIntegratorBase(const LinearFormIntegratorBase<OtherNumber>& other)
        : Parent(other),
          m_v(other.m_v->copy()),
          m_attrs(other.m_attrs)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integrator to move
       */
      LinearFormIntegratorBase(LinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      /**
       * @brief Move constructor from different scalar type.
       * @param[in] other Integrator to move
       */
      template <class OtherNumber>
      LinearFormIntegratorBase(LinearFormIntegratorBase<OtherNumber>&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      /// @brief Virtual destructor
      virtual ~LinearFormIntegratorBase() = default;

      /**
       * @brief Gets the test function.
       * @returns Const reference to the test function
       */
      const FormLanguage::Base& getTestFunction() const
      {
        assert(m_v);
        return *m_v;
      }

      /**
       * @brief Gets the attributes of the integration domain.
       * @returns Set of mesh attributes where integration is performed
       *
       * Returns the mesh attributes (region markers) over which this integrator
       * will perform integration. An empty set means integrate over all regions.
       */
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_attrs;
      }

      /**
       * @brief Restricts integration to a specific mesh region.
       * @param[in] attr Mesh attribute (region marker)
       * @returns Reference to this integrator (for method chaining)
       *
       * Specifies that integration should only occur over elements/faces with
       * the given attribute. Multiple calls to over() are cumulative.
       */
      LinearFormIntegratorBase& over(const Geometry::Attribute& attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      /**
       * @brief Restricts integration to multiple mesh regions.
       * @param[in] a1 First attribute
       * @param[in] a2 Second attribute
       * @param[in] attrs Additional attributes
       * @returns Reference to this integrator (for method chaining)
       */
      template <class A1, class A2, class ... As>
      LinearFormIntegratorBase& over(const A1& a1, const A2& a2, const As&... attrs)
      {
        return over(FlatSet<Geometry::Attribute>{a1, a2, attrs...});
      }

      /**
       * @brief Restricts integration to a set of mesh regions.
       * @param[in] attrs Set of mesh attributes
       * @returns Reference to this integrator (for method chaining)
       *
       * Specifies the mesh attributes over which the integration should take
       * place. This replaces any previously specified attributes.
       */
      LinearFormIntegratorBase& over(const FlatSet<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      Integrator::Type getType() const final override
      {
        return Integrator::Type::Linear;
      }

      virtual const Geometry::Polytope& getPolytope() const = 0;

      virtual LinearFormIntegratorBase& setPolytope(const Geometry::Polytope& polytope) = 0;

      virtual ScalarType integrate(size_t local) = 0;

      virtual LinearFormIntegratorBase* copy() const noexcept override = 0;

      virtual Geometry::Region getRegion() const = 0;

    private:
      std::unique_ptr<FormLanguage::Base> m_v;
      FlatSet<Geometry::Attribute> m_attrs;
  };
}

#endif
