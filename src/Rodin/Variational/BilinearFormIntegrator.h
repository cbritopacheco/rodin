/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BilinearFormIntegrator.h
 * @brief Base classes for bilinear form integrators.
 *
 * This file defines the abstract base classes for bilinear form integrators,
 * which are responsible for computing local contributions to the global system
 * matrix in finite element assembly. Integrators represent specific terms in
 * the weak formulation involving both trial and test functions.
 */
#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <memory>

#include "Rodin/Geometry/Types.h"
#include "Rodin/Geometry/Region.h"
#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Abstract base class for bilinear form integrators.
   *
   * BilinearFormIntegratorBase provides the foundation for implementing specific
   * bilinear form terms in variational formulations. A bilinear form integrator
   * computes the local element matrix contributions:
   * @f[
   *   A^K_{ij} = \int_K a(\phi_j, \psi_i) \, dx
   * @f]
   * where @f$ \phi_j @f$ are trial basis functions and @f$ \psi_i @f$ are test
   * basis functions on element @f$ K @f$.
   *
   * ## Role in Assembly
   * During finite element assembly:
   * 1. The integrator is called for each element in the mesh
   * 2. Local element matrices are computed using quadrature
   * 3. Local matrices are assembled into the global system matrix
   *
   * ## Common Bilinear Forms
   * - **Mass matrix**: @f$ \int_K u \cdot v \, dx @f$
   * - **Stiffness matrix**: @f$ \int_K \nabla u \cdot \nabla v \, dx @f$
   * - **Advection**: @f$ \int_K (\beta \cdot \nabla u) v \, dx @f$
   * - **Reaction**: @f$ \int_K c \, u \cdot v \, dx @f$
   *
   * @tparam Number Scalar type for matrix entries
   * @tparam Derived Derived integrator class (CRTP pattern)
   *
   * @see LinearFormIntegratorBase, BilinearForm
   */
  template <class Number, class Derived>
  class BilinearFormIntegratorBase : public Integrator
  {
    public:
      /// @brief Scalar type for matrix entries
      using ScalarType = Number;

      /// @brief Parent class type
      using Parent = Integrator;

      /**
       * @brief Constructs a bilinear form integrator.
       * @param[in] u Trial function
       * @param[in] v Test function
       *
       * Creates an integrator that will compute contributions involving both
       * the trial function @f$ u @f$ and test function @f$ v @f$.
       */
      template <class Solution, class TrialFES, class TestFES>
      BilinearFormIntegratorBase(
          const TrialFunction<Solution, TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u.copy()), m_v(v.copy())
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Integrator to copy
       */
      BilinearFormIntegratorBase(const BilinearFormIntegratorBase& other)
        : Parent(other),
          m_u(other.m_u->copy()), m_v(other.m_v->copy())
      {}

      /**
       * @brief Copy constructor from different scalar type.
       * @param[in] other Integrator to copy
       */
      template <class OtherNumber, class OtherDerived>
      BilinearFormIntegratorBase(const BilinearFormIntegratorBase<OtherNumber, OtherDerived>& other)
        : Parent(other),
          m_u(other.m_u->copy()), m_v(other.m_v->copy())
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Integrator to move
       */
      BilinearFormIntegratorBase(BilinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v))
      {}

      /**
       * @brief Move constructor from different scalar type.
       * @param[in] other Integrator to move
       */
      template <class OtherNumber, class OtherDerived>
      BilinearFormIntegratorBase(BilinearFormIntegratorBase<OtherNumber, OtherDerived>&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v))
      {}

      /// @brief Virtual destructor
      virtual
      ~BilinearFormIntegratorBase() = default;

      /**
       * @brief Gets the integrator type.
       * @returns Integrator::Type::Bilinear
       */
      Integrator::Type getType() const final override
      {
        return Integrator::Type::Bilinear;
      }

      /**
       * @brief Gets the trial function.
       * @returns Const reference to the trial function
       */
      const FormLanguage::Base& getTrialFunction() const
      {
        assert(m_u);
        return *m_u;
      }

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
       * @brief Creates a copy of this integrator.
       * @returns Pointer to newly allocated copy
       */
      virtual
      BilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      std::unique_ptr<FormLanguage::Base> m_u;
      std::unique_ptr<FormLanguage::Base> m_v;
  };

  template <class Number>
  class LocalBilinearFormIntegratorBase
    : public BilinearFormIntegratorBase<Number, LocalBilinearFormIntegratorBase<Number>>
  {
    public:
      using ScalarType = Number;

      using Parent = BilinearFormIntegratorBase<ScalarType, LocalBilinearFormIntegratorBase>;

      using Parent::Parent;

      /**
       * @brief Copy constructor.
       */
      template <class OtherNumber>
      LocalBilinearFormIntegratorBase(const LocalBilinearFormIntegratorBase<OtherNumber>& other)
        : Parent(other),
          m_attrs(other.m_attrs)
      {}

      /**
       * @brief Move constructor.
       */
      template <class OtherNumber>
      LocalBilinearFormIntegratorBase(LocalBilinearFormIntegratorBase<OtherNumber>&& other)
        : Parent(std::move(other)),
          m_attrs(std::move(other.m_attrs))
      {}

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
      const FlatSet<Geometry::Attribute>& getAttributes() const
      {
        return m_attrs;
      }

      /**
       * @brief Specifies the material reference over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material reference over which the integration should
       * take place.
       */
      LocalBilinearFormIntegratorBase& over(const Geometry::Attribute& attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      template <class A1, class A2, class ... As>
      LocalBilinearFormIntegratorBase& over(const A1& a1, const A2& a2, const As&... attrs)
      {
        return over(FlatSet<Geometry::Attribute>{a1, a2, attrs...});
      }

      /**
       * @brief Specifies the material references over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material references over which the integration should
       * take place.
       */
      LocalBilinearFormIntegratorBase& over(const FlatSet<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      virtual const Geometry::Polytope& getPolytope() const = 0;

      virtual LocalBilinearFormIntegratorBase& setPolytope(const Geometry::Polytope& polytope) = 0;

      virtual ScalarType integrate(size_t tr, size_t te) = 0;

      virtual Geometry::Region getRegion() const = 0;

      virtual
      LocalBilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      FlatSet<Geometry::Attribute> m_attrs;
  };

  template <class Number>
  class GlobalBilinearFormIntegratorBase
    : public BilinearFormIntegratorBase<Number, GlobalBilinearFormIntegratorBase<Number>>
  {
    public:
      using ScalarType = Number;

      using Parent = BilinearFormIntegratorBase<ScalarType, GlobalBilinearFormIntegratorBase<ScalarType>>;

      using Parent::Parent;

      template <class OtherNumber>
      GlobalBilinearFormIntegratorBase(const GlobalBilinearFormIntegratorBase<OtherNumber>& other)
        : Parent(other),
          m_trialAttrs(other.m_trialAttrs),
          m_testAttrs(other.m_testAttrs)
      {}

      /**
       * @brief Move constructor.
       */
      template <class OtherNumber>
      GlobalBilinearFormIntegratorBase(GlobalBilinearFormIntegratorBase<OtherNumber>&& other)
        : Parent(std::move(other)),
          m_trialAttrs(std::move(other.m_trialAttrs)),
          m_testAttrs(std::move(other.m_testAttrs))
      {}

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
      const FlatSet<Geometry::Attribute>& getTrialAttributes() const
      {
        return m_trialAttrs;
      }

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
      const FlatSet<Geometry::Attribute>& getTestAttributes() const
      {
        return m_testAttrs;
      }

      /**
       * @brief Specifies the material reference over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material reference over which the integration should
       * take place.
       */
      GlobalBilinearFormIntegratorBase& setTrialAttributes(
          const FlatSet<Geometry::Attribute>& attrs)
      {
        m_trialAttrs = attrs;
        return *this;
      }

      /**
       * @brief Specifies the material reference over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material reference over which the integration should
       * take place.
       */
      GlobalBilinearFormIntegratorBase& setTestAttributes(
          const FlatSet<Geometry::Attribute>& attrs)
      {
        m_testAttrs = attrs;
        return *this;
      }

      virtual
      GlobalBilinearFormIntegratorBase& setPolytope(const Geometry::Polytope& tau, const Geometry::Polytope& t) = 0;

      virtual ScalarType integrate(size_t tr, size_t te) = 0;

      virtual Geometry::Region getTrialRegion() const = 0;

      virtual Geometry::Region getTestRegion() const = 0;

      virtual
      GlobalBilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      FlatSet<Geometry::Attribute> m_trialAttrs;
      FlatSet<Geometry::Attribute> m_testAttrs;
  };
}

#endif
