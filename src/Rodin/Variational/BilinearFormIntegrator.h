#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <memory>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  /**
   * @brief Abstract base class for bilinear form integrators.
   *
   * This class provides the base functionality for bilinear form integrator
   * objects.
   */
  template <class Number, class Derived>
  class BilinearFormIntegratorBase : public Integrator
  {
    public:
      using ScalarType = Number;

      using Parent = Integrator;

      /**
       * @brief Constructs the object given a TrialFunction and a TestFunction.
       */
      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u.copy()), m_v(v.copy())
      {}

      /**
       * @brief Copy constructor.
       */
      BilinearFormIntegratorBase(const BilinearFormIntegratorBase& other)
        : Parent(other),
          m_u(other.m_u->copy()), m_v(other.m_v->copy())
      {}

      /**
       * @brief Copy constructor.
       */
      template <class OtherNumber, class OtherDerived>
      BilinearFormIntegratorBase(const BilinearFormIntegratorBase<OtherNumber, OtherDerived>& other)
        : Parent(other),
          m_u(other.m_u->copy()), m_v(other.m_v->copy())
      {}

      /**
       * @brief Move constructor.
       */
      BilinearFormIntegratorBase(BilinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v))
      {}

      /**
       * @brief Move constructor.
       */
      template <class OtherNumber, class OtherDerived>
      BilinearFormIntegratorBase(BilinearFormIntegratorBase<OtherNumber, OtherDerived>&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v))
      {}

      virtual
      ~BilinearFormIntegratorBase() = default;

      Integrator::Type getType() const final override
      {
        return Integrator::Type::Bilinear;
      }

      /**
       * @brief Gets a constant reference to trial function object.
       */
      const FormLanguage::Base& getTrialFunction() const
      {
        assert(m_u);
        return *m_u;
      }

      /**
       * @brief Gets a constant reference to test function object.
       */
      const FormLanguage::Base& getTestFunction() const
      {
        assert(m_v);
        return *m_v;
      }

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
      LocalBilinearFormIntegratorBase& over(Geometry::Attribute attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      template <class A1, class A2, class ... As>
      LocalBilinearFormIntegratorBase& over(A1 a1, A2 a2, As... attrs)
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

      virtual Integrator::Region getRegion() const = 0;

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

      virtual Integrator::Region getTrialRegion() const = 0;

      virtual Integrator::Region getTestRegion() const = 0;

      virtual
      GlobalBilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      FlatSet<Geometry::Attribute> m_trialAttrs;
      FlatSet<Geometry::Attribute> m_testAttrs;
  };
}

#endif
