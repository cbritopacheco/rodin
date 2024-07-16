/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>

#include "Rodin/Cast.h"
#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"
#include "Integrator.h"
#include "TestFunction.h"

namespace Rodin::Variational
{
  template <class Number>
  class LinearFormIntegratorBase : public Integrator
  {
    public:
      using ScalarType = Number;

      using Parent = Integrator;

      template <class FES>
      LinearFormIntegratorBase(const TestFunction<FES>& v)
        : m_v(v.copy())
      {}

      LinearFormIntegratorBase(const LinearFormIntegratorBase& other)
        : Parent(other),
          m_v(other.m_v->copy()),
          m_attrs(other.m_attrs)
      {}

      template <class OtherNumber>
      LinearFormIntegratorBase(const LinearFormIntegratorBase<OtherNumber>& other)
        : Parent(other),
          m_v(other.m_v->copy()),
          m_attrs(other.m_attrs)
      {}

      LinearFormIntegratorBase(LinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      template <class OtherNumber>
      LinearFormIntegratorBase(LinearFormIntegratorBase<OtherNumber>&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      virtual ~LinearFormIntegratorBase() = default;

      const FormLanguage::Base& getTestFunction() const
      {
        assert(m_v);
        return *m_v;
      }

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
      LinearFormIntegratorBase& over(Geometry::Attribute attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      template <class A1, class A2, class ... As>
      LinearFormIntegratorBase& over(A1 a1, A2 a2, As... attrs)
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

      virtual
      LinearFormIntegratorBase* copy() const noexcept override = 0;

      virtual Region getRegion() const = 0;

    private:
      std::unique_ptr<FormLanguage::Base> m_v;
      FlatSet<Geometry::Attribute> m_attrs;
  };
}

#endif
