/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>

#include "Rodin/Math/Vector.h"

#include "ForwardDecls.h"
#include "Integrator.h"
#include "TestFunction.h"

namespace Rodin::Variational
{
  class LinearFormIntegratorBase : public Integrator
  {
    public:
      using Parent = Integrator;

      template <class FES>
      LinearFormIntegratorBase(const TestFunction<FES>& v)
        : m_v(v)
      {}

      template <class FES>
      LinearFormIntegratorBase(TestFunction<FES>&& v) = delete;

      LinearFormIntegratorBase(const LinearFormIntegratorBase& other)
        : Parent(other),
          m_v(other.m_v),
          m_attrs(other.m_attrs)
      {}

      LinearFormIntegratorBase(LinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      virtual ~LinearFormIntegratorBase() = default;

      inline
      const FormLanguage::Base& getTestFunction() const
      {
        return m_v.get();
      }

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
      inline
      const std::set<Geometry::Attribute>& getAttributes() const
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
      inline
      LinearFormIntegratorBase& over(Geometry::Attribute attr)
      {
        return over(std::set<Geometry::Attribute>{attr});
      }

      /**
       * @brief Specifies the material references over which to integrate.
       * @returns Reference to self (for method chaining)
       *
       * Specifies the material references over which the integration should
       * take place.
       */
      inline
      LinearFormIntegratorBase& over(const std::set<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      inline
      Integrator::Type getType() const
      final override
      {
        return Integrator::Type::Linear;
      }

      inline
      const Math::Vector getVector() const
      {
        return m_vector;
      }

      inline
      Math::Vector& getVector()
      {
        return m_vector;
      }

      /**
       * @brief Performs the assembly of the element vector for the given
       * element.
       */
      virtual
      void assemble(const Geometry::Polytope& element) = 0;

      virtual
      LinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<const FormLanguage::Base> m_v;
      std::set<Geometry::Attribute> m_attrs;
      Math::Vector m_vector;
  };
}

#endif
