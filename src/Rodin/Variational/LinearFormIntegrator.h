/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Integrator.h"
#include "TestFunction.h"

namespace Rodin::Variational
{
  class LinearFormIntegratorBase : public Integrator
  {
    public:
      using Parent = Integrator;

      LinearFormIntegratorBase(const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& v)
        : m_v(v.copy())
      {}

      LinearFormIntegratorBase(const LinearFormIntegratorBase& other)
        : Parent(other),
          m_v(other.m_v->copy()),
          m_attrs(other.m_attrs)
      {}

      LinearFormIntegratorBase(LinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      virtual ~LinearFormIntegratorBase() = default;

      const ShapeFunctionBase<ShapeFunctionSpaceType::Test>& getTestFunction() const
      {
        assert(m_v);
        return *m_v;
      }

      /**
       * @brief Gets the attributes of the elements being integrated.
       */
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
      LinearFormIntegratorBase& over(const std::set<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      /**
       * @internal
       */
      std::unique_ptr<mfem::LinearFormIntegrator> build() const;

      Integrator::Type getType() const override
      {
        return Integrator::Type::Linear;
      }

      /**
       * @brief Performs the assembly of the element vector for the given
       * element.
       */
      virtual mfem::Vector getVector(
          const Geometry::Simplex& element) const = 0;

      virtual LinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      std::unique_ptr<ShapeFunctionBase<ShapeFunctionSpaceType::Test>> m_v;
      std::set<Geometry::Attribute> m_attrs;
  };
}

#endif
