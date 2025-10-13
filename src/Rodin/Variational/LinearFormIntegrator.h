/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_LINEARFORMINTEGRATOR_H

#include "Rodin/Geometry/Region.h"

#include "ForwardDecls.h"
#include "TestFunction.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  template <class Number>
  class LinearFormIntegratorBase : public Integrator
  {
    public:
      using ScalarType = Number;

      using Parent = Integrator;

      class Scatter
      {
        public:
          Scatter() = default;

          Scatter(const Scatter&) = default;

          Scatter(Scatter&&) = default;

          Scatter& operator=(const Scatter&) = default;

          Scatter& operator=(Scatter&&) = default;

          virtual ~Scatter() = default;

          Scatter& push(Index i, ScalarType v)
          {
            m_data.emplace_back(i, v);
            return *this;
          }

          size_t size() const
          {
            return m_data.size();
          }

          const Index& getIndex(size_t i) const
          {
            return m_data[i].first;
          }

          const ScalarType& getDOF(size_t i) const
          {
            return m_data[i].second;
          }

          void clear()
          {
            m_data.clear();
          }

          auto begin() const
          {
            return m_data.begin();
          }

          auto end() const
          {
            return m_data.end();
          }

        private:
          std::vector<std::pair<Index, ScalarType>> m_data;
      };

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
      LinearFormIntegratorBase& over(const Geometry::Attribute& attr)
      {
        return over(FlatSet<Geometry::Attribute>{attr});
      }

      template <class A1, class A2, class ... As>
      LinearFormIntegratorBase& over(const A1& a1, const A2& a2, const As&... attrs)
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

      Scatter& getScatter()
      {
        return m_scatter;
      }

      const Scatter& getScatter() const
      {
        return m_scatter;
      }

      virtual const Geometry::Polytope& getPolytope() const = 0;

      virtual LinearFormIntegratorBase& setPolytope(const Geometry::Polytope& polytope) = 0;

      virtual ScalarType integrate(size_t local) = 0;

      virtual LinearFormIntegratorBase* copy() const noexcept override = 0;

      virtual Geometry::Region getRegion() const = 0;

    private:
      std::unique_ptr<FormLanguage::Base> m_v;
      FlatSet<Geometry::Attribute> m_attrs;
      Scatter m_scatter;
  };
}

#endif
