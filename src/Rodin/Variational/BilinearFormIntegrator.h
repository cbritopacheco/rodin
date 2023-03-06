#ifndef RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H
#define RODIN_VARIATIONAL_BILINEARFORMINTEGRATOR_H

#include <set>
#include <memory>
#include <mfem.hpp>

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "Integrator.h"

namespace Rodin::Variational
{
  class BilinearFormIntegratorBase : public Integrator
  {
    public:
      using Parent = Integrator;

      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(const TrialFunction<TrialFES>& u, const TestFunction<TestFES>& v)
        : m_u(u), m_v(v)
      {}

      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(TrialFunction<TrialFES>&& u, const TestFunction<TestFES>& v) = delete;

      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(const TrialFunction<TrialFES>& u, TestFunction<TestFES>&& v) = delete;

      template <class TrialFES, class TestFES>
      BilinearFormIntegratorBase(TrialFunction<TrialFES>&& u, TestFunction<TestFES>&& v) = delete;

      BilinearFormIntegratorBase(const BilinearFormIntegratorBase& other)
        : Parent(other),
          m_u(other.m_u), m_v(other.m_v),
          m_attrs(other.m_attrs)
      {}

      BilinearFormIntegratorBase(BilinearFormIntegratorBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)), m_v(std::move(other.m_v)),
          m_attrs(std::move(other.m_attrs))
      {}

      virtual
      ~BilinearFormIntegratorBase() = default;

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
      BilinearFormIntegratorBase& over(Geometry::Attribute attr)
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
      BilinearFormIntegratorBase& over(const std::set<Geometry::Attribute>& attrs)
      {
        assert(attrs.size() > 0);
        m_attrs = attrs;
        return *this;
      }

      inline
      Integrator::Type getType() const
      final override
      {
        return Integrator::Type::Bilinear;
      }

      /**
       * @brief Gets reference to trial function.
       */
      inline
      const FormLanguage::Base& getTrialFunction() const
      {
        return m_u.get();
      }

      /**
       * @brief Gets reference to test function.
       */
      inline
      const FormLanguage::Base& getTestFunction() const
      {
        return m_v.get();
      }

      /**
       * @brief Performs the assembly of the element matrix for the given
       * element.
       */
      virtual
      Math::Matrix getMatrix(const Geometry::Simplex& element) const = 0;

      virtual
      BilinearFormIntegratorBase* copy() const noexcept override = 0;

    private:
      std::reference_wrapper<const FormLanguage::Base> m_u;
      std::reference_wrapper<const FormLanguage::Base> m_v;
      std::set<Geometry::Attribute> m_attrs;
  };
}

#endif
